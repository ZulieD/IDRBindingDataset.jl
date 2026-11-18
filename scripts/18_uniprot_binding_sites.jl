#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using DataFrames
using IDRBindingDataset:
                         UNIPROT_DEFAULT_BATCH_SIZE,
                         _UNIPROT_BINDING_SITES_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_uniprot_binding_sites_table,
                         fetch_uniprot_entries,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "derived", "uniprot_binding_sites"))
const DEFAULT_OUTPUT_FILENAME = "uniprot_binding_sites.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the UniProt binding-site script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Fetch UniProtKB entries of the dataset proteins (REST API, batched requests)
and keep their `Binding site` features that overlap a DisProt disordered
residue (label `1` or `2`). Writes one row per DisProt protein, keyed by
`protein_disprot_id` and `protein_id`.

Run from repository root:
  julia scripts/18_uniprot_binding_sites.jl

Complete example:
  julia scripts/18_uniprot_binding_sites.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --batch_size 90 --output_dir data/derived/uniprot_binding_sites --output_filename uniprot_binding_sites.csv.gz

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip) with `protein_id` and `residue_labels`."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--batch_size"
        help = "UniProt accessions per request."
        arg_type = Int
        default = UNIPROT_DEFAULT_BATCH_SIZE
        "--output_dir"
        help = "Output directory for the CSV and its README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
        "--output_filename"
        help = "Output compressed CSV filename written under `output_dir`."
        arg_type = String
        default = DEFAULT_OUTPUT_FILENAME
    end

    parsed = parse_args(args, settings)
    (
        dataset_path = normpath(parsed["dataset_path"]),
        batch_size = Int(parsed["batch_size"]),
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"])
    )
end

function main()
    options = parse_command_line(ARGS)
    output_path = joinpath(options.output_dir, options.output_filename)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files([output_path, readme_path])

    dataset = read_final_dataset_csv(options.dataset_path)
    fetched = fetch_uniprot_entries(String.(dataset.protein_id);
        batch_size = options.batch_size)
    table = build_uniprot_binding_sites_table(dataset, fetched.entries)

    _write_csv_gz(table, output_path)
    n_found = count(table.uniprot_entry_found)
    n_binding = count(table.binding_in_uniprot)
    _write_artifact_readme(readme_path;
        title = "UniProt Binding Sites On Disordered Residues",
        description = """
UniProtKB `Binding site` features of the dataset proteins, fetched from
`https://rest.uniprot.org/uniprotkb/search`. Only sites overlapping at least one
DisProt residue labelled `1` (disorder) or `2` (binding) are kept.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "18_uniprot_binding_sites.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "UniProt release" => coalesce(fetched.release, "not reported")],
        outputs = [options.output_filename => "$(nrow(table)) proteins",
            "Proteins found in UniProt" => n_found,
            "Proteins with binding sites on disordered residues" => n_binding],
        parameters = ["Accessions per request" => options.batch_size],
        columns = _UNIPROT_BINDING_SITES_COLUMN_DEFINITIONS)

    println("UniProt binding sites written to $(output_path) (proteins=$(nrow(table)), found=$(n_found), with_binding=$(n_binding), release=$(coalesce(fetched.release, "unknown"))).")
    nothing
end

main()
