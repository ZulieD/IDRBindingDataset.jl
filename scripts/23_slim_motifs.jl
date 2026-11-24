#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
                         SLIM_DEFAULT_ELM_TYPES,
                         _SLIM_LABELS_COLUMN_DEFINITIONS,
                         _parse_excluded_instance_logic,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_slim_labels_table,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_INSTANCES_PATH = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "elm", "instances.tsv"))
const DEFAULT_ELMS_INDEX_PATH = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "elm", "elms_index.tsv"))
# ELM-derived data stays in scratch (not versioned), like 13_mapping_elm_motifs.jl.
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "elm", "slim_labels"))
const DEFAULT_OUTPUT_FILENAME = "slim_labels.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the SLiM label script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Short linear motif (SLiM) evidence from ELM for each dataset protein, keyed by
`protein_disprot_id` and `protein_id`. A residue gets `Label_SLiM = 1` when it
lies in an annotated ELM instance of the protein, or in a match of an ELM class
regex lying entirely in a DisProt binding region (labels 2); `NaN` elsewhere.
No disorder information. Used as the fifth source of 24_final_labels.jl.

Input: ELM `instances.tsv` and `elms_index.tsv` downloaded by
11_download_elm.jl. The output stays in data/scratch/elm (ELM-derived data is
not versioned).

Run from repository root:
  julia scripts/23_slim_motifs.jl

Complete example:
  julia scripts/23_slim_motifs.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --instances_path data/scratch/elm/instances.tsv --elms_index_path data/scratch/elm/elms_index.tsv --elm_types LIG,DOC,DEG --output_dir data/scratch/elm/slim_labels

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip)."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--instances_path"
        help = "ELM instances TSV (11_download_elm.jl)."
        arg_type = String
        default = DEFAULT_INSTANCES_PATH
        "--elms_index_path"
        help = "ELM classes index TSV with the regexes (11_download_elm.jl)."
        arg_type = String
        default = DEFAULT_ELMS_INDEX_PATH
        "--elm_types"
        help = "ELM class types used, comma-separated among LIG, DOC, DEG, MOD, CLV, TRG."
        arg_type = String
        default = join(SLIM_DEFAULT_ELM_TYPES, ",")
        "--max_probability"
        help = "Skip ELM classes whose regex chance-match probability (elms_index `Probability`) is above this value (1: keep all)."
        arg_type = Float64
        default = 1.0
        "--exclude_instance_logic"
        help = "Comma-separated InstanceLogic values excluded from the instances."
        arg_type = String
        default = join(DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC, ",")
        "--output_dir"
        help = "Output directory for the CSV and its README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
    end

    parsed = parse_args(args, settings)
    (
        dataset_path = normpath(parsed["dataset_path"]),
        instances_path = normpath(parsed["instances_path"]),
        elms_index_path = normpath(parsed["elms_index_path"]),
        elm_types = Tuple(uppercase.(strip.(split(parsed["elm_types"], ",";
            keepempty = false)))),
        max_probability = parsed["max_probability"],
        excluded_instance_logic = _parse_excluded_instance_logic(
            parsed["exclude_instance_logic"]),
        output_dir = normpath(parsed["output_dir"])
    )
end

"Read an ELM TSV (comment lines start with `#`) with all columns as strings."
function read_elm_tsv(path::AbstractString)
    DataFrame(CSV.File(path; delim = '\t',
        comment = "#", types = String))
end

function main()
    options = parse_command_line(ARGS)
    output_path = joinpath(options.output_dir, DEFAULT_OUTPUT_FILENAME)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files([output_path, readme_path])

    dataset = read_final_dataset_csv(options.dataset_path)
    instances = read_elm_tsv(options.instances_path)
    elms_index = read_elm_tsv(options.elms_index_path)

    table = build_slim_labels_table(dataset, instances, elms_index;
        elm_types = options.elm_types, max_probability = options.max_probability,
        excluded_instance_logic = options.excluded_instance_logic)
    _write_csv_gz(table, output_path)

    with_instances = count(>(0), table.n_instances)
    with_matches = count(>(0), table.n_regex_matches)
    with_slim = count(!ismissing, table.slim_positions)
    _write_artifact_readme(readme_path;
        title = "Short Linear Motif (SLiM) Residue Labels From ELM",
        description = """
One row per dataset protein. `Label_SLiM` has one value per residue of the
UniProt sequence: `1` on residues of an annotated ELM instance of the protein
(`InstanceLogic` not excluded) or of a match of an ELM class regex lying
entirely in a DisProt binding region (labels `2`); `NaN` elsewhere. Only ELM
classes of types $(join(options.elm_types, ", ")) are used. There is no disorder
information in this source. ELM-derived data: kept in scratch, not versioned.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "23_slim_motifs.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "ELM instances" => _relpath_from_repo(options.instances_path),
            "ELM classes index" => _relpath_from_repo(options.elms_index_path)],
        outputs = [DEFAULT_OUTPUT_FILENAME => "$(nrow(table)) proteins",
            "Proteins with an annotated ELM instance" => with_instances,
            "Proteins with a regex match in a DisProt binding region" => with_matches,
            "Proteins with at least one SLiM residue" => with_slim,
            "Annotated instances" => sum(table.n_instances),
            "Regex matches" => sum(table.n_regex_matches)],
        parameters = ["ELM class types" => join(options.elm_types, ","),
            "Maximum regex probability" => options.max_probability,
            "Excluded InstanceLogic" => join(options.excluded_instance_logic, ", ")],
        columns = _SLIM_LABELS_COLUMN_DEFINITIONS)

    println("SLiM labels written to $(output_path) (proteins=$(nrow(table))).")
    println("  proteins with an annotated ELM instance:          $(with_instances)")
    println("  proteins with a regex match in a binding region:  $(with_matches)")
    println("  proteins with at least one SLiM residue:          $(with_slim)")
    println("  annotated instances: $(sum(table.n_instances)), regex matches: $(sum(table.n_regex_matches))")
    nothing
end

main()
