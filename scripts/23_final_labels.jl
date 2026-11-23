#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         FINAL_LABEL_SOURCES,
                         _FINAL_LABELS_COLUMN_DEFINITIONS,
                         _FINAL_LABEL_SOURCES_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_final_labels_table,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_UNIPROT_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "uniprot_binding_sites", "uniprot_binding_sites.csv.gz"))
const DEFAULT_AGGREGATED_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "aggregated_labels", "aggregated_labels.csv.gz"))
const DEFAULT_CHEMICAL_SHIFT_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "chemical_shift_labels", "chemical_shift_labels.csv.gz"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "derived", "final_labels"))
const DEFAULT_OUTPUT_FILENAME = "final_labels.csv.gz"
const DEFAULT_SOURCES_FILENAME = "final_labels_sources.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the final-label script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Combine the four binding evidence sources into final residue labels, one row
per dataset protein keyed by `protein_disprot_id` and `protein_id`:
DisProt labels, UniProt binding sites (18_uniprot_binding_sites.jl), PDB
contacts (`Label_PDB_binding`, 20_aggregated_labels.jl) and chemical shifts
(`Label_chemical_shift_disorder` and `Label_chemical_shift_binding`,
22_chemical_shift_labels.jl).

Per residue:
  Label_aggregate_all_disorder d = mean disorder vote of DisProt (0/1) and NMR
                               (0.5 when DisProt says disordered and NMR ordered)
  Label_aggregate_all_binding  b = x / n over the n binding sources covering the
                               residue (1/4 each when all four cover it), on d > 0
  Label_aggregate_all          d + b in [0, 2]: 0 ordered, 1 disordered,
                               2 disordered binding

Run from repository root:
  julia scripts/23_final_labels.jl

Complete example:
  julia scripts/23_final_labels.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --uniprot_binding_sites_path data/derived/uniprot_binding_sites/uniprot_binding_sites.csv.gz --aggregated_labels_path data/derived/aggregated_labels/aggregated_labels.csv.gz --chemical_shift_labels_path data/derived/chemical_shift_labels/chemical_shift_labels.csv.gz --output_dir data/derived/final_labels

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip)."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--uniprot_binding_sites_path"
        help = "Input UniProt binding sites CSV (gzip) from 18_uniprot_binding_sites.jl."
        arg_type = String
        default = DEFAULT_UNIPROT_PATH
        "--aggregated_labels_path"
        help = "Input aggregated labels CSV (gzip) from 20_aggregated_labels.jl."
        arg_type = String
        default = DEFAULT_AGGREGATED_PATH
        "--chemical_shift_labels_path"
        help = "Input chemical shift labels CSV (gzip) from 22_chemical_shift_labels.jl."
        arg_type = String
        default = DEFAULT_CHEMICAL_SHIFT_PATH
        "--output_dir"
        help = "Output directory for the CSV files and their README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
    end

    parsed = parse_args(args, settings)
    (
        dataset_path = normpath(parsed["dataset_path"]),
        uniprot_binding_sites_path = normpath(parsed["uniprot_binding_sites_path"]),
        aggregated_labels_path = normpath(parsed["aggregated_labels_path"]),
        chemical_shift_labels_path = normpath(parsed["chemical_shift_labels_path"]),
        output_dir = normpath(parsed["output_dir"])
    )
end

"""
    read_string_csv(path, columns)

Read a derived CSV with `columns` forced to `String`.
"""
function read_string_csv(path::AbstractString, columns)
    DataFrame(CSV.File(path; types = Dict(column => String for column in columns)))
end

function main()
    options = parse_command_line(ARGS)
    output_path = joinpath(options.output_dir, DEFAULT_OUTPUT_FILENAME)
    sources_path = joinpath(options.output_dir, DEFAULT_SOURCES_FILENAME)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files([output_path, sources_path, readme_path])

    dataset = read_final_dataset_csv(options.dataset_path)
    uniprot = read_string_csv(options.uniprot_binding_sites_path,
        [:protein_disprot_id, :protein_id, :residues_binding_in_uniprot])
    aggregated = read_string_csv(options.aggregated_labels_path,
        [:protein_disprot_id, :protein_id, :Label_PDB_binding])
    chemical_shift = read_string_csv(options.chemical_shift_labels_path,
        [:protein_disprot_id, :protein_id, :Label_chemical_shift_binding])

    result = build_final_labels_table(dataset, uniprot, aggregated, chemical_shift)
    table, sources = result.table, result.sources
    _write_csv_gz(table, output_path)
    _write_csv_gz(sources, sources_path)

    source_counts = [name => count(s -> name in split(s, ";"), sources.sources)
                     for name in FINAL_LABEL_SOURCES]
    combination_counts = sort!(combine(groupby(sources, :sources), nrow => :n), :n;
        rev = true)
    coverage = [
        "Disordered residues (d > 0)" => sum(sources.n_disordered),
        "DisProt / NMR disorder disagreement (0 < d < 1)" =>
            sum(sources.n_partial_disorder),
        "Disordered for NMR only" => sum(sources.n_nmr_only_disorder),
        "Disordered residues, DisProt only binding source" => sum(sources.n_1_source),
        "Disordered residues, 2 sources" => sum(sources.n_2_sources),
        "Disordered residues, 3 sources" => sum(sources.n_3_sources),
        "Disordered residues, 4 sources" => sum(sources.n_4_sources)]

    _write_artifact_readme(readme_path;
        title = "Final Residue Binding Labels (DisProt, UniProt, PDB, Chemical Shift)",
        description = """
One row per dataset protein. Each vector column has one value per residue of
the UniProt sequence. Disorder `d` is the mean disorder vote of DisProt (`0`/`1`)
and NMR (`Label_chemical_shift_disorder`, when evaluable). On residues with
`d > 0`, binding `b = x / n` where `n` is the number of sources covering the
residue and `x` the sum of their binding votes: DisProt (`1` if labelled binding), UniProt (`1` inside a binding
site, when the protein has UniProt binding sites), PDB (`Label_PDB_binding`,
fraction of used PDB entries in contact with a heteromeric partner) and chemical
shift (`Label_chemical_shift_binding`, fraction of free/bound pairs calling the
residue binding). Each covering source weighs the same (`1/4` when all four
cover the residue). `Label_aggregate_all = d + b` in `[0, 2]`; residues with
`d = 0` get `0`. `$(DEFAULT_SOURCES_FILENAME)` lists, per protein, the sources
available, DisProt/NMR disorder disagreements and how many disordered residues
are covered by 1 to 4 binding sources.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "23_final_labels.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "UniProt binding sites" =>
                _relpath_from_repo(options.uniprot_binding_sites_path),
            "Aggregated labels (PDB)" => _relpath_from_repo(options.aggregated_labels_path),
            "Chemical shift labels" =>
                _relpath_from_repo(options.chemical_shift_labels_path)],
        outputs = vcat(
            [DEFAULT_OUTPUT_FILENAME => "$(nrow(table)) proteins",
                DEFAULT_SOURCES_FILENAME => "$(nrow(sources)) proteins"],
            ["Proteins with source $(name)" => n for (name, n) in source_counts],
            coverage),
        parameters = [
            "Source weighting" => "equal, 1/n over the sources covering the residue",
            "Final label" => "disorder + binding"],
        columns = vcat(_FINAL_LABELS_COLUMN_DEFINITIONS,
            ["($(DEFAULT_SOURCES_FILENAME)) $(column)" => definition
             for (column, definition) in _FINAL_LABEL_SOURCES_COLUMN_DEFINITIONS]))

    println("Final labels written to $(output_path) (proteins=$(nrow(table))).")
    println("Proteins per source:")
    for (name, n) in source_counts
        println("  ", rpad(name, 16), n)
    end
    println("Proteins per source combination:")
    for row in eachrow(combination_counts)
        println("  ", rpad(row.sources, 40), row.n)
    end
    println("Disordered residues per number of covering sources:")
    for (label, n) in coverage
        println("  ", rpad(label, 50), n)
    end
    nothing
end

main()
