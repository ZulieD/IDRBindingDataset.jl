#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using DataFrames
using Dates
using IDRBindingDataset:
                         BMRB_API_URL,
                         BMRB_UNIPROT_MAPPING_FILENAME,
                         _BMRB_MAPPING_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_disprot_bmrb_mapping_table,
                         fetch_bmrb_uniprot_mapping,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_CACHE_DIR = normpath(joinpath(REPO_ROOT, "data", "scratch", "bmrb_cache"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "derived", "disprot_bmrb_mapping"))
const DEFAULT_OUTPUT_FILENAME = "disprot_bmrb_mapping.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the DisProt to BMRB mapping script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Link DisProt proteins of the final dataset to their BMRB NMR entries and
describe each assigned chemical-shift list: matched entity, free / bound /
homo-oligomer state, partner molecules, sample conditions (pH, temperature,
ionic strength), buffer, solvent, spectrometer and shift counts. One row per
(DisProt protein, BMRB entry, shift list), keyed by `protein_disprot_id` and
`protein_id`. Nothing is filtered: the comparison of free and bound states is
done in a later step.

The UniProt -> BMRB mapping and every entry are downloaded from the BMRB API
into `--cache_dir`; a new run only downloads what is missing.

Run from repository root:
  julia scripts/21_disprot_bmrb_mapping.jl

Complete example:
  julia scripts/21_disprot_bmrb_mapping.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --cache_dir data/scratch/bmrb_cache --output_dir data/derived/disprot_bmrb_mapping --output_filename disprot_bmrb_mapping.csv.gz

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip)."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--cache_dir"
        help = "Directory for downloaded BMRB files (mapping and `entry_<id>.json`)."
        arg_type = String
        default = DEFAULT_CACHE_DIR
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
        cache_dir = normpath(parsed["cache_dir"]),
        output_dir = normpath(parsed["output_dir"]),
        output_filename = String(parsed["output_filename"])
    )
end

"Count rows per value of `column`, sorted by decreasing count."
function value_counts(table::AbstractDataFrame, column::Symbol)
    counts = combine(groupby(table, column), nrow => :n)
    sort!(counts, :n; rev = true)
    ["$(column) = $(coalesce(row[column], "missing"))" => row.n for row in eachrow(counts)]
end

function main()
    options = parse_command_line(ARGS)
    output_path = joinpath(options.output_dir, options.output_filename)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files([output_path, readme_path])

    dataset = read_final_dataset_csv(options.dataset_path)
    mapping = fetch_bmrb_uniprot_mapping(; cache_dir = options.cache_dir)
    mapping_path = joinpath(options.cache_dir, BMRB_UNIPROT_MAPPING_FILENAME)
    mapping_date = Dates.format(unix2datetime(mtime(mapping_path)), "yyyy-mm-dd")
    n_mapped = count(id -> haskey(mapping, uppercase(id)), dataset.protein_id)
    println("UniProt accessions with BMRB entries: $(length(mapping)); dataset proteins mapped: $(n_mapped)/$(nrow(dataset)).")

    result = build_disprot_bmrb_mapping_table(dataset, mapping;
        cache_dir = options.cache_dir)
    table = result.table

    _write_csv_gz(table, output_path)
    _write_artifact_readme(readme_path;
        title = "DisProt To BMRB Mapping With NMR Experiment Metadata",
        description = """
One row per (DisProt protein, BMRB entry, assigned chemical-shift list). BMRB
entries come from the BMRB UniProt -> BMRB mapping and are read from the BMRB
API v2 (NMR-STAR JSON). The BMRB entity aligned best to the UniProt sequence is
taken as the protein; the other molecules of the sample assembly are its
partners (`state = bound`), copies of the protein make a homo-oligomer. Sample
conditions, buffer, solvent and spectrometer are those linked to each shift
list. Nothing is filtered: compare free and bound states in a later step.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "21_disprot_bmrb_mapping.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "BMRB API" => BMRB_API_URL,
            "UniProt -> BMRB mapping downloaded on" => mapping_date,
            "BMRB cache" => _relpath_from_repo(options.cache_dir)],
        outputs = vcat(
            [options.output_filename => "$(nrow(table)) rows",
                "Dataset proteins with BMRB entries" => n_mapped,
                "Proteins in the table" => length(unique(table.protein_id))],
            value_counts(table, :status), value_counts(table, :state),
            ["Entries that failed" => join(result.failed_entries, ", ")]),
        columns = _BMRB_MAPPING_COLUMN_DEFINITIONS)

    println("DisProt to BMRB mapping written to $(output_path) (rows=$(nrow(table)), failed entries=$(length(result.failed_entries))).")
    nothing
end

main()
