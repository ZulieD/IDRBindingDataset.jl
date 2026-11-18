#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         DEFAULT_CHECKPOINT_BATCH_SIZE,
                         _PDB_MISSING_RESIDUES_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_dataset_pdb_segments,
                         build_pdb_missing_residues_table,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_SIFTS_CHAIN_UNIPROT_PATH = normpath(joinpath(
    REPO_ROOT, "data", "raw", "sifts", "pdb_chain_uniprot.csv.gz"))
const DEFAULT_SIFTS_DOWNLOAD_DIR = normpath(joinpath(
    REPO_ROOT, "data", "scratch", "sifts_xml"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "derived", "disprot_pdb_mapping"))
const DEFAULT_OUTPUT_FILENAME = "disprot_pdb_mapping.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the DisProt to PDB mapping script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Link DisProt proteins of the final dataset to their PDB structures and count
unobserved (missing) residues. One row per (DisProt protein, PDB entry, UniProt
range), keyed by `protein_disprot_id` and `protein_id`.

1) join `data/dataset/cluster_dataset_labels.csv.gz` with the SIFTS
   chain-level mapping (`pdb_chain_uniprot.csv.gz`) on the UniProt accession,
2) merge chains covering the same UniProt range (`CHAIN`),
3) for each PDB segment, read the SIFTS XML file of the entry and count the
   residues without coordinates, split by DisProt label, keeping the chain(s)
   with the fewest missing residues (`CHAIN_KEPT`).

PDB entries are processed in parallel (start Julia with `--threads`). With
`--checkpoint_dir`, results are saved after every batch and a new run resumes
where the previous one stopped.

Run from repository root:
  julia --threads=8 scripts/17_disprot_pdb_mapping.jl --sifts_xml_dir /path/to/sifts_xml

Complete example:
  julia --threads=8 scripts/17_disprot_pdb_mapping.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --sifts_chain_uniprot_path data/raw/sifts/pdb_chain_uniprot.csv.gz --sifts_xml_dir /alpha/database/sift_new-2025-10-27/sift_files --download_missing_sifts --sifts_download_dir data/scratch/sifts_xml --checkpoint_dir data/scratch/checkpoints/17_disprot_pdb_mapping --batch_size 15 --output_dir data/derived/disprot_pdb_mapping

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip) with `protein_id`, `sequence` and `residue_labels`."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--sifts_chain_uniprot_path"
        help = "Input SIFTS chain-level PDB to UniProt mapping (gzip CSV)."
        arg_type = String
        default = DEFAULT_SIFTS_CHAIN_UNIPROT_PATH
        "--sifts_xml_dir"
        help = "Directory with SIFTS XML files (`<pdb>.xml.gz`, flat or in `ab/` subfolders). Can be repeated."
        arg_type = String
        action = :append_arg
        "--download_missing_sifts"
        help = "Download SIFTS XML files that are not found into `--sifts_download_dir`."
        action = :store_true
        "--sifts_download_dir"
        help = "Directory for downloaded SIFTS XML files (also searched)."
        arg_type = String
        default = DEFAULT_SIFTS_DOWNLOAD_DIR
        "--checkpoint_dir"
        help = "Directory used to save and resume progress (empty: no checkpoint)."
        arg_type = String
        default = ""
        "--batch_size"
        help = "Number of PDB entries per checkpointed batch."
        arg_type = Int
        default = DEFAULT_CHECKPOINT_BATCH_SIZE
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
        sifts_chain_uniprot_path = normpath(parsed["sifts_chain_uniprot_path"]),
        sifts_xml_dirs = normpath.(String.(parsed["sifts_xml_dir"])),
        sifts_download_dir = parsed["download_missing_sifts"] ?
                             normpath(parsed["sifts_download_dir"]) : nothing,
        checkpoint_dir = isempty(parsed["checkpoint_dir"]) ? nothing :
                         normpath(parsed["checkpoint_dir"]),
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
    sifts_chain_uniprot = DataFrame(CSV.File(options.sifts_chain_uniprot_path;
        header = 2, types = String))
    segments = build_dataset_pdb_segments(dataset, sifts_chain_uniprot)
    println("PDB segments: $(nrow(segments)) for $(length(unique(segments.protein_id))) protein(s).")

    result = build_pdb_missing_residues_table(dataset, segments;
        sifts_xml_dirs = options.sifts_xml_dirs,
        sifts_download_dir = options.sifts_download_dir,
        checkpoint_dir = options.checkpoint_dir,
        batch_size = options.batch_size)

    _write_csv_gz(result.table, output_path)
    _write_artifact_readme(readme_path;
        title = "DisProt To PDB Mapping With Missing Residues",
        description = """
One row per (DisProt protein, PDB entry, UniProt range) for dataset proteins
with PDB structures. `protein_disprot_id` and `protein_id` link this table to
the final dataset and to the other derived tables. Missing residues are residues of the deposited sequence without
coordinates, according to the SIFTS XML residue mapping; only residues mapped
to the protein UniProt accession are counted.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "17_disprot_pdb_mapping.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "SIFTS chain mapping" => _relpath_from_repo(options.sifts_chain_uniprot_path),
            "SIFTS XML directories" => join(
                vcat(options.sifts_xml_dirs,
                    something(options.sifts_download_dir, String[])), ", ")],
        outputs = [options.output_filename => "$(nrow(result.table)) rows",
            "PDB entries that failed" => join(result.failed_pdbs, ", ")],
        columns = _PDB_MISSING_RESIDUES_COLUMN_DEFINITIONS)

    println("DisProt to PDB mapping written to $(output_path) (rows=$(nrow(result.table)), failed PDB entries=$(length(result.failed_pdbs))).")
    nothing
end

main()
