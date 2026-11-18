#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         DEFAULT_CHECKPOINT_BATCH_SIZE,
                         PDB_CONTACT_CUTOFF_ANGSTROM,
                         _PDB_CONTACTS_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_pdb_contacts_table,
                         read_final_dataset_csv,
                         read_pdb_quality_table,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_MAPPING_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "disprot_pdb_mapping", "disprot_pdb_mapping.csv.gz"))
const DEFAULT_SCRATCH_DIR = normpath(joinpath(REPO_ROOT, "data", "scratch"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(REPO_ROOT, "data", "derived", "pdb_contacts"))
const DEFAULT_OUTPUT_FILENAME = "pdb_contacts.csv.gz"

"""
    parse_command_line(args)

Parse command-line options for the PDB contacts script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Compute heteromeric interface contacts of DisProt proteins in their PDB
structures. One row per (DisProt protein, PDB entry), keyed by
`protein_disprot_id` and `protein_id`.

For each pair, the protein chains come from `disprot_pdb_mapping.csv.gz`
(17_disprot_pdb_mapping.jl). Partner chains are the other polymer chains that
are not copies of the protein. A protein residue is in contact when one of its
heavy atoms is within `--contact_cutoff` Å of a heavy atom of a partner chain.
Residues are mapped to UniProt positions with SIFTS XML. All contacts are kept
(no disorder filter) and no PDB quality filter is applied: both are applied
when labels are aggregated.

By default only proteins with at least one DisProt binding residue (label `2`)
are processed. Pairs are processed in parallel (start Julia with `--threads`)
and `--checkpoint_dir` allows resuming an interrupted run.

Run from repository root:
  julia --threads=16 scripts/19_pdb_contacts.jl --mmcif_dir /path/to/mmcif --sifts_xml_dir /path/to/sifts_xml

Complete example:
  julia --threads=16 scripts/19_pdb_contacts.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --mapping_path data/derived/disprot_pdb_mapping/disprot_pdb_mapping.csv.gz --mmcif_dir /alpha/database/pdb/mmcif_files --sifts_xml_dir /alpha/database/sift_new-2025-10-27/sift_files --download_missing --pdb_quality_path pdb_information_details.csv --contact_cutoff 6.0 --checkpoint_dir data/scratch/checkpoints/19_pdb_contacts --batch_size 50 --output_dir data/derived/pdb_contacts

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip)."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--mapping_path"
        help = "Input DisProt to PDB mapping CSV (gzip) from 17_disprot_pdb_mapping.jl."
        arg_type = String
        default = DEFAULT_MAPPING_PATH
        "--mmcif_dir"
        help = "Directory with mmCIF files (`<pdb>.cif[.gz]`, flat or in `ab/` subfolders). Can be repeated."
        arg_type = String
        action = :append_arg
        "--sifts_xml_dir"
        help = "Directory with SIFTS XML files (`<pdb>.xml.gz`). Can be repeated."
        arg_type = String
        action = :append_arg
        "--download_missing"
        help = "Download missing mmCIF (RCSB) and SIFTS XML files into `--scratch_dir`."
        action = :store_true
        "--scratch_dir"
        help = "Directory for downloaded files (`mmcif/`, `sifts_xml/`)."
        arg_type = String
        default = DEFAULT_SCRATCH_DIR
        "--pdb_quality_path"
        help = "Optional CSV with PDB `METHOD` and `RESOLUTION` columns, copied into the output (empty: none)."
        arg_type = String
        default = ""
        "--pdb_quality_id_column"
        help = "PDB identifier column of `--pdb_quality_path` (empty: auto-detect)."
        arg_type = String
        default = ""
        "--all_proteins"
        help = "Process every protein with a PDB, not only proteins with DisProt binding residues."
        action = :store_true
        "--contact_cutoff"
        help = "Heavy-atom distance cutoff for interface contacts (Å)."
        arg_type = Float64
        default = PDB_CONTACT_CUTOFF_ANGSTROM
        "--checkpoint_dir"
        help = "Directory used to save and resume progress (empty: no checkpoint)."
        arg_type = String
        default = ""
        "--batch_size"
        help = "Number of protein/PDB pairs per checkpointed batch."
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
    optional_text(key) = isempty(parsed[key]) ? nothing : String(parsed[key])
    scratch_dir = normpath(parsed["scratch_dir"])
    download = parsed["download_missing"]
    (
        dataset_path = normpath(parsed["dataset_path"]),
        mapping_path = normpath(parsed["mapping_path"]),
        mmcif_dirs = normpath.(String.(parsed["mmcif_dir"])),
        sifts_xml_dirs = normpath.(String.(parsed["sifts_xml_dir"])),
        mmcif_download_dir = download ? joinpath(scratch_dir, "mmcif") : nothing,
        sifts_download_dir = download ? joinpath(scratch_dir, "sifts_xml") : nothing,
        pdb_quality_path = optional_text("pdb_quality_path") === nothing ? nothing :
                           normpath(parsed["pdb_quality_path"]),
        pdb_quality_id_column = optional_text("pdb_quality_id_column"),
        only_binding_proteins = !parsed["all_proteins"],
        contact_cutoff = Float64(parsed["contact_cutoff"]),
        checkpoint_dir = optional_text("checkpoint_dir") === nothing ? nothing :
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
    mapping = DataFrame(CSV.File(options.mapping_path;
        types = Dict(:protein_disprot_id => String, :protein_id => String,
            :PDB => String, :CHAIN => String)))
    quality = options.pdb_quality_path === nothing ? Dict{String, Any}() :
              read_pdb_quality_table(options.pdb_quality_path;
        id_column = options.pdb_quality_id_column)

    result = build_pdb_contacts_table(dataset, mapping;
        mmcif_dirs = options.mmcif_dirs, mmcif_download_dir = options.mmcif_download_dir,
        sifts_xml_dirs = options.sifts_xml_dirs,
        sifts_download_dir = options.sifts_download_dir,
        pdb_quality_table = quality,
        only_binding_proteins = options.only_binding_proteins,
        checkpoint_dir = options.checkpoint_dir, batch_size = options.batch_size,
        contact_cutoff = options.contact_cutoff)
    table = result.table

    _write_csv_gz(table, output_path)
    status_counts = combine(groupby(table, :status), nrow => :n)
    sort!(status_counts, :status)
    _write_artifact_readme(readme_path;
        title = "Heteromeric PDB Interface Contacts Of DisProt Proteins",
        description = """
One row per (DisProt protein, PDB entry). Protein chains come from
`disprot_pdb_mapping.csv.gz`; partner chains are the other polymer chains of
the first model that are not copies of the protein (SIFTS chains, or k-mer
identity to the protein sequence). A protein residue is in contact when one of
its heavy atoms (hydrogens removed) is within the contact cutoff of a heavy
atom of a partner chain. Positions are UniProt positions (SIFTS XML mapping),
written as ranges such as `12-40;55`. All contacts are reported, without
disorder or PDB quality filter.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "19_pdb_contacts.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "DisProt to PDB mapping" => _relpath_from_repo(options.mapping_path),
            "mmCIF directories" => join(
                vcat(options.mmcif_dirs,
                    something(options.mmcif_download_dir, String[])), ", "),
            "SIFTS XML directories" => join(
                vcat(options.sifts_xml_dirs,
                    something(options.sifts_download_dir, String[])), ", "),
            "PDB quality table" => something(options.pdb_quality_path, "none")],
        outputs = vcat([options.output_filename => "$(nrow(table)) rows"],
            ["status = $(row.status)" => row.n for row in eachrow(status_counts)],
            ["Pairs that failed" => join(result.failed_items, ", ")]),
        parameters = ["Contact cutoff (Å)" => options.contact_cutoff,
            "Only proteins with DisProt binding residues" => options.only_binding_proteins],
        columns = _PDB_CONTACTS_COLUMN_DEFINITIONS)

    println("PDB contacts written to $(output_path) (rows=$(nrow(table)), failed=$(length(result.failed_items))).")
    nothing
end

main()
