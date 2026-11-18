#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using Plots
using IDRBindingDataset:
                         PDB_LABEL_MAX_RESOLUTION,
                         PDB_LABEL_VOTES,
                         _AGGREGATED_LABELS_COLUMN_DEFINITIONS,
                         _parse_label_vector,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_aggregated_labels_table,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_UNIPROT_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "uniprot_binding_sites", "uniprot_binding_sites.csv.gz"))
const DEFAULT_CONTACTS_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "pdb_contacts", "pdb_contacts.csv.gz"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "derived", "aggregated_labels"))
const DEFAULT_OUTPUT_FILENAME = "aggregated_labels.csv.gz"
const DEFAULT_PDB_LOG_FILENAME = "aggregated_labels_pdb_log.csv.gz"
# Continuous colour scale on [0, 2]: 0 grey, 1 gold, 2 red (white = not covered).
const LABEL_COLORS = cgrad([:grey80, :gold, :red], [0.0, 0.5, 1.0])

"""
    parse_command_line(args)

Parse command-line options for the aggregated-label script.
"""
function parse_command_line(args)
    settings = ArgParseSettings(
        description = """
Aggregate DisProt labels, UniProt binding sites (18_uniprot_binding_sites.jl)
and heteromeric PDB contacts (19_pdb_contacts.jl) into residue-level binding
labels. One row per DisProt protein with at least one used PDB entry, keyed by
`protein_disprot_id` and `protein_id`.

PDB entries are filtered here: status (`ok`, plus homo-oligomers with
`--keep_pdb_without_partner`) and quality (`--max_resolution`,
`--drop_no_resolution`, `--drop_no_info`). PDB contacts are only counted on
DisProt disordered residues (labels `1`/`2`).

Run from repository root:
  julia scripts/20_aggregated_labels.jl

Complete example:
  julia scripts/20_aggregated_labels.jl --dataset_path data/dataset/cluster_dataset_labels.csv.gz --uniprot_binding_sites_path data/derived/uniprot_binding_sites/uniprot_binding_sites.csv.gz --pdb_contacts_path data/derived/pdb_contacts/pdb_contacts.csv.gz --max_resolution 5.0 --pdb_vote fraction --plot_dir data/scratch/plots/aggregated_labels --output_dir data/derived/aggregated_labels

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
        "--pdb_contacts_path"
        help = "Input PDB contacts CSV (gzip) from 19_pdb_contacts.jl."
        arg_type = String
        default = DEFAULT_CONTACTS_PATH
        "--max_resolution"
        help = "Use PDB entries with resolution strictly below this value (Å)."
        arg_type = Float64
        default = PDB_LABEL_MAX_RESOLUTION
        "--drop_no_resolution"
        help = "Do not use PDB entries without resolution (e.g. NMR)."
        action = :store_true
        "--drop_no_info"
        help = "Do not use PDB entries absent from the PDB quality table."
        action = :store_true
        "--keep_pdb_without_partner"
        help = "Also use PDB entries without heteromeric partner (observed residues, no contact)."
        action = :store_true
        "--pdb_vote"
        help = "PDB contribution to the final rate: fraction, any or majority."
        arg_type = String
        default = "fraction"
        "--plot_dir"
        help = "Directory for one PNG per protein (empty: no plots)."
        arg_type = String
        default = ""
        "--output_dir"
        help = "Output directory for the CSV files and their README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
    end

    parsed = parse_args(args, settings)
    vote = Symbol(parsed["pdb_vote"])
    vote in PDB_LABEL_VOTES ||
        throw(ArgumentError("--pdb_vote must be one of $(join(PDB_LABEL_VOTES, ", "))."))
    (
        dataset_path = normpath(parsed["dataset_path"]),
        uniprot_binding_sites_path = normpath(parsed["uniprot_binding_sites_path"]),
        pdb_contacts_path = normpath(parsed["pdb_contacts_path"]),
        max_resolution = Float64(parsed["max_resolution"]),
        keep_no_resolution = !parsed["drop_no_resolution"],
        keep_no_info = !parsed["drop_no_info"],
        keep_pdb_without_partner = parsed["keep_pdb_without_partner"],
        vote = vote,
        plot_dir = isempty(parsed["plot_dir"]) ? nothing : normpath(parsed["plot_dir"]),
        output_dir = normpath(parsed["output_dir"])
    )
end

"""
    plot_protein_labels(row, output_path)

Draw the DisProt, UniProt, aggregated PDB and final label tracks of one
aggregated-labels row on a common `[0, 2]` scale (white = not covered).
"""
function plot_protein_labels(row, output_path::AbstractString)
    pdb_n = _parse_label_vector(row.pdb_N)
    pdb_p = _parse_label_vector(row.pdb_proba)
    pdb_track = [n == 0 ? NaN : p > 0 ? 1.0 + p : 0.0 for (n, p) in zip(pdb_n, pdb_p)]
    tracks = [
        "DisProt" => _parse_label_vector(row.disprot),
        (row.has_uniprot ? "UniProt binding" : "UniProt (no annotation)") =>
            _parse_label_vector(row.uniprot),
        "PDB aggregated ($(row.n_pdb) PDB)" => pdb_track,
        "Final label (1 + rate)" => _parse_label_vector(row.final_label)
    ]
    matrix = permutedims(reduce(hcat, last.(tracks)))
    n_tracks, n_positions = size(matrix)
    figure = heatmap(1:n_positions, 1:n_tracks, matrix;
        c = LABEL_COLORS, clims = (0, 2), colorbar = true,
        colorbar_title = "label (0 grey · 1 gold · 2 red)",
        yflip = true, yticks = (1:n_tracks, first.(tracks)), ytickfontsize = 8,
        xlabel = "UniProt position",
        title = "$(row.protein_disprot_id) / $(row.protein_id)   (white = not covered)",
        size = (clamp(3n_positions, 1000, 3200), 160 + 40n_tracks),
        left_margin = 45Plots.mm, bottom_margin = 8Plots.mm, right_margin = 8Plots.mm,
        background_color_inside = :white, framestyle = :box)
    hline!(figure, collect(0.5:1:(n_tracks + 0.5)); color = :white, lw = 2, label = false)
    savefig(figure, output_path)
    nothing
end

function main()
    options = parse_command_line(ARGS)
    output_path = joinpath(options.output_dir, DEFAULT_OUTPUT_FILENAME)
    pdb_log_path = joinpath(options.output_dir, DEFAULT_PDB_LOG_FILENAME)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files([output_path, pdb_log_path, readme_path])

    dataset = read_final_dataset_csv(options.dataset_path)
    uniprot = DataFrame(CSV.File(options.uniprot_binding_sites_path;
        types = Dict(:protein_disprot_id => String, :protein_id => String,
            :residues_binding_in_uniprot => String)))
    contacts = DataFrame(CSV.File(options.pdb_contacts_path;
        types = Dict(:protein_disprot_id => String, :protein_id => String,
            :PDB => String, :status => String, :quality_category => String,
            :resolution => Float64, :resolved_positions => String,
            :contact_positions => String)))

    result = build_aggregated_labels_table(dataset, uniprot, contacts;
        max_resolution = options.max_resolution,
        keep_no_resolution = options.keep_no_resolution,
        keep_no_info = options.keep_no_info,
        keep_pdb_without_partner = options.keep_pdb_without_partner, vote = options.vote)
    table, pdb_log = result.table, result.pdb_log

    _write_csv_gz(table, output_path)
    _write_csv_gz(pdb_log, pdb_log_path)

    if options.plot_dir !== nothing
        mkpath(options.plot_dir)
        for row in eachrow(table)
            try
                plot_protein_labels(row,
                    joinpath(options.plot_dir, "$(row.protein_disprot_id).png"))
            catch err
                println("Plot failed for $(row.protein_disprot_id): $(sprint(showerror, err))")
            end
        end
    end

    reason_counts = sort!(combine(groupby(pdb_log, :reason), nrow => :n), :reason)
    _write_artifact_readme(readme_path;
        title = "Residue Binding Labels Aggregated Per DisProt Protein",
        description = """
One row per DisProt protein with at least one used PDB entry. Each vector column
has one value per residue of the UniProt sequence (`NaN` = not covered).
PDB evidence counts, for each residue, the used PDB entries where it is
observed (`pdb_N`) and the fraction where it contacts a heteromeric partner
(`pdb_proba`); contacts are only counted on DisProt disordered residues.
`$(DEFAULT_PDB_LOG_FILENAME)` records, for every protein/PDB pair of the
contacts table, whether it was used and why.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "20_aggregated_labels.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "UniProt binding sites" =>
                _relpath_from_repo(options.uniprot_binding_sites_path),
            "PDB contacts" => _relpath_from_repo(options.pdb_contacts_path)],
        outputs = vcat(
            [DEFAULT_OUTPUT_FILENAME => "$(nrow(table)) proteins",
                DEFAULT_PDB_LOG_FILENAME => "$(nrow(pdb_log)) protein/PDB pairs"],
            ["PDB pairs: $(row.reason)" => row.n for row in eachrow(reason_counts)]),
        parameters = ["Maximum resolution (Å, strict)" => options.max_resolution,
            "Use PDB without resolution" => options.keep_no_resolution,
            "Use PDB absent from quality table" => options.keep_no_info,
            "Use PDB without heteromeric partner" => options.keep_pdb_without_partner,
            "PDB vote" => options.vote],
        columns = _AGGREGATED_LABELS_COLUMN_DEFINITIONS)

    println("Aggregated labels written to $(output_path) (proteins=$(nrow(table)), PDB pairs used=$(count(pdb_log.used))/$(nrow(pdb_log))).")
    nothing
end

main()
