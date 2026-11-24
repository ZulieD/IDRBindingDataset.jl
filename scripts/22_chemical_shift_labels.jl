#!/usr/bin/env julia

import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using ArgParse
using CSV
using DataFrames
using IDRBindingDataset:
                         CHEMICAL_SHIFT_CRITERIA,
                         CHEMICAL_SHIFT_DEFAULT_SETTINGS,
                         _CHEMICAL_SHIFT_LABELS_COLUMN_DEFINITIONS,
                         _relpath_from_repo,
                         _write_artifact_readme,
                         _write_csv_gz,
                         build_chemical_shift_labels,
                         read_final_dataset_csv,
                         remove_files

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_DATASET_PATH = normpath(joinpath(
    REPO_ROOT, "data", "dataset", "cluster_dataset_labels.csv.gz"))
const DEFAULT_MAPPING_PATH = normpath(joinpath(
    REPO_ROOT, "data", "derived", "disprot_bmrb_mapping", "disprot_bmrb_mapping.csv.gz"))
const DEFAULT_CACHE_DIR = normpath(joinpath(REPO_ROOT, "data", "scratch", "bmrb_cache"))
const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    REPO_ROOT, "data", "derived", "chemical_shift_labels"))
const OUTPUT_FILENAMES = (
    labels = "chemical_shift_labels.csv.gz",
    lists = "chemical_shift_lists.csv.gz",
    comparisons = "chemical_shift_comparisons.csv.gz"
)

"""
    parse_command_line(args)

Parse command-line options for the chemical-shift label script.
"""
function parse_command_line(args)
    defaults = CHEMICAL_SHIFT_DEFAULT_SETTINGS
    settings = ArgParseSettings(
        description = """
Build residue labels from BMRB chemical shifts for DisProt proteins: disorder
from free-state shifts, binding from comparable free/bound pairs.

Input: `disprot_bmrb_mapping.csv.gz` (21_disprot_bmrb_mapping.jl) and the BMRB
entries cached by that script. Outputs, keyed by `protein_disprot_id` and
`protein_id`:
- `chemical_shift_labels.csv.gz`: one row per DisProt protein with
  `Label_chemical_shift`, `Label_chemical_shift_disorder`,
  `Label_chemical_shift_binding` (one value per UniProt residue);
- `chemical_shift_lists.csv.gz`: every shift list, used or not, and why;
- `chemical_shift_comparisons.csv.gz`: every bound × free comparison with each
  criterion (pH, temperature, ionic strength, buffer, solvent, 15N window,
  field, spectrometer).

Run from repository root:
  julia scripts/22_chemical_shift_labels.jl

Complete example:
  julia scripts/22_chemical_shift_labels.jl --mapping_path data/derived/disprot_bmrb_mapping/disprot_bmrb_mapping.csv.gz --cache_dir data/scratch/bmrb_cache --require pH,temperature,ionic,buffer,solvent --max_delta_pH 0.5 --max_delta_temperature 10 --csp_k_sd 1.5 --output_dir data/derived/chemical_shift_labels

The script activates and instantiates `scripts/Project.toml` internally.
""",
    )

    @add_arg_table! settings begin
        "--dataset_path"
        help = "Input final dataset CSV (gzip)."
        arg_type = String
        default = DEFAULT_DATASET_PATH
        "--mapping_path"
        help = "Input DisProt to BMRB mapping CSV (gzip) from 21_disprot_bmrb_mapping.jl."
        arg_type = String
        default = DEFAULT_MAPPING_PATH
        "--cache_dir"
        help = "BMRB cache directory (`entry_<id>.json`), filled by 21_disprot_bmrb_mapping.jl."
        arg_type = String
        default = DEFAULT_CACHE_DIR
        "--min_identity"
        help = "Minimum BMRB entity / UniProt identity; mutated residues are not used."
        arg_type = Float64
        default = defaults.min_identity
        "--binding_flank"
        help = "No mutation accepted on DisProt binding residues extended by this many residues on each side."
        arg_type = Int
        default = defaults.binding_flank
        "--max_indels"
        help = "Maximum insertions/deletions in the entity alignment."
        arg_type = Int
        default = defaults.max_indels
        "--min_aligned"
        help = "Minimum aligned residues."
        arg_type = Int
        default = defaults.min_aligned
        "--lenient_numbering"
        help = "Accept shift lists whose numbering does not match the entity on 100 % of atoms."
        action = :store_true
        "--trim_ends"
        help = "Residues not analysed at each end of the NMR construct."
        arg_type = Int
        default = defaults.trim_ends
        "--random_coil"
        help = "Random coil table: ncidp (Tamiola 2010) or wishart (Wishart 1995)."
        arg_type = String
        default = string(defaults.random_coil)
        "--secondary_shift_threshold"
        help = "|smoothed Δδ(Cα) − Δδ(Cβ)| below this value (ppm) -> disordered."
        arg_type = Float64
        default = defaults.secondary_shift_threshold
        "--hn_min"
        help = "Lower bound of the amide 1H window of disordered residues (ppm)."
        arg_type = Float64
        default = defaults.hn_min
        "--hn_max"
        help = "Upper bound of the amide 1H window of disordered residues (ppm)."
        arg_type = Float64
        default = defaults.hn_max
        "--min_run"
        help = "Minimum disordered run with Cα/Cβ shifts."
        arg_type = Int
        default = defaults.min_run
        "--min_run_hn"
        help = "Minimum disordered run with amide 1H only."
        arg_type = Int
        default = defaults.min_run_hn
        "--require"
        help = "Comparability criteria that must not differ, comma-separated among $(join(CHEMICAL_SHIFT_CRITERIA, ", "))."
        arg_type = String
        default = join(defaults.required_criteria, ",")
        "--strict_metadata"
        help = "Also reject free/bound pairs with unknown (not deposited) criteria."
        action = :store_true
        "--max_delta_pH"
        help = "Maximum |ΔpH| between free and bound states."
        arg_type = Float64
        default = defaults.max_delta_pH
        "--max_delta_temperature"
        help = "Maximum |ΔT| (K) between free and bound states."
        arg_type = Float64
        default = defaults.max_delta_temperature
        "--max_delta_ionic"
        help = "Maximum ionic strength difference (mM)."
        arg_type = Float64
        default = defaults.max_delta_ionic
        "--n15_window_min"
        help = "15N windows narrower than this (ppm) must be identical in both states."
        arg_type = Float64
        default = defaults.n15_window_min
        "--csp_k_sd"
        help = "CSP threshold = mean + k·SD."
        arg_type = Float64
        default = defaults.csp_k_sd
        "--csp_min"
        help = "Minimum CSP threshold (ppm)."
        arg_type = Float64
        default = defaults.csp_min
        "--exclude_glycine_csp"
        help = "Exclude glycines from CSP."
        action = :store_true
        "--no_lost_peaks"
        help = "Do not count peaks lost in the bound state as binding."
        action = :store_true
        "--min_run_binding"
        help = "Minimum run of consecutive binding residues."
        arg_type = Int
        default = defaults.min_run_binding
        "--min_common_residues"
        help = "Residues with H/N in both free and bound states needed to compare them."
        arg_type = Int
        default = defaults.min_common_residues
        "--binding_run_mode"
        help = "separate (perturbed and lost runs counted separately) or combined."
        arg_type = String
        default = string(defaults.binding_run_mode)
        "--p_disorder_threshold"
        help = "Residue disordered when p_disorder ≥ this value."
        arg_type = Float64
        default = defaults.p_disorder_threshold
        "--output_dir"
        help = "Output directory for the CSV files and their README."
        arg_type = String
        default = DEFAULT_OUTPUT_DIR
    end

    parsed = parse_args(args, settings)
    required = Tuple(Symbol.(strip.(split(parsed["require"], ","; keepempty = false))))
    run_settings = merge(defaults,
        (min_identity = parsed["min_identity"], binding_flank = parsed["binding_flank"],
            max_indels = parsed["max_indels"],
            min_aligned = parsed["min_aligned"],
            strict_numbering = !parsed["lenient_numbering"],
            trim_ends = parsed["trim_ends"], random_coil = Symbol(parsed["random_coil"]),
            secondary_shift_threshold = parsed["secondary_shift_threshold"],
            hn_min = parsed["hn_min"], hn_max = parsed["hn_max"],
            min_run = parsed["min_run"], min_run_hn = parsed["min_run_hn"],
            required_criteria = required, strict_metadata = parsed["strict_metadata"],
            max_delta_pH = parsed["max_delta_pH"],
            max_delta_temperature = parsed["max_delta_temperature"],
            max_delta_ionic = parsed["max_delta_ionic"],
            n15_window_min = parsed["n15_window_min"], csp_k_sd = parsed["csp_k_sd"],
            csp_min = parsed["csp_min"], exclude_glycine_csp = parsed["exclude_glycine_csp"],
            use_lost_peaks = !parsed["no_lost_peaks"],
            min_run_binding = parsed["min_run_binding"],
            min_common_residues = parsed["min_common_residues"],
            binding_run_mode = Symbol(parsed["binding_run_mode"]),
            p_disorder_threshold = parsed["p_disorder_threshold"]))
    (
        dataset_path = normpath(parsed["dataset_path"]),
        mapping_path = normpath(parsed["mapping_path"]),
        cache_dir = normpath(parsed["cache_dir"]),
        settings = run_settings,
        output_dir = normpath(parsed["output_dir"])
    )
end

"Count rows per value of the given columns."
function value_counts(table::AbstractDataFrame, columns::Vector{Symbol}, prefix)
    nrow(table) == 0 && return Pair{String, Int}[]
    counts = sort!(combine(groupby(table, columns), nrow => :n), columns)
    ["$(prefix): " * join(("$(c) = $(row[c])" for c in columns), ", ") => row.n
     for row in eachrow(counts)]
end

function main()
    options = parse_command_line(ARGS)
    output_paths = map(name -> joinpath(options.output_dir, name), OUTPUT_FILENAMES)
    readme_path = joinpath(options.output_dir, "README.md")
    remove_files(vcat(collect(output_paths), [readme_path]))

    dataset = read_final_dataset_csv(options.dataset_path)
    mapping = DataFrame(CSV.File(options.mapping_path;
        types = Dict(:protein_disprot_id => String, :protein_id => String,
            :bmrb_id => String, :shift_list_id => String, :status => String,
            :state => String)))

    result = build_chemical_shift_labels(dataset, mapping; cache_dir = options.cache_dir,
        settings = options.settings)
    for name in keys(OUTPUT_FILENAMES)
        _write_csv_gz(result[name], output_paths[name])
    end

    comparisons = result.comparisons
    rejected = comparisons[.!comparisons.accepted, :]
    criterion_counts = ["rejected comparisons with different $(c)" =>
                            count(==("different"), rejected[!, c])
                        for c in CHEMICAL_SHIFT_CRITERIA]
    _write_artifact_readme(readme_path;
        title = "Residue Labels From BMRB Chemical Shifts",
        description = """
Disorder is called on free-state shift lists from secondary shifts
Δδ(Cα) − Δδ(Cβ) (random coil: $(options.settings.random_coil)) and the amide
1H window; binding on pairs of a bound list with its most similar comparable
free list, from chemical-shift perturbations (CSP) and lost peaks in runs of
consecutive residues. `$(OUTPUT_FILENAMES.lists)` gives, for every shift list,
whether it was used and why; `$(OUTPUT_FILENAMES.comparisons)` gives every
bound × free comparison with each criterion (`same`, `different`, `unknown`),
`accepted` (no required criterion differs) and `selected` (pair used), and
the buffering species of both samples (buffers are compared on their species
only: additives and salts are ignored, concentrations may differ). Shift lists
need an entity identical to UniProt on the DisProt binding region ±
$(options.settings.binding_flank) residues and at least
$(options.settings.min_identity) identical overall; mutated residues elsewhere
are not used.
""",
        script_relpath = _relpath_from_repo(joinpath(@__DIR__, "22_chemical_shift_labels.jl")),
        inputs = ["Dataset" => _relpath_from_repo(options.dataset_path),
            "DisProt to BMRB mapping" => _relpath_from_repo(options.mapping_path),
            "BMRB cache" => _relpath_from_repo(options.cache_dir)],
        outputs = vcat(
            [name => "$(nrow(result[key])) rows" for (key, name) in pairs(OUTPUT_FILENAMES)],
            value_counts(result.lists, [:state, :used], "shift lists"),
            value_counts(result.lists, [:state, :reason], "shift lists"),
            ["accepted comparisons" => count(comparisons.accepted)],
            criterion_counts,
            ["Proteins that failed" => join(result.failed_proteins, ", ")]),
        parameters = [string(k) => v for (k, v) in pairs(options.settings)],
        columns = _CHEMICAL_SHIFT_LABELS_COLUMN_DEFINITIONS)

    println("Chemical shift labels written to $(options.output_dir) (proteins=$(nrow(result.labels)), shift lists used=$(count(result.lists.used))/$(nrow(result.lists)), failed=$(length(result.failed_proteins))).")
    nothing
end

main()
