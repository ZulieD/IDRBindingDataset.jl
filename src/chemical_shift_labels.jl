using DataFrames
using Statistics

# ---------------------------------------------------------------------------
# Residue labels from BMRB chemical shifts (free and bound states).
#
# Input: `disprot_bmrb_mapping.csv.gz` (script 21) and the cached BMRB entries.
# 1. Every shift list of the protein entity is renumbered on the UniProt
#    sequence (entity alignment, numbering and residue-type checks).
# 2. Free state: disorder is called from secondary shifts Δδ(Cα) − Δδ(Cβ) and
#    the amide 1H window (residues of IDRs have random-coil shifts).
# 3. Bound state: each bound list is compared with the most similar comparable
#    free list (pH, temperature, ionic strength, buffer, solvent, 15N window);
#    binding residues have a chemical-shift perturbation (CSP) above threshold
#    or a peak that disappears, in runs of consecutive residues.
# 4. Per residue: p_disorder = disordered calls / evaluable free lists,
#    p_binding = binding calls / pairs where the residue is observed.
#
# References:
# - secondary shifts / CSI: Wishart, Sykes & Richards 1992, Biochemistry 31:1647
# - Δδ(Cα) − Δδ(Cβ) (cancels 13C referencing errors): Metzler et al. 1993,
#   Biochemistry 32:13818 (PMID 8268157); for IDPs (SSP): Marsh et al. 2006,
#   Protein Sci 15:2795 (PMID 17088319)
# - random coil shifts: ncIDP, Tamiola, Acar & Mulder 2010, J Am Chem Soc
#   132:18000 (doi:10.1021/ja105656t; values from github.com/ktamiola/ncIDP,
#   private/original-data-csv/Tamiola.csv); Wishart et al. 1995, J Biomol NMR 5:67
# - CSP with 15N weight α = 0.14 (0.2 for Gly): Williamson 2013, Prog Nucl
#   Magn Reson Spectrosc 73:1
# Other thresholds are project choices, exposed as settings.
# ---------------------------------------------------------------------------

"""
    CHEMICAL_SHIFT_DEFAULT_SETTINGS

Default settings of `build_chemical_shift_labels`:
- `min_identity`, `max_indels`, `min_aligned`: BMRB entity / UniProt alignment
  required (1.0 = identical; with 1.0 a mutated or modified residue rejects the list)
- `strict_numbering`: reject a list whose residue numbering does not match the
  entity sequence on 100 % of atoms
- `trim_ends`: first and last residues of the construct not analysed (flexible ends)
- `random_coil`: `:ncidp` (Tamiola 2010) or `:wishart` (Wishart 1995)
- `secondary_shift_threshold`: |smoothed Δδ(Cα) − Δδ(Cβ)| below it (ppm) -> disordered
- `hn_min`, `hn_max`: amide 1H window of disordered residues (ppm)
- `min_run`, `min_run_hn`: minimum disordered run (Cα/Cβ available / 1H only)
- `max_delta_pH`, `max_delta_temperature` (K), `max_delta_ionic` (mM),
  `max_delta_field` (MHz), `buffer_tolerance` (relative), `window_tolerance`
  (relative), `n15_window_min` (ppm; wider 15N windows may differ): comparability
- `required_criteria`: criteria that must not differ (`:pH`, `:temperature`,
  `:ionic`, `:buffer`, `:solvent`, `:window`, `:field`, `:spectrometer`)
- `strict_metadata`: also reject pairs with unknown (not deposited) criteria
- `csp_k_sd`, `csp_min`: CSP threshold = mean + k·SD (after trimming), at least `csp_min`
- `min_class_n`: residues needed to compute a class-specific (disordered /
  ordered) CSP threshold
- `exclude_glycine_csp`, `use_lost_peaks`, `min_run_binding`,
  `binding_run_mode` (`:separate` / `:combined`): binding calls
- `min_common_residues`: residues with H/N in both states needed for a pair
- `p_disorder_threshold`: residue disordered when p_disorder ≥ this value
"""
const CHEMICAL_SHIFT_DEFAULT_SETTINGS = (
    min_identity = 1.0, max_indels = 0, min_aligned = 15, strict_numbering = true,
    trim_ends = 5, random_coil = :ncidp, secondary_shift_threshold = 1.0,
    hn_min = 7.8, hn_max = 8.6, min_run = 3, min_run_hn = 5,
    max_delta_pH = 0.5, max_delta_temperature = 5.0, max_delta_ionic = 10.0,
    max_delta_field = 1.0, buffer_tolerance = 0.0, window_tolerance = 0.01,
    n15_window_min = 30.0,
    required_criteria = (:pH, :temperature, :ionic, :buffer, :solvent, :window),
    strict_metadata = false, csp_k_sd = 1.5, csp_min = 0.05, min_class_n = 5,
    exclude_glycine_csp = false, use_lost_peaks = true, min_run_binding = 5,
    binding_run_mode = :separate, min_common_residues = 5, p_disorder_threshold = 0.5)

const CHEMICAL_SHIFT_CRITERIA = (:pH, :temperature, :ionic, :field, :spectrometer,
    :buffer, :solvent, :window)
const CSP_NITROGEN_WEIGHT = 0.14          # Williamson 2013
const CSP_NITROGEN_WEIGHT_GLYCINE = 0.20  # Williamson 2013

# ncIDP random coil shifts (Tamiola, Acar & Mulder 2010). For residue X:
# (previous, current, next) = correction when X is at i-1, value of X itself,
# correction when X is at i+1: δrc(i) = current(X_i) + previous(X_i-1) + next(X_i+1).
const NCIDP_CA = Dict(
    'A' => (0.161, 52.599, -0.187), 'C' => (-0.006, 58.327, 0.047),
    'E' => (0.171, 56.65, -0.071), 'D' => (0.29, 54.331, -0.095),
    'G' => (0.0, 45.236, 0.0), 'F' => (-0.163, 57.934, -0.247),
    'I' => (-0.087, 61.247, -0.129), 'H' => (0.069, 55.964, -0.248),
    'K' => (0.195, 56.412, -0.142), 'M' => (0.028, 55.591, -0.061),
    'L' => (0.083, 55.26, -0.184), 'N' => (0.293, 53.231, -0.067),
    'Q' => (0.125, 55.84, -0.098), 'P' => (0.021, 63.18, -2.423),
    'S' => (0.155, 58.352, -0.151), 'R' => (0.029, 56.088, -0.258),
    'T' => (0.143, 61.926, -0.122), 'W' => (0.014, 57.5, -0.201),
    'V' => (0.075, 62.347, -0.214), 'Y' => (-0.245, 57.761, -0.148))
const NCIDP_CB = Dict(
    'A' => (-0.135, 19.102, 0.018), 'C' => (0.008, 28.085, -0.17),
    'E' => (-0.119, 30.225, 0.065), 'D' => (-0.182, 41.089, 0.079),
    'G' => (0.0, 0.0, 0.0), 'F' => (0.059, 39.66, 0.023),
    'I' => (-0.111, 38.563, 0.093), 'H' => (-0.093, 29.719, 0.041),
    'K' => (-0.102, 32.921, 0.029), 'M' => (-0.139, 32.69, -0.205),
    'L' => (-0.143, 42.212, -0.053), 'N' => (-0.205, 38.79, -0.013),
    'Q' => (-0.064, 29.509, 0.004), 'P' => (-0.151, 32.072, -0.426),
    'S' => (-0.153, 63.766, 0.087), 'R' => (-0.038, 30.691, 0.046),
    'T' => (-0.072, 69.794, 0.086), 'W' => (-0.046, 29.38, 0.069),
    'V' => (-0.147, 32.674, 0.044), 'Y' => (0.017, 38.75, 0.019))
# Wishart et al. 1995 random coil (Cα, Cβ) without neighbour corrections.
const WISHART_RANDOM_COIL = Dict(
    'A' => (52.5, 19.1), 'R' => (56.0, 30.9), 'N' => (53.1, 38.9), 'D' => (54.2, 41.1),
    'C' => (58.2, 28.0), 'Q' => (55.7, 29.4), 'E' => (56.6, 29.9), 'G' => (45.1, NaN),
    'H' => (55.0, 29.0), 'I' => (61.1, 38.8), 'L' => (55.1, 42.4), 'K' => (56.2, 33.1),
    'M' => (55.4, 32.9), 'F' => (57.7, 39.6), 'P' => (63.3, 32.1), 'S' => (58.3, 63.8),
    'T' => (61.8, 69.8), 'W' => (57.5, 29.6), 'Y' => (57.9, 38.8), 'V' => (62.2, 32.9))

# ---------------------------------------------------------------------------
# Residue shift tables
# ---------------------------------------------------------------------------
"""
    _random_coil(sequence, position, residue, mode)

Random coil (Cα, Cβ) of `residue` at `position` of `sequence` (Cβ = `NaN` for
glycine); `:ncidp` adds the i-1 and i+1 neighbour corrections. Returns
`nothing` for unknown residue types.
"""
function _random_coil(sequence::AbstractString, position::Integer, residue::Char,
        mode::Symbol)
    mode === :wishart && return get(WISHART_RANDOM_COIL, residue, nothing)
    haskey(NCIDP_CA, residue) || return nothing
    characters = collect(sequence)
    ca = NCIDP_CA[residue][2]
    cb = residue == 'G' ? NaN : NCIDP_CB[residue][2]
    if position > 1 && haskey(NCIDP_CA, characters[position - 1])
        ca += NCIDP_CA[characters[position - 1]][1]
        residue == 'G' || (cb += NCIDP_CB[characters[position - 1]][1])
    end
    if position < length(characters) && haskey(NCIDP_CA, characters[position + 1])
        ca += NCIDP_CA[characters[position + 1]][3]
        residue == 'G' || (cb += NCIDP_CB[characters[position + 1]][3])
    end
    (ca, cb)
end

"""
    _residue_shift_table(table, target_entity_id, single_polymer, alignment,
        entity_sequence, protein_sequence, settings)

Turn `Atom_chem_shift` rows (one per atom) of the protein entity into one row per
UniProt position with `H`, `N`, `CA`, `CB` and `C` shifts. Returns
`(table, reason)`: `table` is `nothing` when the list is rejected, with
`reason` in `no_atom_columns`, `no_numbering_column`, `numbering_mismatch`,
`mutated_residues` or `no_backbone_shift`.
"""
function _residue_shift_table(table::AbstractDataFrame, target_entity_id::AbstractString,
        single_polymer::Bool, alignment, entity_sequence::AbstractString,
        protein_sequence::AbstractString, settings)
    columns = names(table)
    ("Atom_ID" in columns && "Val" in columns) || return (nothing, "no_atom_columns")
    has_entity = "Entity_ID" in columns
    has_residue_type = "Comp_ID" in columns
    keep(row) = !has_entity || _bmrb_value(row.Entity_ID) == target_entity_id ||
                (single_polymer && _bmrb_value(row.Entity_ID) === nothing)
    entity = collect(entity_sequence)
    protein = collect(protein_sequence)

    candidates = [c for c in ("Seq_ID", "Comp_index_ID") if c in columns]
    isempty(candidates) && return (nothing, "no_numbering_column")
    best_column, best_fraction = first(candidates), -1.0
    for column in candidates
        n_ok = n_atoms = 0
        for row in eachrow(table)
            (keep(row) && has_residue_type) || continue
            position = tryparse(Int, row[column])
            position === nothing && continue
            n_atoms += 1
            residue = _one_letter_code(row.Comp_ID)
            n_ok += 1 <= position <= length(entity) && entity[position] == residue
        end
        fraction = n_atoms > 0 ? n_ok / n_atoms : 0.0
        fraction > best_fraction && ((best_column, best_fraction) = (column, fraction))
    end
    if settings.strict_numbering && has_residue_type && best_fraction < 1.0
        return (nothing, "numbering_mismatch")
    end

    shifts = Dict{Int, Dict{String, Float64}}()
    residue_names = Dict{Int, String}()
    n_entity = length(entity)
    for row in eachrow(table)
        keep(row) || continue
        position = tryparse(Int, row[best_column])
        position === nothing && continue
        # Construct ends (expression tags included) are always flexible.
        (position <= settings.trim_ends || position > n_entity - settings.trim_ends) &&
            continue
        uniprot_position = get(alignment.mapping, position, 0)
        uniprot_position == 0 && continue   # tag, linker, outside the alignment
        atom = uppercase(strip(row.Atom_ID))
        atom = atom == "HN" ? "H" : atom
        atom in BMRB_BACKBONE_ATOMS || continue
        value = tryparse(Float64, row.Val)
        value === nothing && continue
        residue_shifts = get!(shifts, uniprot_position, Dict{String, Float64}())
        haskey(residue_shifts, atom) || (residue_shifts[atom] = value)
        residue_names[uniprot_position] = has_residue_type ? uppercase(row.Comp_ID) :
                                          string(protein[uniprot_position])
    end
    isempty(shifts) && return (nothing, "no_backbone_shift")

    positions = sort!(collect(keys(shifts)))
    result = DataFrame(pos = positions, aa = [protein[p] for p in positions],
        comp = [residue_names[p] for p in positions])
    for atom in BMRB_BACKBONE_ATOMS
        result[!, atom] = Union{Missing, Float64}[get(shifts[p], atom, missing)
                                                  for p in positions]
    end
    if has_residue_type
        n_mutated = count(i -> _one_letter_code(result.comp[i]) != result.aa[i],
            1:nrow(result))
        n_mutated > 0 && settings.min_identity >= 1.0 &&
            return (nothing, "mutated_residues")
    end
    (result, "ok")
end

# ---------------------------------------------------------------------------
# Disorder from free-state shifts
# ---------------------------------------------------------------------------
"""
    _secondary_shifts(residues, protein_sequence, mode)

Combined secondary shift Δδ(Cα) − Δδ(Cβ) per residue (raw and smoothed over
i-1, i, i+1). Glycines are excluded (no Cβ, very variable Cα).
"""
function _secondary_shifts(residues::AbstractDataFrame, protein_sequence::AbstractString,
        mode::Symbol)
    n = nrow(residues)
    raw = Vector{Union{Missing, Float64}}(missing, n)
    for i in 1:n
        residues.comp[i] == "GLY" && continue
        residue = _one_letter_code(residues.comp[i])
        residue == 'X' && (residue = residues.aa[i])
        coil = _random_coil(protein_sequence, residues.pos[i], residue, mode)
        coil === nothing && continue
        ca, cb = residues.CA[i], residues.CB[i]
        (ismissing(ca) || ismissing(cb) || isnan(coil[2])) && continue
        raw[i] = (ca - coil[1]) - (cb - coil[2])
    end
    by_position = Dict(residues.pos[i] => raw[i] for i in 1:n if !ismissing(raw[i]))
    smooth = Vector{Union{Missing, Float64}}(missing, n)
    for i in 1:n
        haskey(by_position, residues.pos[i]) || continue
        p = residues.pos[i]
        smooth[i] = mean(by_position[q]
        for q in (p - 1, p, p + 1) if haskey(by_position, q))
    end
    (raw, smooth)
end

"""
    _call_disorder!(residues, protein_sequence, settings)::String

Add `sec_smooth` and `disorder` (`true`/`false`/`missing`) to a free-state
residue table. With Cα/Cβ: |smoothed secondary shift| < threshold and amide 1H
(if known) in the disordered window; with 1H only: 1H in the window. Glycines
inherit the call of their neighbours when these agree; disordered runs shorter
than `min_run` (`min_run_hn` with 1H only) are set to ordered. Returns the
method used (`secondary_shifts` or `hn_window`).
"""
function _call_disorder!(residues::DataFrame, protein_sequence::AbstractString, settings)
    _, smooth = _secondary_shifts(residues, protein_sequence, settings.random_coil)
    residues.sec_smooth = smooth
    has_carbon = any(!ismissing, smooth)
    n = nrow(residues)
    call = Vector{Union{Missing, Bool}}(missing, n)
    for i in 1:n
        h = residues.H[i]
        hn_ok = ismissing(h) ? missing : settings.hn_min <= h <= settings.hn_max
        if has_carbon
            ismissing(smooth[i]) && continue
            call[i] = abs(smooth[i]) < settings.secondary_shift_threshold &&
                      (ismissing(hn_ok) || hn_ok)
        else
            ismissing(hn_ok) || (call[i] = hn_ok)
        end
    end
    if has_carbon
        index = Dict(residues.pos[i] => i for i in 1:n)
        for i in 1:n
            (residues.comp[i] == "GLY" && ismissing(call[i])) || continue
            neighbours = [call[index[q]]
                          for q in (residues.pos[i] - 1, residues.pos[i] + 1)
                          if haskey(index, q) && !ismissing(call[index[q]])]
            (!isempty(neighbours) && all(==(first(neighbours)), neighbours)) &&
                (call[i] = first(neighbours))
        end
    end
    min_run = has_carbon ? settings.min_run : settings.min_run_hn
    i = 1
    while i <= n
        if coalesce(call[i], false)
            j = i
            # A one-residue gap (e.g. proline without amide H) does not break a run.
            while j < n && coalesce(call[j + 1], false) &&
                  residues.pos[j + 1] - residues.pos[j] <= 2
                j += 1
            end
            j - i + 1 < min_run && (call[i:j] .= false)
            i = j + 1
        else
            i += 1
        end
    end
    residues.disorder = call
    has_carbon ? "secondary_shifts" : "hn_window"
end

# ---------------------------------------------------------------------------
# Comparability of free and bound experiments
# ---------------------------------------------------------------------------
"""
    _numeric_criterion(a, b, tolerance)::Symbol

`:same` within `tolerance`, `:different` otherwise, `:unknown` if a value is missing.
"""
function _numeric_criterion(a, b, tolerance::Real)::Symbol
    (ismissing(a) || ismissing(b)) && return :unknown
    abs(a - b) <= tolerance ? :same : :different
end

"""
    _compare_experiments(free, bound, settings)

Compare the metadata of a free and a bound shift list (`pH`, `temperature_K`,
`ionic_strength_mM`, `fields_MHz`, `spectrometers`, `buffer_mM`, `solvents`,
`n15_sweep_width_ppm`). Returns `(accepted, criteria, n_unknown)` where
`criteria` maps each of `CHEMICAL_SHIFT_CRITERIA` to `:same`, `:different` or
`:unknown`. Only `settings.required_criteria` decide `accepted`.
"""
function _compare_experiments(free, bound, settings)
    criteria = Dict{Symbol, Symbol}()
    criteria[:pH] = _numeric_criterion(free.pH, bound.pH, settings.max_delta_pH)
    criteria[:temperature] = _numeric_criterion(free.temperature_K, bound.temperature_K,
        settings.max_delta_temperature)
    criteria[:ionic] = _numeric_criterion(free.ionic_strength_mM, bound.ionic_strength_mM,
        settings.max_delta_ionic)
    criteria[:field] = (isempty(free.fields_MHz) || isempty(bound.fields_MHz)) ? :unknown :
                       any(
        x -> any(y -> abs(x - y) <= settings.max_delta_field,
            bound.fields_MHz),
        free.fields_MHz) ? :same : :different
    criteria[:spectrometer] = (isempty(free.spectrometers) ||
                               isempty(bound.spectrometers)) ?
                              :unknown :
                              isempty(intersect(lowercase.(free.spectrometers),
        lowercase.(bound.spectrometers))) ? :different : :same
    # Buffer: same components, concentrations within the relative tolerance.
    criteria[:buffer] = if isempty(free.buffer_mM) || isempty(bound.buffer_mM)
        :unknown
    elseif Set(keys(free.buffer_mM)) != Set(keys(bound.buffer_mM))
        :different
    else
        any(
            k -> !ismissing(free.buffer_mM[k]) && !ismissing(bound.buffer_mM[k]) &&
                 abs(free.buffer_mM[k] - bound.buffer_mM[k]) >
                 settings.buffer_tolerance *
                 max(free.buffer_mM[k], bound.buffer_mM[k]),
            collect(keys(free.buffer_mM))) ? :different : :same
    end
    # Solvent: both samples in H2O (e.g. 90 % H2O / 10 % D2O).
    criteria[:solvent] = (isempty(free.solvents) || isempty(bound.solvents)) ? :unknown :
                         all(s -> occursin("h2o", lowercase(s)),
        vcat(free.solvents, bound.solvents)) ? :same : :different
    # 15N window: a narrow window folds peaks (aliasing) and must then be identical.
    criteria[:window] = if ismissing(free.n15_sweep_width_ppm) ||
                           ismissing(bound.n15_sweep_width_ppm)
        :unknown
    else
        a, b = free.n15_sweep_width_ppm, bound.n15_sweep_width_ppm
        min(a, b) >= settings.n15_window_min ? :same :
        abs(a - b) <= settings.window_tolerance * max(a, b) ? :same : :different
    end
    required = [criteria[c] for c in settings.required_criteria]
    n_unknown = count(==(:unknown), required)
    accepted = !any(==(:different), required) &&
               (!settings.strict_metadata || n_unknown == 0)
    (accepted = accepted, criteria = criteria, n_unknown = n_unknown)
end

"""
    _condition_distance(free, bound)::Float64

Distance between sample conditions used to choose among comparable free lists
(10 × |ΔpH| + |ΔT|, 5 per unknown value).
"""
function _condition_distance(free, bound)::Float64
    d = (ismissing(free.pH) || ismissing(bound.pH)) ? 5.0 : 10 * abs(free.pH - bound.pH)
    d + ((ismissing(free.temperature_K) || ismissing(bound.temperature_K)) ? 5.0 :
     abs(free.temperature_K - bound.temperature_K))
end

# ---------------------------------------------------------------------------
# Chemical-shift perturbations
# ---------------------------------------------------------------------------
"""
    _keep_runs(flags, positions, min_run)::BitVector

Keep only `true` runs spanning at least `min_run` positions; a gap of one
position (residue without peak, e.g. proline) does not break a run.
"""
function _keep_runs(flags::AbstractVector{Bool}, positions::AbstractVector{<:Integer},
        min_run::Integer)::BitVector
    kept = falses(length(flags))
    i, n = 1, length(flags)
    while i <= n
        if flags[i]
            j = i
            while j < n && flags[j + 1] && positions[j + 1] - positions[j] <= 2
                j += 1
            end
            positions[j] - positions[i] + 1 >= min_run && (kept[i:j] .= true)
            i = j + 1
        else
            i += 1
        end
    end
    kept
end

"""
    _csp_threshold(values, settings)::Float64

Mean + k·SD of `values` after iteratively removing values above mean + 3 SD,
with a floor of `settings.csp_min`.
"""
function _csp_threshold(values, settings)::Float64
    kept = Float64.(collect(values))
    for _ in 1:20
        length(kept) < 3 && break
        mu, sigma = mean(kept), std(kept)
        trimmed = kept[kept .<= mu + 3sigma]
        length(trimmed) == length(kept) && break
        kept = trimmed
    end
    length(kept) < 2 && return settings.csp_min
    max(mean(kept) + settings.csp_k_sd * std(kept), settings.csp_min)
end

"""
    _chemical_shift_pair(free, bound, settings)

Compare one free and one bound residue table (`pos`, `comp`, `H`, `N`; the free
one also has `disorder`). Global referencing offsets (median ΔH, ΔN) are
removed; CSP = √[(ΔH² + (α·ΔN)²)/2]. Thresholds are computed separately on
residues disordered / ordered in the free state (global threshold for classes
smaller than `min_class_n`). A peak present in the free state and absent in
the bound state, inside the assigned region, is `lost`. Binding residues are
runs of at least `min_run_binding` perturbed (and/or lost) residues.

Returns `nothing` when fewer than `min_common_residues` residues have H/N in
both states, otherwise `(table, threshold_disorder, threshold_order, n_common)`.
"""
function _chemical_shift_pair(free::AbstractDataFrame, bound::AbstractDataFrame, settings)
    bound_hn = Dict(bound.pos[i] => (bound.H[i], bound.N[i])
    for i in 1:nrow(bound) if !ismissing(bound.H[i]) && !ismissing(bound.N[i]))
    isempty(bound_hn) && return nothing
    pair = DataFrame(pos = Int[], comp = String[], H_free = Float64[], N_free = Float64[],
        H_bound = Union{Missing, Float64}[], N_bound = Union{Missing, Float64}[],
        disorder_free = Union{Missing, Bool}[])
    for i in 1:nrow(free)
        (ismissing(free.H[i]) || ismissing(free.N[i])) && continue
        h, n = get(bound_hn, free.pos[i], (missing, missing))
        push!(pair, (free.pos[i], free.comp[i], free.H[i], free.N[i], h, n,
            free.disorder[i]))
    end
    common = .!ismissing.(pair.H_bound)
    n_common = count(common)
    n_common < settings.min_common_residues && return nothing

    offset_h = median(Float64.(pair.H_bound[common]) .- pair.H_free[common])
    offset_n = median(Float64.(pair.N_bound[common]) .- pair.N_free[common])
    pair.csp = Union{Missing, Float64}[common[i] && !(settings.exclude_glycine_csp &&
                                         pair.comp[i] == "GLY") ?
                                       sqrt(((pair.H_bound[i] - pair.H_free[i] - offset_h)^2 +
                                             ((pair.comp[i] == "GLY" ?
                                               CSP_NITROGEN_WEIGHT_GLYCINE :
                                               CSP_NITROGEN_WEIGHT) *
                                              (pair.N_bound[i] - pair.N_free[i] -
                                               offset_n))^2) / 2) : missing
                                       for i in 1:nrow(pair)]
    has_csp = .!ismissing.(pair.csp)
    count(has_csp) < settings.min_common_residues && return nothing
    disordered = coalesce.(pair.disorder_free .== true, false)
    ordered = coalesce.(pair.disorder_free .== false, false)
    threshold_all = _csp_threshold(pair.csp[has_csp], settings)
    threshold_disorder = count(has_csp .& disordered) >= settings.min_class_n ?
                         _csp_threshold(pair.csp[has_csp .& disordered], settings) :
                         threshold_all
    threshold_order = count(has_csp .& ordered) >= settings.min_class_n ?
                      _csp_threshold(pair.csp[has_csp .& ordered], settings) :
                      threshold_all
    thresholds = [disordered[i] ? threshold_disorder :
                  ordered[i] ? threshold_order :
                  threshold_all for i in 1:nrow(pair)]
    perturbed = [has_csp[i] && pair.csp[i] > thresholds[i] for i in 1:nrow(pair)]
    low, high = extrema(keys(bound_hn))
    lost = [settings.use_lost_peaks && !common[i] && low <= pair.pos[i] <= high
            for i in 1:nrow(pair)]
    pair.binding = if settings.binding_run_mode === :combined
        _keep_runs(perturbed .| lost, pair.pos, settings.min_run_binding)
    else
        _keep_runs(perturbed, pair.pos, settings.min_run_binding) .|
        _keep_runs(lost, pair.pos, settings.min_run_binding)
    end
    pair.lost = lost
    (table = pair, threshold_disorder = threshold_disorder,
        threshold_order = threshold_order, n_common = n_common)
end

# ---------------------------------------------------------------------------
# Per-protein computation
# ---------------------------------------------------------------------------
"""
    _empty_chemical_shift_tables()

Empty output tables (schemas) of `build_chemical_shift_labels`.
"""
function _empty_chemical_shift_tables()
    labels = DataFrame(protein_disprot_id = String[], protein_id = String[],
        n_free_lists = Int[], n_pairs = Int[], free_lists = String[],
        pairs = Union{Missing, String}[], disprot = String[], cs_n_free = String[],
        cs_n_pairs = String[], Label_chemical_shift = String[],
        Label_chemical_shift_disorder = String[], Label_chemical_shift_binding = String[])
    lists = DataFrame(protein_disprot_id = String[], protein_id = String[],
        bmrb_id = String[], shift_list_id = Union{Missing, String}[], state = String[],
        used = Bool[], reason = String[], n_residues = Union{Missing, Int}[],
        disorder_method = Union{Missing, String}[], n_disordered = Union{Missing, Int}[],
        n_binding = Union{Missing, Int}[])
    # `Pair{Symbol, Any}` keeps each column type (a plain `vcat` would promote them).
    columns = Pair{Symbol, Any}[:protein_disprot_id => String[], :protein_id => String[],
        :bound_bmrb_id => String[], :bound_shift_list_id => String[],
        :free_bmrb_id => String[], :free_shift_list_id => String[]]
    append!(columns, [c => String[] for c in CHEMICAL_SHIFT_CRITERIA])
    append!(columns,
        Pair{Symbol, Any}[:accepted => Bool[], :selected => Bool[],
            :n_common => Union{Missing, Int}[], :n_binding => Union{Missing, Int}[],
            :csp_threshold_disorder => Union{Missing, Float64}[],
            :csp_threshold_order => Union{Missing, Float64}[]])
    comparisons = DataFrame(columns)
    (labels = labels, lists = lists, comparisons = comparisons)
end

"""
    _shift_list_experiment(row, json, protein_sequence, settings)

Rebuild one usable shift list (`row` of the BMRB mapping table) from its parsed
entry `json`: residue table and comparability metadata. Returns
`(experiment, reason)` with `experiment === nothing` when rejected.
"""
function _shift_list_experiment(row, json, protein_sequence::AbstractString, settings)
    saveframes = parse_bmrb_entry(json)
    entities = _bmrb_entities(saveframes)
    target = _bmrb_target_entity(entities, protein_sequence)
    target === nothing && return (nothing, "no_polypeptide")
    target_index, alignment = target
    entity = entities[target_index, :]
    shift_list = findfirst(
        s -> s.category == "assigned_chemical_shifts" &&
             _bmrb_tag(s, "ID") == string(row.shift_list_id),
        saveframes)
    shift_list === nothing && return (nothing, "shift_list_not_found")
    shift_list = saveframes[shift_list]
    atom_tables = [loop.table
                   for loop in shift_list.loops if loop.category == "Atom_chem_shift"]
    isempty(atom_tables) && return (nothing, "no_atom_columns")
    single_polymer = count(==("polymer"), entities.type) == 1
    residues, reason = _residue_shift_table(vcat(atom_tables...; cols = :union),
        entity.entity_id, single_polymer, alignment, entity.sequence, protein_sequence,
        settings)
    residues === nothing && return (nothing, reason)

    condition_id = _bmrb_tag(shift_list, "Sample_condition_list_ID")
    condition = get(_bmrb_sample_conditions(saveframes), something(condition_id, "?"),
        (pH = missing, temperature_K = missing, ionic_strength_mM = missing))
    metadata = _bmrb_shift_list_metadata(shift_list,
        _bmrb_entry_metadata(saveframes, entities))
    experiment = (bmrb_id = string(row.bmrb_id), shift_list_id = string(row.shift_list_id),
        state = String(row.state), residues = residues, pH = condition.pH,
        temperature_K = condition.temperature_K,
        ionic_strength_mM = condition.ionic_strength_mM, fields_MHz = metadata.fields_MHz,
        spectrometers = metadata.spectrometers, buffer_mM = metadata.buffer_mM,
        solvents = metadata.solvents, n15_sweep_width_ppm = metadata.n15_sweep_width_ppm)
    (experiment, "ok")
end

"""
    _mapping_rejection_reason(row, settings)

Reason why a BMRB mapping row cannot be used before reading its shifts, or
`nothing` when it can.
"""
function _mapping_rejection_reason(row, settings)
    coalesce(row.status, "") == "ok" || return "status:$(coalesce(row.status, "missing"))"
    state = coalesce(row.state, "")
    state == "homo_oligomer" && return "homo_oligomer"
    state in ("free", "bound") || return "state:$(state)"
    coalesce(row.entity_identity, 0.0) < settings.min_identity &&
        return "low_entity_identity"
    coalesce(row.entity_n_aligned, 0) < settings.min_aligned &&
        return "few_aligned_residues"
    coalesce(row.entity_indels, typemax(Int)) > settings.max_indels && return "indels"
    nothing
end

"""
    _protein_chemical_shift_tables(protein_disprot_id, protein_id, sequence,
        residue_labels, mapping_rows, entry_json, settings)

Labels, shift-list log and free/bound comparisons of one protein. `entry_json`
returns the parsed BMRB entry for an entry ID.
"""
function _protein_chemical_shift_tables(protein_disprot_id::AbstractString,
        protein_id::AbstractString, sequence::AbstractString,
        residue_labels::AbstractVector{<:Integer}, mapping_rows::AbstractDataFrame,
        entry_json::Function, settings)
    tables = _empty_chemical_shift_tables()
    ids = (String(protein_disprot_id), String(protein_id))
    log_list(row, used, reason; n_residues = missing, method = missing,
        n_disordered = missing, n_binding = missing) = push!(tables.lists,
        (ids..., string(row.bmrb_id),
            ismissing(row.shift_list_id) ? missing :
            string(row.shift_list_id),
            String(coalesce(row.state, "missing")), used, reason, n_residues, method,
            n_disordered, n_binding))

    free, bound = NamedTuple[], NamedTuple[]
    for row in eachrow(mapping_rows)
        reason = _mapping_rejection_reason(row, settings)
        if reason !== nothing
            log_list(row, false, reason)
            continue
        end
        experiment, reason = _shift_list_experiment(row, entry_json(string(row.bmrb_id)),
            sequence, settings)
        if experiment === nothing
            log_list(row, false, reason)
            continue
        end
        push!(experiment.state == "free" ? free : bound, (
            row = row, experiment = experiment))
    end

    # Free state: disorder calls.
    free_used = NamedTuple[]
    for item in free
        method = _call_disorder!(item.experiment.residues, sequence, settings)
        calls = item.experiment.residues.disorder
        if all(ismissing, calls)
            log_list(item.row, false, "no_evaluable_residue";
                n_residues = nrow(item.experiment.residues), method = method)
            continue
        end
        log_list(item.row, true, "used"; n_residues = nrow(item.experiment.residues),
            method = method, n_disordered = count(x -> coalesce(x, false), calls))
        push!(free_used, item)
    end

    # Bound state: comparable free list -> CSP -> binding.
    cs_pairs = NamedTuple[]
    for item in bound
        b = item.experiment
        if isempty(free_used)
            log_list(item.row, false, "no_free_state"; n_residues = nrow(b.residues))
            continue
        end
        candidates = []
        comparison_rows = Int[]
        for free_item in free_used
            f = free_item.experiment
            comparison = _compare_experiments(f, b, settings)
            push!(tables.comparisons,
                (ids..., b.bmrb_id, b.shift_list_id, f.bmrb_id, f.shift_list_id,
                    [string(comparison.criteria[c]) for c in CHEMICAL_SHIFT_CRITERIA]...,
                    comparison.accepted, false, missing, missing, missing, missing))
            comparison.accepted || continue
            push!(candidates, (free = free_item, comparison = comparison))
            push!(comparison_rows, nrow(tables.comparisons))
        end
        if isempty(candidates)
            log_list(item.row, false, "no_comparable_free"; n_residues = nrow(b.residues))
            continue
        end
        # Fewest unknown criteria, then closest conditions, then most residues.
        best = argmin([c.comparison.n_unknown * 100 +
                       _condition_distance(c.free.experiment, b) -
                       1e-4 * nrow(c.free.experiment.residues) for c in candidates])
        f = candidates[best].free.experiment
        row_index = comparison_rows[best]
        tables.comparisons.selected[row_index] = true
        result = _chemical_shift_pair(f.residues, b.residues, settings)
        if result === nothing
            log_list(item.row, false, "few_common_residues"; n_residues = nrow(b.residues))
            continue
        end
        n_binding = count(result.table.binding)
        tables.comparisons.n_common[row_index] = result.n_common
        tables.comparisons.n_binding[row_index] = n_binding
        tables.comparisons.csp_threshold_disorder[row_index] = result.threshold_disorder
        tables.comparisons.csp_threshold_order[row_index] = result.threshold_order
        log_list(
            item.row, true, "used"; n_residues = nrow(b.residues), n_binding = n_binding)
        push!(cs_pairs, (free = f, bound = b, result = result))
    end

    isempty(free_used) && return tables

    # Per-residue aggregation.
    n_positions = length(sequence)
    n_disordered, n_ordered = zeros(Int, n_positions), zeros(Int, n_positions)
    for item in free_used, i in 1:nrow(item.experiment.residues)

        call = item.experiment.residues.disorder[i]
        ismissing(call) && continue
        position = item.experiment.residues.pos[i]
        call ? (n_disordered[position] += 1) : (n_ordered[position] += 1)
    end
    n_pairs, n_pairs_binding = zeros(Int, n_positions), zeros(Int, n_positions)
    for pair in cs_pairs, r in eachrow(pair.result.table)
        # A pair is evidence only where the residue is observed (CSP or lost peak).
        (!ismissing(r.csp) || r.lost) || continue
        n_pairs[r.pos] += 1
        r.binding && (n_pairs_binding[r.pos] += 1)
    end
    n_free = n_disordered .+ n_ordered
    p_disorder = [n > 0 ? d / n : NaN for (d, n) in zip(n_disordered, n_free)]
    p_binding = [n > 0 ? b / n : NaN for (b, n) in zip(n_pairs_binding, n_pairs)]
    disordered = [!isnan(p) && p >= settings.p_disorder_threshold for p in p_disorder]
    label_binding = [disordered[p] && n_pairs[p] > 0 ? p_binding[p] : NaN
                     for p in 1:n_positions]
    label = [n_free[p] == 0 ? NaN :
             disordered[p] ? 1.0 + (isnan(label_binding[p]) ? 0.0 : label_binding[p]) : 0.0
             for p in 1:n_positions]

    list_name(e) = "$(e.bmrb_id):$(e.shift_list_id)"
    push!(tables.labels,
        (ids..., length(free_used), length(cs_pairs),
            join([list_name(item.experiment) for item in free_used], ";"),
            isempty(cs_pairs) ? missing :
            join(["$(list_name(p.free))>$(list_name(p.bound))" for p in cs_pairs], ";"),
            _format_label_vector(Float64.(residue_labels)),
            _format_label_vector(Float64.(n_free)), _format_label_vector(Float64.(n_pairs)),
            _format_label_vector(label), _format_label_vector(p_disorder),
            _format_label_vector(label_binding)))
    tables
end

"""
    build_chemical_shift_labels(dataset, bmrb_mapping; cache_dir,
        settings = CHEMICAL_SHIFT_DEFAULT_SETTINGS, get_request = _default_bmrb_get_request,
        logger = println)

Build chemical-shift residue labels for the dataset proteins of `bmrb_mapping`
(output of `build_disprot_bmrb_mapping_table`), reading the BMRB entries from
`cache_dir` (downloaded when missing).

Returns `(labels, lists, comparisons, failed_proteins)`:
- `labels`: one row per DisProt protein with at least one usable free list;
  vectors have one value per UniProt residue (`NaN` = no NMR data):
  - `Label_chemical_shift`: `0` ordered, `1 + p_binding` disordered
    (`p_disorder ≥ p_disorder_threshold`; `1` when no pair covers the residue)
  - `Label_chemical_shift_disorder`: `p_disorder`, fraction of free lists calling the residue disordered
  - `Label_chemical_shift_binding`: `p_binding` on disordered residues covered by a pair
- `lists`: every BMRB shift list of the protein, with `used` and `reason`
- `comparisons`: every bound × usable free list comparison, with each criterion
  (`same` / `different` / `unknown`), `accepted` and `selected`
"""
function build_chemical_shift_labels(dataset::AbstractDataFrame,
        bmrb_mapping::AbstractDataFrame; cache_dir::AbstractString,
        settings = CHEMICAL_SHIFT_DEFAULT_SETTINGS,
        get_request::Function = _default_bmrb_get_request, logger::Function = println)
    settings.random_coil in (:ncidp, :wishart) ||
        throw(ArgumentError("`random_coil` must be :ncidp or :wishart."))
    settings.binding_run_mode in (:separate, :combined) ||
        throw(ArgumentError("`binding_run_mode` must be :separate or :combined."))
    all(in(CHEMICAL_SHIFT_CRITERIA), settings.required_criteria) ||
        throw(ArgumentError("`required_criteria` must be among $(CHEMICAL_SHIFT_CRITERIA)."))
    _require_columns(
        dataset, [:protein_id, :protein_disprot_id, :sequence, :residue_labels],
        "dataset")
    _require_columns(bmrb_mapping,
        [:protein_disprot_id, :protein_id, :bmrb_id, :status, :state, :shift_list_id,
            :entity_identity, :entity_n_aligned, :entity_indels], "bmrb_mapping")

    labels_by_protein = _dataset_residue_labels_by_protein(dataset)
    sequences = Dict(String(r.protein_id) => uppercase(string(r.sequence))
    for r in eachrow(dataset) if !ismissing(r.sequence))
    tables = _empty_chemical_shift_tables()
    failed = String[]
    entry_cache = Dict{String, Any}()
    entry_json(bmrb_id) = get!(entry_cache, bmrb_id) do
        _cached_bmrb_json(joinpath(cache_dir, "entry_$(bmrb_id).json"),
            "$(BMRB_API_URL)/entry/$(bmrb_id)"; get_request = get_request)
    end

    groups = groupby(DataFrame(bmrb_mapping), :protein_id; sort = true)
    logger("Proteins with BMRB entries: $(length(groups))")
    for (index, rows) in enumerate(groups)
        protein_id = String(first(rows.protein_id))
        (haskey(sequences, protein_id) && haskey(labels_by_protein, protein_id)) || continue
        empty!(entry_cache)
        try
            protein = _protein_chemical_shift_tables(
                String(first(rows.protein_disprot_id)),
                protein_id, sequences[protein_id], labels_by_protein[protein_id], rows,
                entry_json, settings)
            for name in keys(tables)
                append!(tables[name], protein[name])
            end
        catch err
            push!(failed, protein_id)
            logger("Chemical shift labels failed for $(protein_id): $(sprint(showerror, err))")
        end
        index % 50 == 0 && logger("  $(index)/$(length(groups)) proteins done")
    end
    sort!(tables.lists, [:protein_disprot_id, :protein_id, :bmrb_id, :shift_list_id])
    sort!(tables.comparisons,
        [:protein_disprot_id, :protein_id, :bound_bmrb_id,
            :bound_shift_list_id, :free_bmrb_id, :free_shift_list_id])
    (labels = tables.labels, lists = tables.lists, comparisons = tables.comparisons,
        failed_proteins = failed)
end

"""
    _CHEMICAL_SHIFT_LABELS_COLUMN_DEFINITIONS

Data dictionary of the `labels` table of `build_chemical_shift_labels`.
"""
const _CHEMICAL_SHIFT_LABELS_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "n_free_lists" => "free-state shift lists used for disorder",
    "n_pairs" => "free/bound pairs used for binding",
    "free_lists" => "free lists used (`bmrb_id:shift_list_id`, `;`-separated)",
    "pairs" => "free/bound pairs used (`free>bound`)",
    "disprot" => "DisProt residue labels (`0`/`1`/`2`), for comparison",
    "cs_n_free" => "free lists where the residue disorder could be evaluated",
    "cs_n_pairs" => "pairs where the residue is observed (CSP or lost peak)",
    "Label_chemical_shift" => "`0` ordered, `1 + p_binding` disordered (`p_disorder ≥ p_disorder_threshold`; `1` when no pair covers the residue); `NaN` without NMR data",
    "Label_chemical_shift_disorder" => "`p_disorder`: fraction of free lists calling the residue disordered (secondary shifts and amide 1H window)",
    "Label_chemical_shift_binding" => "`p_binding`: fraction of pairs calling the residue binding (CSP or lost peak runs), on NMR-disordered residues covered by a pair; `NaN` otherwise"
]
