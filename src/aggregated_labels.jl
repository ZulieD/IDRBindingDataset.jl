using DataFrames

# ---------------------------------------------------------------------------
# Residue-level binding labels aggregated per DisProt protein.
#
# Combines three sources for each dataset protein:
# - DisProt residue labels (`residue_labels` of the final dataset),
# - UniProt binding sites (`uniprot_binding_sites.csv.gz`, script 18),
# - heteromeric PDB contacts (`pdb_contacts.csv.gz`, script 19).
#
# PDB evidence: N = number of used PDB entries where a residue is observed,
# x = number of those where it contacts a partner chain (contacts are only
# kept on DisProt disordered residues, labels 1 or 2), p_pdb = x / N.
# Final rate per residue = (DisProt 2 + UniProt 2 + PDB vote) / number of
# methods covering the residue.
#
# Thresholds are project choices (not taken from a reference), kept as named
# constants and exposed as keyword arguments.
# ---------------------------------------------------------------------------
const PDB_LABEL_MAX_RESOLUTION = 5.0      # PDB used if resolution < this (Å)
const PDB_LABEL_VOTES = (:fraction, :any, :majority)

"""
    _passes_pdb_quality(category, resolution; max_resolution, keep_no_resolution, keep_no_info)

Return `true` when one PDB entry passes the quality filter. `category` is the
`quality_category` of the contacts table (`resolved`, `no_resolution`, `absent`).
"""
function _passes_pdb_quality(
        category::AbstractString,
        resolution;
        max_resolution::Real = PDB_LABEL_MAX_RESOLUTION,
        keep_no_resolution::Bool = true,
        keep_no_info::Bool = true
)::Bool
    category == "resolved" && return !ismissing(resolution) && resolution < max_resolution
    category == "no_resolution" ? keep_no_resolution : keep_no_info
end

"""
    _pdb_vote(p_pdb, vote)::Float64

Contribution of the aggregated PDB evidence to the numerator of the final rate:
`:fraction` -> `p_pdb`; `:any` -> `1` if `p_pdb > 0`; `:majority` -> `1` if
`p_pdb >= 0.5`.
"""
function _pdb_vote(p_pdb::Real, vote::Symbol)::Float64
    vote === :fraction && return Float64(p_pdb)
    vote === :any && return Float64(p_pdb > 0)
    vote === :majority && return Float64(p_pdb >= 0.5)
    throw(ArgumentError("Unknown PDB vote $(repr(vote)); use one of $(PDB_LABEL_VOTES)."))
end

"""
    _final_binding_rates(disprot, uniprot, has_uniprot, pdb_n, pdb_p; vote = :fraction)

Combine methods per residue. For each position, `N` counts the methods covering
it (DisProt if labelled, UniProt if the protein has UniProt annotations, PDB if
observed in at least one PDB) and `x` the binding evidence (DisProt `2`, UniProt
`2`, PDB vote). Returns `(n_methods, x_methods, rate, final_label)` where
`final_label` is `1 + rate` when `rate > 0`, else the DisProt label clamped to
`[0, 1]` (`NaN` when uncovered).
"""
function _final_binding_rates(
        disprot::AbstractVector{Float64},
        uniprot::AbstractVector{Float64},
        has_uniprot::Bool,
        pdb_n::AbstractVector{<:Integer},
        pdb_p::AbstractVector{Float64};
        vote::Symbol = :fraction
)
    n_positions = length(disprot)
    n_methods = zeros(Int, n_positions)
    x_methods = zeros(n_positions)
    rate = fill(NaN, n_positions)
    final_label = fill(NaN, n_positions)
    for position in 1:n_positions
        n, x = 0, 0.0
        if !isnan(disprot[position])
            n += 1
            x += disprot[position] == 2.0
        end
        if has_uniprot
            n += 1
            x += uniprot[position] == 2.0
        end
        if pdb_n[position] > 0
            n += 1
            x += _pdb_vote(pdb_p[position], vote)
        end
        n_methods[position], x_methods[position] = n, x
        n == 0 && continue
        rate[position] = x / n
        final_label[position] = rate[position] > 0 ? 1.0 + rate[position] :
                                isnan(disprot[position]) ? 0.0 : min(disprot[position], 1.0)
    end
    (n_methods = n_methods, x_methods = x_methods, rate = rate, final_label = final_label)
end

"""
    _source_labels(disprot, in_idr, covered, score)

Per-source labels where disorder comes from DisProt (`in_idr`) restricted to
residues `covered` by the source and labelled in DisProt (`NaN` elsewhere).
`score[p]` in `[0, 1]` is the binding strength given by the source.

Returns `(label, disorder, binding)`:
- `label`: `0` ordered, `1 + score` disordered
- `disorder`: `0` / `1`
- `binding`: `score` on disordered residues, `NaN` otherwise
"""
function _source_labels(
        disprot::AbstractVector{Float64},
        in_idr::AbstractVector{Bool},
        covered::AbstractVector{Bool},
        score::AbstractVector{Float64}
)
    n_positions = length(disprot)
    label, disorder, binding = fill(NaN, n_positions), fill(NaN, n_positions),
    fill(NaN, n_positions)
    for position in 1:n_positions
        (covered[position] && !isnan(disprot[position])) || continue
        disorder[position] = Float64(in_idr[position])
        if in_idr[position]
            binding[position] = score[position]
            label[position] = 1.0 + score[position]
        else
            label[position] = 0.0
        end
    end
    (label = label, disorder = disorder, binding = binding)
end

"""
    _format_label_vector(values)::String

Encode a numeric vector as `"[0.0,1.5,NaN]"` (3 decimal digits).
"""
function _format_label_vector(values::AbstractVector{<:Real})::String
    "[" *
    join((isnan(v) ? "NaN" : string(round(Float64(v); digits = 3)) for v in values),
        ",") * "]"
end

"""
    _parse_label_vector(text)::Vector{Float64}

Parse a vector written by `_format_label_vector` (`"[0.0,1.5,NaN]"`).
"""
function _parse_label_vector(text::AbstractString)::Vector{Float64}
    inner = strip(strip(text), ['[', ']'])
    isempty(inner) && return Float64[]
    [parse(Float64, strip(token)) for token in split(inner, ",")]
end

"""
    _empty_aggregated_labels_table()

Return the empty output table (schema) of `build_aggregated_labels_table`.
"""
function _empty_aggregated_labels_table()::DataFrame
    DataFrame(protein_disprot_id = String[], protein_id = String[], n_pdb = Int[],
        pdb_ids = String[], has_uniprot = Bool[], disprot = String[], uniprot = String[],
        pdb_N = String[], pdb_proba = String[], rate = String[], final_label = String[],
        order_disorder = String[], disorder_binding = String[], Label_PDB = String[],
        Label_PDB_disorder = String[], Label_PDB_binding = String[],
        Label_Uniprot = String[], Label_Uniprot_disorder = String[],
        Label_Uniprot_binding = String[])
end

"""
    _protein_aggregated_labels_row(protein_disprot_id, protein_id, residue_labels,
        uniprot_positions, pdb_rows; vote)

Aggregate one protein. `pdb_rows` are the contacts-table rows of the PDB
entries kept for this protein. Returns `nothing` when no PDB entry has an
observed residue.
"""
function _protein_aggregated_labels_row(
        protein_disprot_id::AbstractString,
        protein_id::AbstractString,
        residue_labels::AbstractVector{<:Integer},
        uniprot_positions::AbstractVector{<:Integer},
        pdb_rows;
        vote::Symbol
)
    n_positions = length(residue_labels)
    disprot = Float64.(residue_labels)
    in_idr = [label in (_RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING)
              for label in residue_labels]

    pdb_n, pdb_x = zeros(Int, n_positions), zeros(Int, n_positions)
    used_pdbs = String[]
    for row in pdb_rows
        resolved = filter(p -> 1 <= p <= n_positions, _parse_residue_ranges(row.resolved_positions))
        isempty(resolved) && continue
        contacts = Set(p
        for p in _parse_residue_ranges(row.contact_positions)
        if 1 <= p <= n_positions && in_idr[p])
        for position in resolved
            pdb_n[position] += 1
            position in contacts && (pdb_x[position] += 1)
        end
        push!(used_pdbs, uppercase(row.PDB))
    end
    isempty(used_pdbs) && return nothing

    has_uniprot = !isempty(uniprot_positions)
    uniprot = has_uniprot ? zeros(n_positions) : fill(NaN, n_positions)
    for position in uniprot_positions
        1 <= position <= n_positions && (uniprot[position] = 2.0)
    end

    pdb_p = [n > 0 ? x / n : NaN for (n, x) in zip(pdb_n, pdb_x)]
    rates = _final_binding_rates(disprot, uniprot, has_uniprot, pdb_n, pdb_p; vote = vote)
    order_disorder = Float64.(in_idr)
    disorder_binding = [in_idr[p] ? rates.rate[p] : NaN for p in 1:n_positions]
    pdb_labels = _source_labels(disprot, in_idr, pdb_n .> 0,
        [isnan(p) ? 0.0 : p for p in pdb_p])
    uniprot_labels = _source_labels(disprot, in_idr, fill(has_uniprot, n_positions),
        [Float64(u == 2.0) for u in uniprot])

    (protein_disprot_id = String(protein_disprot_id), protein_id = String(protein_id),
        n_pdb = length(used_pdbs), pdb_ids = join(sort!(used_pdbs), ";"),
        has_uniprot = has_uniprot, disprot = _format_label_vector(disprot),
        uniprot = _format_label_vector(uniprot),
        pdb_N = _format_label_vector(Float64.(pdb_n)), pdb_proba = _format_label_vector(pdb_p),
        rate = _format_label_vector(rates.rate),
        final_label = _format_label_vector(rates.final_label),
        order_disorder = _format_label_vector(order_disorder),
        disorder_binding = _format_label_vector(disorder_binding),
        Label_PDB = _format_label_vector(pdb_labels.label),
        Label_PDB_disorder = _format_label_vector(pdb_labels.disorder),
        Label_PDB_binding = _format_label_vector(pdb_labels.binding),
        Label_Uniprot = _format_label_vector(uniprot_labels.label),
        Label_Uniprot_disorder = _format_label_vector(uniprot_labels.disorder),
        Label_Uniprot_binding = _format_label_vector(uniprot_labels.binding))
end

"""
    build_aggregated_labels_table(dataset, uniprot_binding_sites, pdb_contacts;
        max_resolution = PDB_LABEL_MAX_RESOLUTION, keep_no_resolution = true,
        keep_no_info = true, keep_pdb_without_partner = false, vote = :fraction)

Aggregate DisProt labels, UniProt binding sites (`build_uniprot_binding_sites_table`)
and PDB contacts (`build_pdb_contacts_table`) into one row per DisProt protein.

PDB entries are used when their `status` is `ok` (or `no_heteromeric_partner`
with `keep_pdb_without_partner = true`, counting observed residues without
contact) and they pass the quality filter (`resolution < max_resolution`;
entries without resolution or absent from the quality table are kept unless
`keep_no_resolution` / `keep_no_info` are `false`). PDB contacts are only
counted on DisProt disordered residues (labels `1`/`2`).

Only proteins with at least one used PDB entry are returned. Returns
`(table, pdb_log)` where `pdb_log` gives, for every contacts row, whether the
entry was used and why (`used`, `status`, `low_quality`, `no_resolved_position`).
"""
function build_aggregated_labels_table(
        dataset::AbstractDataFrame,
        uniprot_binding_sites::AbstractDataFrame,
        pdb_contacts::AbstractDataFrame;
        max_resolution::Real = PDB_LABEL_MAX_RESOLUTION,
        keep_no_resolution::Bool = true,
        keep_no_info::Bool = true,
        keep_pdb_without_partner::Bool = false,
        vote::Symbol = :fraction
)
    vote in PDB_LABEL_VOTES ||
        throw(ArgumentError("Unknown PDB vote $(repr(vote)); use one of $(PDB_LABEL_VOTES)."))
    _require_columns(dataset, [:protein_id, :protein_disprot_id, :residue_labels], "dataset")
    _require_columns(uniprot_binding_sites, [:protein_id, :residues_binding_in_uniprot],
        "uniprot_binding_sites")
    _require_columns(pdb_contacts,
        [:protein_disprot_id, :protein_id, :PDB, :status, :quality_category, :resolution,
            :resolved_positions, :contact_positions], "pdb_contacts")

    labels_by_protein = _dataset_residue_labels_by_protein(dataset)
    uniprot_by_protein = Dict(String(row.protein_id) =>
                                  _parse_residue_ranges(row.residues_binding_in_uniprot)
    for row in eachrow(uniprot_binding_sites))
    usable_statuses = keep_pdb_without_partner ? ("ok", "no_heteromeric_partner") : ("ok",)

    pdb_log = DataFrame(protein_disprot_id = String[], protein_id = String[],
        PDB = String[], used = Bool[], reason = String[])
    table = _empty_aggregated_labels_table()

    for protein_rows in groupby(DataFrame(pdb_contacts), :protein_id; sort = true)
        protein_id = String(first(protein_rows.protein_id))
        protein_disprot_id = String(first(protein_rows.protein_disprot_id))
        labels = get(labels_by_protein, protein_id, nothing)
        labels === nothing && continue

        kept = DataFrameRow[]
        for row in eachrow(protein_rows)
            reason = if !(row.status in usable_statuses)
                "status:$(row.status)"
            elseif !_passes_pdb_quality(row.quality_category, row.resolution;
                max_resolution = max_resolution, keep_no_resolution = keep_no_resolution,
                keep_no_info = keep_no_info)
                "low_quality"
            elseif ismissing(row.resolved_positions)
                "no_resolved_position"
            else
                "used"
            end
            push!(pdb_log, (protein_disprot_id, protein_id, String(row.PDB),
                reason == "used", reason))
            reason == "used" && push!(kept, row)
        end

        aggregated = _protein_aggregated_labels_row(protein_disprot_id, protein_id, labels,
            get(uniprot_by_protein, protein_id, Int[]), kept; vote = vote)
        aggregated === nothing || push!(table, aggregated)
    end

    sort!(table, [:protein_disprot_id, :protein_id])
    sort!(pdb_log, [:protein_disprot_id, :protein_id, :PDB])
    (table = table, pdb_log = pdb_log)
end

"""
    _AGGREGATED_LABELS_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_aggregated_labels_table`.
Vector columns are written as `[v1,v2,...]` with one value per residue of the
UniProt sequence (`NaN` = not covered).
"""
const _AGGREGATED_LABELS_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "n_pdb" => "number of PDB entries used (status and quality filters passed, at least one observed residue)",
    "pdb_ids" => "PDB entries used (`;`-separated)",
    "has_uniprot" => "`true` when the protein has UniProt binding sites on disordered residues",
    "disprot" => "DisProt residue labels (`0`/`1`/`2`)",
    "uniprot" => "UniProt track: `2` in a binding site, `0` elsewhere, `NaN` if no UniProt annotation",
    "pdb_N" => "number of used PDB entries where the residue is observed",
    "pdb_proba" => "fraction of those PDB entries where the residue contacts a heteromeric partner (`p_pdb`); contacts only counted on DisProt disordered residues",
    "rate" => "binding evidence across methods: (DisProt `2` + UniProt `2` + PDB vote) / number of covering methods",
    "final_label" => "`1 + rate` when `rate > 0`, otherwise the DisProt label clamped to `[0, 1]`",
    "order_disorder" => "`1` for DisProt disordered residues (labels `1`/`2`), `0` otherwise",
    "disorder_binding" => "`rate` on disordered residues, `NaN` otherwise",
    "Label_PDB" => "`0` ordered, `1 + p_pdb` disordered; residues observed in at least one used PDB",
    "Label_PDB_disorder" => "`0`/`1` DisProt disorder on residues observed in at least one used PDB",
    "Label_PDB_binding" => "`p_pdb` on disordered residues observed in at least one used PDB",
    "Label_Uniprot" => "`0` ordered, `1` disordered, `2` disordered in a UniProt binding site; `NaN` if no UniProt annotation",
    "Label_Uniprot_disorder" => "`0`/`1` DisProt disorder; `NaN` if no UniProt annotation",
    "Label_Uniprot_binding" => "`0`/`1` UniProt binding site on disordered residues; `NaN` if no UniProt annotation"
]
