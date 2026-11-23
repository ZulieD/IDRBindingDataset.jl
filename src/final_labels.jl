using DataFrames

# ---------------------------------------------------------------------------
# Final residue-level labels combining the four evidence sources.
#
# Disorder, per residue: d = mean of the disorder votes of
#   - DisProt: 1 if labelled disordered or binding (1/2), else 0,
#   - chemical shift: `Label_chemical_shift_disorder` of script 22 (fraction of
#     free lists calling the residue disordered), when evaluable.
#   e.g. disordered in DisProt but ordered by NMR: d = 0.5.
# Binding, on residues with d > 0: b = x / n over the sources covering the
# residue (each source weighs 1/n, i.e. 1/4 when all four cover it):
#   - DisProt: always covers the residue; 1 if labelled binding (2), else 0,
#   - UniProt: covers the protein when it has binding sites (script 18);
#     1 inside a binding site, else 0,
#   - PDB: `Label_PDB_binding` of script 20 (fraction of used PDB entries where
#     the residue contacts a heteromeric partner; DisProt disordered residues),
#   - chemical shift: `Label_chemical_shift_binding` of script 22 (fraction of
#     free/bound pairs calling the residue binding; NMR disordered residues).
# Final label = d + b in [0, 2]: 0 ordered, 1 disordered, 2 disordered binding.
# Residues ordered for every disorder source (d = 0) get 0 and no binding value.
# ---------------------------------------------------------------------------

"""
    FINAL_LABEL_SOURCES

Names of the evidence sources combined by `build_final_labels_table`, in the
order used by `_final_residue_labels`.
"""
const FINAL_LABEL_SOURCES = ("disprot", "uniprot", "pdb", "chemical_shift")

"""
    _final_residue_labels(residue_labels, source_bindings; nmr_disorder = nothing)

Combine DisProt `residue_labels` (`0`/`1`/`2`) with the binding vectors of the
other sources and the NMR disorder vector (`NaN` = residue not covered by that
source; `nothing` = source absent for this protein).

Returns `(label, disorder, binding, n_sources)`:
- `disorder`: `d`, mean of the DisProt (`0`/`1`) and NMR disorder votes, in `[0, 1]`,
- `binding`: `b = x / n`, mean binding vote of the covering sources on residues
  with `d > 0`, `NaN` where `d = 0`,
- `label`: `d + b` in `[0, 2]` (`0` where `d = 0`),
- `n_sources`: number of binding sources covering each residue with `d > 0` (`0` elsewhere).
"""
function _final_residue_labels(
        residue_labels::AbstractVector{<:Integer},
        source_bindings::AbstractVector;
        nmr_disorder::Union{Nothing, AbstractVector{<:Real}} = nothing
)
    n_positions = length(residue_labels)
    for vector in vcat(source_bindings, [nmr_disorder])
        vector === nothing && continue
        length(vector) == n_positions || throw(DimensionMismatch(
            "Source vector has $(length(vector)) values for a sequence of $(n_positions) residues."))
    end

    label, disorder = zeros(n_positions), zeros(n_positions)
    binding_rate = fill(NaN, n_positions)
    n_sources = zeros(Int, n_positions)
    for position in 1:n_positions
        disprot = residue_labels[position]
        disorder_votes = Float64(disprot in
                                 (_RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING))
        n_disorder = 1
        if nmr_disorder !== nothing && !isnan(nmr_disorder[position])
            disorder_votes += clamp(Float64(nmr_disorder[position]), 0.0, 1.0)
            n_disorder += 1
        end
        d = disorder_votes / n_disorder
        disorder[position] = d
        d > 0 || continue

        n = 1
        x = Float64(disprot == _RESIDUE_LABEL_BINDING)
        for binding in source_bindings
            binding === nothing && continue
            value = binding[position]
            isnan(value) && continue
            n += 1
            x += clamp(Float64(value), 0.0, 1.0)
        end
        b = x / n
        binding_rate[position] = b
        label[position] = d + b
        n_sources[position] = n
    end
    (label = label, disorder = disorder, binding = binding_rate, n_sources = n_sources)
end

"""
    _uniprot_binding_vector(positions, n_positions)

UniProt binding vector: `nothing` when the protein has no UniProt binding site,
otherwise `1` inside a binding site and `0` elsewhere.
"""
function _uniprot_binding_vector(positions::AbstractVector{<:Integer}, n_positions::Integer)
    isempty(positions) && return nothing
    vector = zeros(n_positions)
    for position in positions
        1 <= position <= n_positions && (vector[position] = 1.0)
    end
    vector
end

"""
    _binding_vectors_by_protein(table, column, name)

Map `protein_id` to the parsed binding vector of `column` (written by
`_format_label_vector`).
"""
function _binding_vectors_by_protein(
        table::AbstractDataFrame, column::Symbol, name::AbstractString
)::Dict{String, Vector{Float64}}
    _require_columns(table, [:protein_id, column], name)
    Dict(String(row.protein_id) => _parse_label_vector(string(row[column]))
    for row in eachrow(table) if !ismissing(row[column]))
end

"""
    build_final_labels_table(dataset, uniprot_binding_sites, aggregated_labels,
        chemical_shift_labels)

Combine DisProt labels, UniProt binding sites (script 18), PDB contacts
(`Label_PDB_binding` of script 20) and chemical shifts
(`Label_chemical_shift_disorder` and `Label_chemical_shift_binding` of script
22) into final residue labels, one row per dataset protein.

Disorder `d` is the mean of the DisProt (`0`/`1`) and NMR disorder votes.
On residues with `d > 0`, binding `b = x / n` where `n` is the number of
sources covering the residue and `x` the sum of their binding votes (DisProt
and UniProt vote `0`/`1`, PDB and chemical shift give a fraction); each source
covering the residue weighs the same (`1/4` when all four do).

Returns `(table, sources)`:
- `table`: `protein_disprot_id`, `protein_id`, `Label_aggregate_all` (`d + b`
  in `[0, 2]`), `Label_aggregate_all_disorder` (`d` in `[0, 1]`),
  `Label_aggregate_all_binding` (`b` in `[0, 1]`, `NaN` where `d = 0`),
- `sources`: one row per protein with the sources available (`;`-separated),
  the number of residues with `d > 0`, with `0 < d < 1`, disordered for NMR
  only, and covered by 1, 2, 3 or 4 binding sources.
"""
function build_final_labels_table(
        dataset::AbstractDataFrame,
        uniprot_binding_sites::AbstractDataFrame,
        aggregated_labels::AbstractDataFrame,
        chemical_shift_labels::AbstractDataFrame
)
    _require_columns(dataset, [:protein_id, :protein_disprot_id, :residue_labels], "dataset")
    _require_columns(uniprot_binding_sites, [:protein_id, :residues_binding_in_uniprot],
        "uniprot_binding_sites")

    uniprot_by_protein = Dict(String(row.protein_id) =>
                                  _parse_residue_ranges(row.residues_binding_in_uniprot)
    for row in eachrow(uniprot_binding_sites))
    pdb_by_protein = _binding_vectors_by_protein(aggregated_labels, :Label_PDB_binding,
        "aggregated_labels")
    nmr_by_protein = _binding_vectors_by_protein(chemical_shift_labels,
        :Label_chemical_shift_binding, "chemical_shift_labels")
    nmr_disorder_by_protein = _binding_vectors_by_protein(chemical_shift_labels,
        :Label_chemical_shift_disorder, "chemical_shift_labels")

    table = DataFrame(protein_disprot_id = String[], protein_id = String[],
        Label_aggregate_all = String[], Label_aggregate_all_disorder = String[],
        Label_aggregate_all_binding = String[])
    sources = DataFrame(protein_disprot_id = String[], protein_id = String[],
        sources = String[], n_disordered = Int[], n_partial_disorder = Int[],
        n_nmr_only_disorder = Int[], n_1_source = Int[],
        n_2_sources = Int[], n_3_sources = Int[], n_4_sources = Int[])

    for row in eachrow(dataset)
        ismissing(row.residue_labels) && continue
        protein_id = String(row.protein_id)
        residue_labels = _decode_residue_labels(string(row.residue_labels))
        n_positions = length(residue_labels)

        bindings = Union{Nothing, Vector{Float64}}[
            _uniprot_binding_vector(get(uniprot_by_protein, protein_id, Int[]), n_positions),
            get(pdb_by_protein, protein_id, nothing),
            get(nmr_by_protein, protein_id, nothing)]
        nmr_disorder = get(nmr_disorder_by_protein, protein_id, nothing)
        labels = _final_residue_labels(residue_labels, bindings;
            nmr_disorder = nmr_disorder)

        covered = [binding !== nothing && any(!isnan, binding) for binding in bindings]
        covered[end] |= nmr_disorder !== nothing && any(!isnan, nmr_disorder)
        available = vcat(["disprot"],
            [name
             for (name, is_covered) in zip(FINAL_LABEL_SOURCES[2:end], covered)
             if is_covered])
        in_disprot_idr = [label in (_RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING)
                          for label in residue_labels]

        push!(table,
            (String(row.protein_disprot_id), protein_id,
                _format_label_vector(labels.label), _format_label_vector(labels.disorder),
                _format_label_vector(labels.binding)))
        push!(sources,
            (String(row.protein_disprot_id), protein_id, join(available, ";"),
                count(>(0), labels.disorder), count(d -> 0 < d < 1, labels.disorder),
                count(p -> labels.disorder[p] > 0 && !in_disprot_idr[p], 1:n_positions),
                count(==(1), labels.n_sources),
                count(==(2), labels.n_sources), count(==(3), labels.n_sources),
                count(==(4), labels.n_sources)))
    end

    sort!(table, [:protein_disprot_id, :protein_id])
    sort!(sources, [:protein_disprot_id, :protein_id])
    (table = table, sources = sources)
end

"""
    _FINAL_LABELS_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_final_labels_table`. Vector
columns are written as `[v1,v2,...]` with one value per residue of the UniProt
sequence.
"""
const _FINAL_LABELS_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "Label_aggregate_all" => "`Label_aggregate_all_disorder + Label_aggregate_all_binding` in `[0, 2]`: `0` ordered, `1` disordered, `2` disordered binding (all sources agree); `0` where no source calls the residue disordered",
    "Label_aggregate_all_disorder" => "`d` in `[0, 1]`: mean disorder vote of DisProt (`0`/`1`) and NMR (`Label_chemical_shift_disorder`, when evaluable); e.g. `0.5` when DisProt calls the residue disordered and NMR ordered",
    "Label_aggregate_all_binding" => "`b` in `[0, 1]` on residues with `d > 0`: mean binding vote of the covering sources (DisProt `0`/`1`, UniProt `0`/`1`, `Label_PDB_binding`, `Label_chemical_shift_binding`); `NaN` where `d = 0`"
]

"""
    _FINAL_LABEL_SOURCES_COLUMN_DEFINITIONS

Data dictionary of the per-protein source table returned by `build_final_labels_table`.
"""
const _FINAL_LABEL_SOURCES_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein",
    "protein_id" => "UniProt accession of the dataset protein",
    "sources" => "sources with at least one covered residue (`disprot`, `uniprot`, `pdb`, `chemical_shift`; `;`-separated)",
    "n_disordered" => "residues called disordered by at least one source (`d > 0`)",
    "n_partial_disorder" => "residues where DisProt and NMR disagree on disorder (`0 < d < 1`)",
    "n_nmr_only_disorder" => "residues ordered in DisProt but disordered for NMR",
    "n_1_source" => "disordered residues with DisProt as only binding source",
    "n_2_sources" => "disordered residues covered by 2 binding sources",
    "n_3_sources" => "disordered residues covered by 3 binding sources",
    "n_4_sources" => "disordered residues covered by all 4 binding sources"
]
