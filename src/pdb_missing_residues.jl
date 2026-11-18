using DataFrames

# ---------------------------------------------------------------------------
# PDB segments of dataset proteins and their missing (unobserved) residues.
#
# One "PDB segment" is one (protein_id, PDB, SP_BEG, SP_END) combination from
# the SIFTS chain-level mapping `pdb_chain_uniprot.csv.gz`; all PDB chains that
# cover the same UniProt range are merged in `CHAIN`. SIFTS columns
# RES_BEG/RES_END/PDB_BEG/PDB_END are not used (and often missing), so they
# are dropped.
# ---------------------------------------------------------------------------
const PDB_SEGMENT_KEY_COLUMNS = [:protein_id, :PDB, :SP_BEG, :SP_END]

"""
    _parse_sifts_position(value)

Parse one SIFTS position value into `Int`, returning `missing` for empty or
non-integer values.
"""
function _parse_sifts_position(value)::Union{Missing, Int}
    ismissing(value) && return missing
    value isa Integer && return Int(value)
    parsed = tryparse(Int, strip(string(value)))
    parsed === nothing ? missing : parsed
end

"""
    build_dataset_pdb_segments(dataset, sifts_chain_uniprot)::DataFrame

Join dataset proteins (`protein_id`, UniProt accession) with the SIFTS chain-level mapping
(`SP_PRIMARY`) and merge chains covering the same UniProt range.

Returns one row per `(protein_id, PDB, SP_BEG, SP_END)` with columns
`protein_disprot_id`, `protein_id`, `PDB` (lowercase entry code), `CHAIN` (`;`-separated, sorted),
`SP_BEG`, `SP_END`. Proteins without PDB mapping are not included. SIFTS rows
with non-integer `SP_BEG`/`SP_END` are dropped.
"""
function build_dataset_pdb_segments(
        dataset::AbstractDataFrame,
        sifts_chain_uniprot::AbstractDataFrame
)::DataFrame
    _require_columns(dataset, [:protein_id, :protein_disprot_id], "dataset")
    _require_columns(sifts_chain_uniprot, [:PDB, :CHAIN, :SP_PRIMARY, :SP_BEG, :SP_END],
        "sifts_chain_uniprot")

    sifts = DataFrame(
        protein_id = String.(strip.(string.(sifts_chain_uniprot.SP_PRIMARY))),
        PDB = _pdb_entry_id.(sifts_chain_uniprot.PDB),
        CHAIN = String.(strip.(string.(sifts_chain_uniprot.CHAIN))),
        SP_BEG = _parse_sifts_position.(sifts_chain_uniprot.SP_BEG),
        SP_END = _parse_sifts_position.(sifts_chain_uniprot.SP_END)
    )
    dropmissing!(sifts, [:SP_BEG, :SP_END])
    sifts.SP_BEG = Int.(sifts.SP_BEG)
    sifts.SP_END = Int.(sifts.SP_END)

    proteins = unique(DataFrame(protein_id = String.(dataset.protein_id),
        protein_disprot_id = String.(dataset.protein_disprot_id)))
    joined = innerjoin(proteins, sifts; on = :protein_id)

    segments = combine(
        groupby(joined, vcat(PDB_SEGMENT_KEY_COLUMNS, [:protein_disprot_id])),
        :CHAIN => (chains -> _format_chain_list(sort!(unique(chains)))) => :CHAIN)
    sort!(segments, PDB_SEGMENT_KEY_COLUMNS)
    select!(segments, :protein_disprot_id, :protein_id, :PDB, :CHAIN, :SP_BEG, :SP_END)
end

"""
    _count_missing_residues(records, labels, accession, sp_beg, sp_end)

Count residues of one PDB chain without coordinates, restricted to residues
mapped to UniProt `accession` inside `sp_beg:sp_end`. DisProt `labels` (full
sequence, `0`/`1`/`2`) split the count by class.

Returns `(n_missing, n_missing_disorder, n_missing_disorder_binding)`, where
`n_missing_disorder` counts label `1` and `n_missing_disorder_binding` label `2`.
"""
function _count_missing_residues(
        records::AbstractVector{SIFTSResidueRecord},
        labels::AbstractVector{<:Integer},
        accession::AbstractString,
        sp_beg::Integer,
        sp_end::Integer
)
    target = uppercase(strip(accession))
    counted = Set{Int}()
    n_missing = n_disorder = n_binding = 0
    for record in records
        record.is_missing || continue
        record.uniprot_id == target || continue
        position = record.uniprot_position
        (sp_beg <= position <= sp_end && 1 <= position <= length(labels)) || continue
        position in counted && continue
        push!(counted, position)
        n_missing += 1
        labels[position] == _RESIDUE_LABEL_DISORDERED && (n_disorder += 1)
        labels[position] == _RESIDUE_LABEL_BINDING && (n_binding += 1)
    end
    (n_missing = n_missing, n_missing_disorder = n_disorder,
        n_missing_disorder_binding = n_binding)
end

"""
    _select_best_chains(chain_counts)

Keep the chains with the fewest missing residues. `chain_counts` is a vector of
`chain => counts` pairs (counts as returned by `_count_missing_residues`).
Returns `(chains, counts)`; ties are all kept, in input order. Returns
`(String[], nothing)` for an empty input.
"""
function _select_best_chains(chain_counts::AbstractVector)
    isempty(chain_counts) && return (String[], nothing)
    best = minimum(last(pair).n_missing for pair in chain_counts)
    kept = [pair for pair in chain_counts if last(pair).n_missing == best]
    (String[first(pair) for pair in kept], last(first(kept)))
end

"""
    _segment_missing_residue_row(segment, labels, chain_records)

Compute the missing-residue columns of one PDB segment row. `chain_records` is
the per-chain output of `_read_sifts_residue_records`, or `nothing` when the
SIFTS XML file is unavailable (counts are then `missing` and every chain is
kept).
"""
function _segment_missing_residue_row(
        segment,
        labels::AbstractVector{<:Integer},
        chain_records::Union{Nothing, AbstractDict}
)
    chains = _parse_chain_list(segment.CHAIN)
    if chain_records === nothing
        return (CHAIN_KEPT = _format_chain_list(chains),
            N_MISSING = missing, N_MISSING_DISORDER = missing,
            N_MISSING_DISORDER_BINDING = missing)
    end

    chain_counts = Pair{String, Any}[]
    for chain in chains
        haskey(chain_records, chain) || continue
        push!(chain_counts,
            chain => _count_missing_residues(chain_records[chain], labels,
                segment.protein_id, segment.SP_BEG, segment.SP_END))
    end
    kept, counts = _select_best_chains(chain_counts)
    counts === nothing && return (CHAIN_KEPT = missing,
        N_MISSING = missing, N_MISSING_DISORDER = missing,
        N_MISSING_DISORDER_BINDING = missing)
    (CHAIN_KEPT = _format_chain_list(kept), N_MISSING = counts.n_missing,
        N_MISSING_DISORDER = counts.n_missing_disorder,
        N_MISSING_DISORDER_BINDING = counts.n_missing_disorder_binding)
end

"""
    _empty_missing_residue_table()

Return the empty result table (schema) of the missing-residue computation.
"""
function _empty_missing_residue_table()::DataFrame
    DataFrame(protein_id = String[], PDB = String[], SP_BEG = Int[], SP_END = Int[],
        CHAIN_KEPT = Union{Missing, String}[], N_MISSING = Union{Missing, Int}[],
        N_MISSING_DISORDER = Union{Missing, Int}[],
        N_MISSING_DISORDER_BINDING = Union{Missing, Int}[])
end

"""
    _pdb_entry_missing_residue_table(pdb, segments, labels_by_protein;
        sifts_xml_dirs, sifts_download_dir = nothing, logger = println)

Compute missing-residue columns for all segments of one PDB entry, parsing its
SIFTS XML file only once.
"""
function _pdb_entry_missing_residue_table(
        pdb::AbstractString,
        segments::AbstractDataFrame,
        labels_by_protein::AbstractDict;
        sifts_xml_dirs,
        sifts_download_dir::Union{Nothing, AbstractString} = nothing,
        logger::Function = println
)::DataFrame
    path = _sifts_xml_path(pdb; sifts_xml_dirs = sifts_xml_dirs,
        download_dir = sifts_download_dir, logger = logger)
    path === nothing && logger("SIFTS XML not found for $(pdb); N_MISSING* left missing.")
    chain_records = path === nothing ? nothing : _read_sifts_residue_records(path)

    table = _empty_missing_residue_table()
    for segment in eachrow(segments)
        labels = get(labels_by_protein, segment.protein_id, nothing)
        labels === nothing && continue
        valid = 1 <= segment.SP_BEG <= segment.SP_END <= length(labels)
        valid ||
            logger("Invalid SP_BEG/SP_END for $(segment.protein_id)/$(pdb); row skipped.")
        valid || continue
        push!(table,
            merge(
                (protein_id = segment.protein_id, PDB = segment.PDB,
                    SP_BEG = segment.SP_BEG, SP_END = segment.SP_END),
                _segment_missing_residue_row(segment, labels, chain_records)))
    end
    table
end

"""
    _dataset_residue_labels_by_protein(dataset)

Return `Dict(protein_id => Vector{Int})` from the compact `residue_labels`
strings of the dataset.
"""
function _dataset_residue_labels_by_protein(dataset::AbstractDataFrame)::Dict{
        String, Vector{Int}}
    _require_columns(dataset, [:protein_id, :residue_labels], "dataset")
    Dict(String(row.protein_id) => _decode_residue_labels(string(row.residue_labels))
    for row in eachrow(dataset) if !ismissing(row.residue_labels))
end

"""
    build_pdb_missing_residues_table(dataset, segments;
        sifts_xml_dirs,
        sifts_download_dir = nothing,
        checkpoint_dir = nothing,
        batch_size = DEFAULT_CHECKPOINT_BATCH_SIZE,
        logger = println
    )

Add missing-residue statistics to the PDB `segments` returned by
`build_dataset_pdb_segments`, using one SIFTS XML file per PDB entry.

Added columns:
- `COMPLETE`: the segment covers the whole UniProt sequence
- `CHAIN_KEPT`: chains of `CHAIN` with the fewest missing residues (`;`-separated)
- `N_MISSING`: residues of the segment without coordinates in the kept chain(s)
- `N_MISSING_DISORDER`: same, restricted to DisProt label `1`
- `N_MISSING_DISORDER_BINDING`: same, restricted to DisProt label `2`

Only residues mapped by SIFTS to the segment UniProt accession are counted
(expression tags and other unmapped residues are ignored). PDB entries are
processed in parallel and can be checkpointed (see `run_checkpointed_batches`).
Returns `(table, failed_pdbs)`.
"""
function build_pdb_missing_residues_table(
        dataset::AbstractDataFrame,
        segments::AbstractDataFrame;
        sifts_xml_dirs,
        sifts_download_dir::Union{Nothing, AbstractString} = nothing,
        checkpoint_dir::Union{Nothing, AbstractString} = nothing,
        batch_size::Integer = DEFAULT_CHECKPOINT_BATCH_SIZE,
        logger::Function = println
)
    _require_columns(dataset, [:protein_id, :sequence, :residue_labels], "dataset")
    _require_columns(segments,
        [:protein_disprot_id, :protein_id, :PDB, :CHAIN, :SP_BEG, :SP_END], "segments")

    labels_by_protein = _dataset_residue_labels_by_protein(dataset)
    segments_by_pdb = groupby(segments, :PDB)
    pdbs = sort!(unique(String.(segments.PDB)))

    result = run_checkpointed_batches(pdbs,
        (missing_residues = _empty_missing_residue_table(),);
        checkpoint_dir = checkpoint_dir, batch_size = batch_size,
        logger = logger) do pdb
        (missing_residues = _pdb_entry_missing_residue_table(pdb,
            segments_by_pdb[(PDB = pdb,)], labels_by_protein;
            sifts_xml_dirs = sifts_xml_dirs, sifts_download_dir = sifts_download_dir,
            logger = logger),)
    end

    sequence_lengths = Dict(String(row.protein_id) => length(string(row.sequence))
    for row in eachrow(dataset) if !ismissing(row.sequence))
    table = leftjoin(DataFrame(segments),
        unique(result.tables.missing_residues, PDB_SEGMENT_KEY_COLUMNS);
        on = PDB_SEGMENT_KEY_COLUMNS)
    table.COMPLETE = [row.SP_BEG == 1 &&
                          row.SP_END == get(sequence_lengths, row.protein_id, -1)
                      for row in eachrow(table)]
    select!(table, :protein_disprot_id, :protein_id, :PDB, :CHAIN, :SP_BEG, :SP_END,
        :COMPLETE, :CHAIN_KEPT,
        :N_MISSING, :N_MISSING_DISORDER, :N_MISSING_DISORDER_BINDING)
    sort!(table, PDB_SEGMENT_KEY_COLUMNS)
    (table = table, failed_pdbs = result.failed_items)
end

"""
    _PDB_MISSING_RESIDUES_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_pdb_missing_residues_table`.
"""
const _PDB_MISSING_RESIDUES_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "PDB" => "PDB entry code (lowercase), from SIFTS `pdb_chain_uniprot.csv.gz`",
    "CHAIN" => "PDB author chains covering the same UniProt range (`;`-separated)",
    "SP_BEG" => "first UniProt position covered by the PDB chain(s) (1-based)",
    "SP_END" => "last UniProt position covered by the PDB chain(s) (1-based, inclusive)",
    "COMPLETE" => "`true` when `SP_BEG:SP_END` covers the whole UniProt sequence",
    "CHAIN_KEPT" => "chains of `CHAIN` with the fewest missing residues (`;`-separated); all chains when the SIFTS XML file is unavailable",
    "N_MISSING" => "residues of `SP_BEG:SP_END` without coordinates in the kept chain(s) (SIFTS `Not_Observed`); `missing` when the SIFTS XML file is unavailable",
    "N_MISSING_DISORDER" => "same as `N_MISSING`, restricted to DisProt label `1`",
    "N_MISSING_DISORDER_BINDING" => "same as `N_MISSING`, restricted to DisProt label `2`"
]
