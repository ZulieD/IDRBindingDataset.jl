using CSV
using DataFrames

# ---------------------------------------------------------------------------
# Heteromeric interface contacts of dataset proteins in their PDB structures.
#
# For one protein and one PDB entry: the protein chain(s) are the SIFTS chains
# of `disprot_pdb_mapping.csv.gz`; partner chains are the other polymer chains
# that are not copies of the protein (homo-oligomers are excluded). A protein
# residue is in contact when one of its heavy atoms is within the cutoff of a
# heavy atom of a partner chain. Residues are mapped to UniProt positions with
# SIFTS XML. All contacts are reported; filters (disorder, resolution, ...) are
# applied when labels are aggregated.
#
# Thresholds are project choices (not taken from a reference), kept as named
# constants and exposed as keyword arguments.
# ---------------------------------------------------------------------------
const PDB_CONTACT_CUTOFF_ANGSTROM = 6.0     # heavy-atom distance for a contact
const PDB_CONTACT_HOMO_KMER_LENGTH = 5      # k-mer length for homo-chain detection
const PDB_CONTACT_HOMO_KMER_MIN_FRACTION = 0.8
const PDB_CONTACT_MIN_MAPPING_IDENTITY = 0.9  # below: warning (doubtful mapping)
const RCSB_MMCIF_DOWNLOAD_URL = "https://files.rcsb.org/download/%s.cif"
const MMCIF_EXTENSIONS = (".cif", ".cif.gz")

"""
    _AMINO_ACID_THREE_TO_ONE

Three-letter to one-letter amino-acid codes, including common modified residues
(`MSE` -> `M`, `SEC` -> `U`, `PYL` -> `O`).
"""
const _AMINO_ACID_THREE_TO_ONE = Dict(
    "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D', "CYS" => 'C',
    "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G', "HIS" => 'H', "ILE" => 'I',
    "LEU" => 'L', "LYS" => 'K', "MET" => 'M', "PHE" => 'F', "PRO" => 'P',
    "SER" => 'S', "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
    "MSE" => 'M', "SEC" => 'U', "PYL" => 'O')

"""
    _one_letter_code(name)::Char

Return the one-letter amino-acid code of residue `name`, or `'X'`.
"""
_one_letter_code(name::AbstractString)::Char = get(
    _AMINO_ACID_THREE_TO_ONE, uppercase(name), 'X')

# ---------------------------------------------------------------------------
# PDB quality (method / resolution)
# ---------------------------------------------------------------------------
"""
    _parse_resolution(value)

Parse a resolution value (`2.1`, `"2.10"`, `"1.9, 2.3"`) into the smallest
`Float64`; `missing` when no number is present (e.g. `"NOT APPLICABLE"`).
"""
function _parse_resolution(value)::Union{Missing, Float64}
    ismissing(value) && return missing
    value isa Real && return isnan(value) ? missing : Float64(value)
    numbers = [parse(Float64, m.match) for m in eachmatch(r"\d+(?:\.\d+)?", string(value))]
    isempty(numbers) ? missing : minimum(numbers)
end

"""
    read_pdb_quality_table(path; id_column = nothing)

Read a PDB information CSV with `RESOLUTION` and `METHOD` columns and return
`Dict(pdb_entry => (method, resolution))`. The PDB identifier column is
`id_column`, or is detected among common names (`PDB`, `pdb_id`, `IDCODE`, ...).
"""
function read_pdb_quality_table(
        path::AbstractString;
        id_column::Union{Nothing, AbstractString} = nothing
)::Dict{
        String, NamedTuple{(:method, :resolution), Tuple{String, Union{Missing, Float64}}}}
    table = CSV.read(path, DataFrame; stringtype = String)
    _require_columns(table, [:RESOLUTION, :METHOD], "PDB quality table")
    candidates = ("PDB", "pdb", "PDB_ID", "pdb_id", "PDBID", "pdbid", "IDCODE", "idcode",
        "entry_id", "ENTRY_ID", "structureId", "ID", "id")
    column = if id_column === nothing
        index = findfirst(in(names(table)), candidates)
        index === nothing ? "" : candidates[index]
    else
        String(id_column)
    end
    column in names(table) ||
        throw(ArgumentError("PDB identifier column not found in $(path); set `id_column`."))

    quality = Dict{
        String, NamedTuple{(:method, :resolution), Tuple{String, Union{Missing, Float64}}}}()
    for row in eachrow(table)
        id = _normalize_optional_text(row[column])
        ismissing(id) && continue
        entry = _pdb_entry_id(id)
        method = something(
            _normalize_optional_text(row.METHOD; uppercase_output = true,
                missing_value = nothing),
            "missing")
        resolution = _parse_resolution(row.RESOLUTION)
        if haskey(quality, entry)  # keep the first known values
            previous = quality[entry]
            resolution = ismissing(previous.resolution) ? resolution : previous.resolution
            method = previous.method == "missing" ? method : previous.method
        end
        quality[entry] = (method = method, resolution = resolution)
    end
    quality
end

"""
    _pdb_quality_category(quality)

Return `:resolved` (resolution known), `:no_resolution` (e.g. NMR) or `:absent`
(`quality === nothing`, entry not in the quality table).
"""
function _pdb_quality_category(quality)::Symbol
    quality === nothing && return :absent
    ismissing(quality.resolution) ? :no_resolution : :resolved
end

# ---------------------------------------------------------------------------
# Structures (mmCIF) and chains
# ---------------------------------------------------------------------------
"""
    _mmcif_path(pdb; mmcif_dirs, download_dir = nothing, logger = println)

Return the mmCIF path of `pdb` found in `mmcif_dirs`, downloading it from RCSB
into `download_dir` when needed. Returns `nothing` when unavailable.
"""
function _mmcif_path(
        pdb::AbstractString;
        mmcif_dirs,
        download_dir::Union{Nothing, AbstractString} = nothing,
        logger::Function = println
)::Union{Nothing, String}
    search_dirs = download_dir === nothing ? collect(mmcif_dirs) :
                  vcat(collect(mmcif_dirs), [download_dir])
    path = _find_pdb_entry_file(search_dirs, pdb, MMCIF_EXTENSIONS)
    (path !== nothing || download_dir === nothing) && return path

    entry = _pdb_entry_id(pdb)
    destination = joinpath(download_dir, entry * ".cif")
    try
        download_to_file_with_retries(
            replace(RCSB_MMCIF_DOWNLOAD_URL, "%s" => uppercase(entry)),
            destination; logger = logger, log_success = false)
        return destination
    catch err
        logger("mmCIF download failed for $(entry): $(sprint(showerror, err))")
        return nothing
    end
end

"""
    _is_hydrogen(atom)::Bool

Return `true` for hydrogen and deuterium atoms.
"""
_is_hydrogen(atom)::Bool = uppercase(strip(atom.element)) in ("H", "D")

"""
    _load_first_model_heavy_atoms(path)

Read one mmCIF file and return the residues of its first model with hydrogens
removed (`Vector{MIToS.PDB.PDBResidue}`).
"""
function _load_first_model_heavy_atoms(path::AbstractString)::Vector{MIToS.PDB.PDBResidue}
    residues = MIToS.PDB.read_file(path, MIToS.PDB.MMCIFFile)
    isempty(residues) && return MIToS.PDB.PDBResidue[]
    model = first(residues).id.model
    output = MIToS.PDB.PDBResidue[]
    for residue in residues
        residue.id.model == model || continue
        atoms = filter(!_is_hydrogen, residue.atoms)
        isempty(atoms) || push!(output, MIToS.PDB.PDBResidue(residue.id, atoms))
    end
    output
end

"""
    _is_polymer_residue(residue)::Bool

Return `true` for polymer residues (`ATOM` records, plus modified amino acids).
"""
_is_polymer_residue(residue)::Bool = residue.id.group == "ATOM" ||
                                     haskey(_AMINO_ACID_THREE_TO_ONE, uppercase(residue.id.name))

"""
    _polymer_chains(residues)::Vector{String}

Return the chains having at least one `ATOM` residue, in file order.
"""
_polymer_chains(residues)::Vector{String} = unique(String[residue.id.chain
                                                          for residue in residues
                                                          if residue.id.group == "ATOM"])

"""
    _chain_sequence(residues, chain)::String

Return the one-letter sequence of the polymer residues of `chain`.
"""
function _chain_sequence(residues, chain::AbstractString)::String
    String([_one_letter_code(residue.id.name)
            for residue in residues
            if residue.id.chain == chain && _is_polymer_residue(residue)])
end

"""
    _kmer_fraction(query, reference; k = PDB_CONTACT_HOMO_KMER_LENGTH)::Float64

Fraction of the `k`-mers of `query` present in `reference`.
"""
function _kmer_fraction(
        query::AbstractString,
        reference::AbstractString;
        k::Integer = PDB_CONTACT_HOMO_KMER_LENGTH
)::Float64
    isempty(query) && return 0.0
    length(query) < k && return occursin(query, reference) ? 1.0 : 0.0
    reference_chars = collect(reference)
    query_chars = collect(query)
    reference_kmers = Set(String(reference_chars[i:(i + k - 1)])
    for i in 1:(length(reference_chars) - k + 1))
    n_kmers = length(query_chars) - k + 1
    count(i -> String(query_chars[i:(i + k - 1)]) in reference_kmers, 1:n_kmers) / n_kmers
end

"""
    _homologous_chains(chain_sequences, target_chains, protein_sequence;
        min_fraction = PDB_CONTACT_HOMO_KMER_MIN_FRACTION)::Set{String}

Return the chains that are copies of the target protein: the `target_chains`
(SIFTS) plus any chain whose k-mer fraction against `protein_sequence` is at
least `min_fraction`. `chain_sequences` maps chain => sequence.
"""
function _homologous_chains(
        chain_sequences::AbstractDict,
        target_chains,
        protein_sequence::AbstractString;
        min_fraction::Real = PDB_CONTACT_HOMO_KMER_MIN_FRACTION
)::Set{String}
    homologous = Set{String}(c for c in target_chains if haskey(chain_sequences, c))
    for (chain, sequence) in chain_sequences
        chain in homologous && continue
        _kmer_fraction(sequence, protein_sequence) >= min_fraction &&
            push!(homologous, chain)
    end
    homologous
end

# ---------------------------------------------------------------------------
# Contacts
# ---------------------------------------------------------------------------
"""
    _residue_center_and_radius(residue)

Return the atom centroid of `residue` and the largest atom distance to it.
"""
function _residue_center_and_radius(residue)
    coordinates = [atom.coordinates for atom in residue.atoms]
    center = sum(coordinates) / length(coordinates)
    radius = maximum(sqrt(sum(abs2, c - center)) for c in coordinates)
    (center, radius)
end

"""
    _interface_residue_indices(query, partner, cutoff)::Set{Int}

Indices of `query` residues with at least one atom within `cutoff` Å of an atom
of a `partner` residue. A bounding-sphere test skips distant residue pairs
before the exact `MIToS.PDB.contact` check.
"""
function _interface_residue_indices(query, partner, cutoff::Real)::Set{Int}
    partner_residues = [residue for residue in partner if !isempty(residue.atoms)]
    hits = Set{Int}()
    isempty(partner_residues) && return hits
    partner_spheres = [_residue_center_and_radius(residue) for residue in partner_residues]
    for (index, residue) in enumerate(query)
        isempty(residue.atoms) && continue
        center, radius = _residue_center_and_radius(residue)
        for (partner_residue, (partner_center, partner_radius)) in zip(partner_residues, partner_spheres)
            sqrt(sum(abs2, center - partner_center)) > radius + partner_radius + cutoff &&
                continue
            if MIToS.PDB.contact(residue, partner_residue, Float64(cutoff))
                push!(hits, index)
                break
            end
        end
    end
    hits
end

"""
    _chain_contact_positions(query, pdb_to_uniprot, hits, protein_sequence)

Map one chain's residues to UniProt positions and split interface `hits`.
Returns `(resolved, contacts, n_unmapped_contacts, n_mapped, n_identical)`:
`resolved`/`contacts` are UniProt positions of observed / in-contact residues,
`n_unmapped_contacts` counts contacts without UniProt mapping and `n_identical`
the mapped residues whose amino acid matches `protein_sequence`.
"""
function _chain_contact_positions(
        query,
        pdb_to_uniprot::AbstractDict{String, Int},
        hits::AbstractSet{Int},
        protein_sequence::AbstractString
)
    sequence = collect(protein_sequence)
    resolved, contacts = Set{Int}(), Set{Int}()
    n_unmapped = n_mapped = n_identical = 0
    for (index, residue) in enumerate(query)
        position = get(pdb_to_uniprot, residue.id.number, 0)
        if !(1 <= position <= length(sequence))
            index in hits && (n_unmapped += 1)
            continue
        end
        n_mapped += 1
        _one_letter_code(residue.id.name) == sequence[position] && (n_identical += 1)
        push!(resolved, position)
        index in hits && push!(contacts, position)
    end
    (resolved = resolved, contacts = contacts, n_unmapped_contacts = n_unmapped,
        n_mapped = n_mapped, n_identical = n_identical)
end

# ---------------------------------------------------------------------------
# Contacts table (one row per protein and PDB entry)
# ---------------------------------------------------------------------------
"""
    _format_position_ranges(positions)::Union{Missing, String}

Encode residue positions as compact `;`-separated ranges (`[1, 2, 3, 7]` ->
`"1-3;7"`), readable with `_parse_residue_ranges`. Returns `missing` when
`positions` is empty.
"""
function _format_position_ranges(positions)::Union{Missing, String}
    sorted = sort!(unique(collect(Int, positions)))
    isempty(sorted) && return missing
    ranges = String[]
    range_start = range_stop = first(sorted)
    for position in Iterators.drop(sorted, 1)
        if position == range_stop + 1
            range_stop = position
            continue
        end
        push!(ranges, range_start == range_stop ? string(range_start) :
                      "$(range_start)-$(range_stop)")
        range_start = range_stop = position
    end
    push!(ranges, range_start == range_stop ? string(range_start) :
                  "$(range_start)-$(range_stop)")
    join(ranges, ";")
end

"""
    _optional_chain_list(chains)::Union{Missing, String}

`_format_chain_list(chains)`, or `missing` for an empty list.
"""
_optional_chain_list(chains)::Union{Missing, String} = isempty(chains) ? missing :
                                                       _format_chain_list(chains)

"""
    PDB_CONTACT_STATUSES

Values of the `status` column of the contacts table:
- `ok`: contacts computed against at least one heteromeric partner chain
- `no_heteromeric_partner`: only copies of the protein (or no other polymer);
  observed positions are reported, there are no contacts
- `no_mmcif` / `no_sifts`: structure or SIFTS XML file unavailable
- `no_target_chain`: none of the SIFTS chains is a polymer chain of the mmCIF
- `no_mapped_residue`: no observed residue mapped to the protein UniProt sequence
"""
const PDB_CONTACT_STATUSES = ("ok", "no_heteromeric_partner", "no_mmcif", "no_sifts",
    "no_target_chain", "no_mapped_residue")

"""
    _empty_pdb_contacts_table()

Return the empty contacts table (schema) of `build_pdb_contacts_table`.
"""
function _empty_pdb_contacts_table()::DataFrame
    DataFrame(protein_disprot_id = String[], protein_id = String[], PDB = String[],
        method = String[], resolution = Union{Missing, Float64}[],
        quality_category = String[], status = String[],
        chains = Union{Missing, String}[], partner_chains = Union{Missing, String}[],
        homo_chains_excluded = Union{Missing, String}[], n_mapped = Int[],
        identity = Float64[], n_resolved = Int[],
        resolved_positions = Union{Missing, String}[], n_contacts = Int[],
        contact_positions = Union{Missing, String}[], n_contacts_unmapped = Int[])
end

"""
    _protein_pdb_contacts_row(protein_disprot_id, protein_id, pdb, target_chains,
        sequence; quality, mmcif_dirs, mmcif_download_dir, sifts_xml_dirs,
        sifts_download_dir, contact_cutoff, min_mapping_identity, logger)

Compute the contacts-table row of one protein in one PDB entry.
"""
function _protein_pdb_contacts_row(
        protein_disprot_id::AbstractString,
        protein_id::AbstractString,
        pdb::AbstractString,
        target_chains::AbstractVector{<:AbstractString},
        sequence::AbstractString;
        quality,
        mmcif_dirs,
        mmcif_download_dir::Union{Nothing, AbstractString},
        sifts_xml_dirs,
        sifts_download_dir::Union{Nothing, AbstractString},
        contact_cutoff::Real,
        min_mapping_identity::Real,
        logger::Function
)
    entry = _pdb_entry_id(pdb)
    row(status; chains = String[], partners = String[], homo = String[],
        n_mapped = 0, identity = NaN, resolved = Int[], contacts = Int[],
        n_unmapped = 0) = (
        protein_disprot_id = String(protein_disprot_id), protein_id = String(protein_id),
        PDB = entry, method = quality === nothing ? "missing" : quality.method,
        resolution = quality === nothing ? missing : quality.resolution,
        quality_category = string(_pdb_quality_category(quality)), status = status,
        chains = _optional_chain_list(chains),
        partner_chains = _optional_chain_list(partners),
        homo_chains_excluded = _optional_chain_list(homo), n_mapped = n_mapped,
        identity = identity, n_resolved = length(resolved),
        resolved_positions = _format_position_ranges(resolved),
        n_contacts = length(contacts), contact_positions = _format_position_ranges(contacts),
        n_contacts_unmapped = n_unmapped)

    mmcif = _mmcif_path(entry; mmcif_dirs = mmcif_dirs,
        download_dir = mmcif_download_dir, logger = logger)
    mmcif === nothing && return row("no_mmcif")
    sifts = _sifts_xml_path(entry; sifts_xml_dirs = sifts_xml_dirs,
        download_dir = sifts_download_dir, logger = logger)
    sifts === nothing && return row("no_sifts")

    residues = _load_first_model_heavy_atoms(mmcif)
    polymer = _polymer_chains(residues)
    chains = [chain for chain in target_chains if chain in polymer]
    isempty(chains) && return row("no_target_chain")

    chain_sequences = Dict(chain => _chain_sequence(residues, chain) for chain in polymer)
    homologous = _homologous_chains(chain_sequences, chains, sequence)
    partners = [chain for chain in polymer if !(chain in homologous)]
    homo_excluded = sort!([chain for chain in homologous if !(chain in chains)])
    partner_residues = [residue
                        for residue in residues
                        if residue.id.chain in partners && _is_polymer_residue(residue)]
    sifts_records = _read_sifts_residue_records(sifts)

    resolved, contacts = Set{Int}(), Set{Int}()
    n_mapped = n_identical = n_unmapped = 0
    for chain in chains
        query = [residue
                 for residue in residues
                 if residue.id.chain == chain && residue.id.name != "HOH"]
        mapping = _pdb_to_uniprot_positions(
            get(sifts_records, chain, SIFTSResidueRecord[]), protein_id)
        hits = isempty(partner_residues) ? Set{Int}() :
               _interface_residue_indices(query, partner_residues, contact_cutoff)
        chain_result = _chain_contact_positions(query, mapping, hits, sequence)
        union!(resolved, chain_result.resolved)
        union!(contacts, chain_result.contacts)
        n_mapped += chain_result.n_mapped
        n_identical += chain_result.n_identical
        n_unmapped += chain_result.n_unmapped_contacts
    end

    identity = n_mapped == 0 ? NaN : n_identical / n_mapped
    if n_mapped > 0 && identity < min_mapping_identity
        logger("$(protein_id)/$(entry): low PDB/UniProt identity ($(round(identity; digits = 2)))")
    end
    status = n_mapped == 0 ? "no_mapped_residue" :
             isempty(partners) ? "no_heteromeric_partner" : "ok"
    row(status; chains = chains, partners = partners, homo = homo_excluded,
        n_mapped = n_mapped, identity = identity, resolved = resolved,
        contacts = contacts, n_unmapped = n_unmapped)
end

"""
    build_pdb_contacts_table(dataset, disprot_pdb_mapping;
        mmcif_dirs = String[], mmcif_download_dir = nothing,
        sifts_xml_dirs = String[], sifts_download_dir = nothing,
        pdb_quality_table = Dict(), only_binding_proteins = true,
        checkpoint_dir = nothing, batch_size = DEFAULT_CHECKPOINT_BATCH_SIZE,
        contact_cutoff = PDB_CONTACT_CUTOFF_ANGSTROM,
        min_mapping_identity = PDB_CONTACT_MIN_MAPPING_IDENTITY,
        logger = println)

Compute heteromeric interface contacts for every (protein, PDB entry) pair of
`disprot_pdb_mapping` (output of `build_pdb_missing_residues_table`). With
`only_binding_proteins = true`, only proteins with at least one DisProt binding
residue (label `2`) are processed.

Returns `(table, failed_items)`, with one row per protein and PDB entry
(columns in `_PDB_CONTACTS_COLUMN_DEFINITIONS`). PDB method and resolution come
from `pdb_quality_table` (see `read_pdb_quality_table`); no quality filter is
applied here. Pairs are processed in parallel and can be checkpointed (see
`run_checkpointed_batches`).
"""
function build_pdb_contacts_table(
        dataset::AbstractDataFrame,
        disprot_pdb_mapping::AbstractDataFrame;
        mmcif_dirs = String[],
        mmcif_download_dir::Union{Nothing, AbstractString} = nothing,
        sifts_xml_dirs = String[],
        sifts_download_dir::Union{Nothing, AbstractString} = nothing,
        pdb_quality_table::AbstractDict = Dict{String, Any}(),
        only_binding_proteins::Bool = true,
        checkpoint_dir::Union{Nothing, AbstractString} = nothing,
        batch_size::Integer = DEFAULT_CHECKPOINT_BATCH_SIZE,
        contact_cutoff::Real = PDB_CONTACT_CUTOFF_ANGSTROM,
        min_mapping_identity::Real = PDB_CONTACT_MIN_MAPPING_IDENTITY,
        logger::Function = println
)
    _require_columns(
        dataset, [:protein_id, :protein_disprot_id, :sequence, :residue_labels],
        "dataset")
    _require_columns(disprot_pdb_mapping, [:protein_disprot_id, :protein_id, :PDB, :CHAIN],
        "disprot_pdb_mapping")

    labels_by_protein = _dataset_residue_labels_by_protein(dataset)
    sequence_by_protein = Dict(String(row.protein_id) => string(row.sequence)
    for row in eachrow(dataset) if !ismissing(row.sequence))

    pairs = Dict{String, NamedTuple}()
    for row in eachrow(disprot_pdb_mapping)
        protein_id = String(row.protein_id)
        haskey(sequence_by_protein, protein_id) || continue
        labels = get(labels_by_protein, protein_id, Int[])
        only_binding_proteins && !(_RESIDUE_LABEL_BINDING in labels) && continue
        entry = _pdb_entry_id(row.PDB)
        key = "$(protein_id)|$(entry)"
        previous = get(pairs, key, nothing)
        chains = _parse_chain_list(row.CHAIN)
        pairs[key] = (protein_disprot_id = String(row.protein_disprot_id),
            protein_id = protein_id, PDB = entry,
            chains = previous === nothing ? chains : unique(vcat(previous.chains, chains)))
    end
    items = sort!(collect(keys(pairs)))
    logger("Protein/PDB pairs to process: $(length(items)) ($(length(unique(p.protein_id for p in values(pairs)))) proteins)")

    result = run_checkpointed_batches(items, (contacts = _empty_pdb_contacts_table(),);
        checkpoint_dir = checkpoint_dir, batch_size = batch_size,
        logger = logger) do item
        pair = pairs[item]
        # `local`: the outer function also has a `table` variable, which a
        # closure would otherwise share between threads.
        local pair_table = _empty_pdb_contacts_table()
        push!(pair_table,
            _protein_pdb_contacts_row(pair.protein_disprot_id, pair.protein_id,
                pair.PDB, pair.chains, sequence_by_protein[pair.protein_id];
                quality = get(pdb_quality_table, pair.PDB, nothing),
                mmcif_dirs = mmcif_dirs, mmcif_download_dir = mmcif_download_dir,
                sifts_xml_dirs = sifts_xml_dirs, sifts_download_dir = sifts_download_dir,
                contact_cutoff = contact_cutoff, min_mapping_identity = min_mapping_identity,
                logger = logger))
        (contacts = pair_table,)
    end

    table = unique(result.tables.contacts, [:protein_id, :PDB])
    sort!(table, [:protein_disprot_id, :protein_id, :PDB])
    (table = table, failed_items = result.failed_items)
end

"""
    _PDB_CONTACTS_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_pdb_contacts_table`.
"""
const _PDB_CONTACTS_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "PDB" => "PDB entry code (lowercase)",
    "method" => "experimental method from the PDB quality table (`missing` if absent)",
    "resolution" => "resolution in Å from the PDB quality table (empty if absent or not applicable)",
    "quality_category" => "`resolved` (resolution known), `no_resolution` (e.g. NMR) or `absent` (not in the quality table)",
    "status" => "`ok`, `no_heteromeric_partner`, `no_mmcif`, `no_sifts`, `no_target_chain` or `no_mapped_residue` (see `PDB_CONTACT_STATUSES`)",
    "chains" => "protein chains of the PDB entry used (SIFTS chains present in the mmCIF, `;`-separated)",
    "partner_chains" => "heteromeric partner chains (polymer chains that are not copies of the protein)",
    "homo_chains_excluded" => "other copies of the protein excluded from partners (k-mer identity)",
    "n_mapped" => "observed protein residues mapped to the UniProt sequence (all protein chains)",
    "identity" => "fraction of mapped residues whose amino acid matches the UniProt sequence (`NaN` if none)",
    "n_resolved" => "number of distinct UniProt positions observed in the structure",
    "resolved_positions" => "UniProt positions observed in the structure, as ranges (`12-40;55`)",
    "n_contacts" => "number of distinct UniProt positions in contact with a partner chain",
    "contact_positions" => "UniProt positions in contact with a partner chain (heavy atoms within the cutoff), as ranges; all residues, not only disordered ones",
    "n_contacts_unmapped" => "residues in contact without UniProt mapping (e.g. tags)"
]
