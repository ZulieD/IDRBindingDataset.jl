using DataFrames

# ---------------------------------------------------------------------------
# Residue-level SIFTS XML mapping (PDB residue <-> UniProt residue).
#
# SIFTS XML files (one per PDB entry, `<pdb>.xml.gz`) map every residue of the
# deposited sequence to its UniProt position and flag residues without
# coordinates in the structure (`missing`). Source:
# https://www.ebi.ac.uk/pdbe/docs/sifts/ (accessed on 17-11-2026), read with
# `MIToS.SIFTS` (https://diegozea.github.io/MIToS.jl/stable/SIFTS/).
# ---------------------------------------------------------------------------
const SIFTS_XML_EXTENSIONS = (".xml.gz", ".xml")

"""
    SIFTSResidueRecord

Minimal, plain-data view of one `MIToS.SIFTS.SIFTSResidue`:
- `chain`: PDB author chain identifier
- `pdb_number`: PDB author residue number (with insertion code), `""` if absent
- `uniprot_id`: mapped UniProt accession, `""` if unmapped
- `uniprot_position`: mapped UniProt position, `0` if unmapped
- `is_missing`: `true` when the residue has no coordinates in the structure
"""
const SIFTSResidueRecord = NamedTuple{
    (:chain, :pdb_number, :uniprot_id, :uniprot_position, :is_missing),
    Tuple{String, String, String, Int, Bool}
}

"""
    _pdb_entry_id(value)::String

Normalize a PDB identifier to its lowercase 4-character entry code
(e.g. `"1ABC_A"` -> `"1abc"`).
"""
function _pdb_entry_id(value)::String
    text = lowercase(strip(string(value)))
    occursin(r"^[0-9][a-z0-9]{3}", text) ? text[1:4] : text
end

"""
    _find_pdb_entry_file(dirs, pdb, extensions)

Return the first existing `<pdb><extension>` file found in `dirs`, looking both
directly inside each directory and inside the PDB-style middle-letters
subdirectory (e.g. `ab/1abc.xml.gz`). Returns `nothing` when absent.
"""
function _find_pdb_entry_file(
        dirs,
        pdb::AbstractString,
        extensions
)::Union{Nothing, String}
    entry = _pdb_entry_id(pdb)
    for dir in dirs, extension in extensions

        for path in (joinpath(dir, entry * extension),
            joinpath(dir, entry[2:3], entry * extension))
            isfile(path) && return path
        end
    end
    nothing
end

"""
    _sifts_xml_path(pdb; sifts_xml_dirs, download_dir = nothing, logger = println)

Return the path of the SIFTS XML file of `pdb`. When it is not found in
`sifts_xml_dirs` and `download_dir` is given, the file is downloaded there with
`MIToS.SIFTS.downloadsifts` (atomic write). Returns `nothing` when the file is
unavailable.
"""
function _sifts_xml_path(
        pdb::AbstractString;
        sifts_xml_dirs,
        download_dir::Union{Nothing, AbstractString} = nothing,
        logger::Function = println
)::Union{Nothing, String}
    search_dirs = download_dir === nothing ? collect(sifts_xml_dirs) :
                  vcat(collect(sifts_xml_dirs), [download_dir])
    path = _find_pdb_entry_file(search_dirs, pdb, SIFTS_XML_EXTENSIONS)
    path === nothing || return path
    download_dir === nothing && return nothing

    entry = _pdb_entry_id(pdb)
    destination = joinpath(download_dir, entry * ".xml.gz")
    try
        mkpath(download_dir)
        tmp_path = tempname(download_dir) * ".xml.gz"
        MIToS.SIFTS.downloadsifts(entry; filename = tmp_path)
        mv(tmp_path, destination; force = true)
        return destination
    catch err
        logger("SIFTS XML download failed for $(entry): $(sprint(showerror, err))")
        return nothing
    end
end

"""
    _sifts_residue_record(residue)::SIFTSResidueRecord

Convert one `MIToS.SIFTS.SIFTSResidue` into a `SIFTSResidueRecord`.
"""
function _sifts_residue_record(residue::MIToS.SIFTS.SIFTSResidue)::SIFTSResidueRecord
    pdb = residue.PDB
    uniprot = residue.UniProt
    chain = ismissing(pdb) ? "" : pdb.chain
    pdb_number = ismissing(pdb) ? "" : pdb.number
    uniprot_id = ismissing(uniprot) ? "" : uppercase(uniprot.id)
    position = ismissing(uniprot) ? nothing : tryparse(Int, uniprot.number)
    (
        chain = chain,
        pdb_number = pdb_number,
        uniprot_id = position === nothing ? "" : uniprot_id,
        uniprot_position = position === nothing ? 0 : position,
        is_missing = residue.missing
    )
end

"""
    _read_sifts_residue_records(sifts_xml_path)

Parse one SIFTS XML file once and return its residues grouped by PDB author
chain: `Dict{String, Vector{SIFTSResidueRecord}}`.
"""
function _read_sifts_residue_records(
        sifts_xml_path::AbstractString
)::Dict{String, Vector{SIFTSResidueRecord}}
    residues = MIToS.SIFTS.read_file(sifts_xml_path, MIToS.SIFTS.SIFTSXML)
    by_chain = Dict{String, Vector{SIFTSResidueRecord}}()
    for residue in residues
        record = _sifts_residue_record(residue)
        isempty(record.chain) && continue
        push!(get!(by_chain, record.chain, SIFTSResidueRecord[]), record)
    end
    by_chain
end

"""
    _pdb_to_uniprot_positions(records, accession)

Map PDB author residue numbers of observed residues to UniProt positions of
`accession`: `Dict{String, Int}`.
"""
function _pdb_to_uniprot_positions(
        records::AbstractVector{SIFTSResidueRecord},
        accession::AbstractString
)::Dict{String, Int}
    target = uppercase(strip(accession))
    mapping = Dict{String, Int}()
    for record in records
        record.is_missing && continue
        isempty(record.pdb_number) && continue
        record.uniprot_id == target || continue
        mapping[record.pdb_number] = record.uniprot_position
    end
    mapping
end

"""
    _format_chain_list(chains)::String

Encode a list of PDB chain identifiers as one `;`-separated string.
"""
_format_chain_list(chains)::String = join((string(chain) for chain in chains), ";")

"""
    _parse_chain_list(value)::Vector{String}

Parse a chain list written by `_format_chain_list` (`"A;B"`). The legacy
`["A", "B"]` form is also accepted. Returns an empty vector for `missing`.
"""
function _parse_chain_list(value)::Vector{String}
    ismissing(value) && return String[]
    text = strip(string(value))
    quoted = [String(something(m.captures[1])) for m in eachmatch(r"\"([^\"]+)\"", text)]
    isempty(quoted) || return quoted
    [String(strip(token))
     for token in split(strip(text, ['[', ']']), r"[;,]")
     if !isempty(strip(token))]
end
