using CodecZlib
using DataFrames
using HTTP
using JSON

# ---------------------------------------------------------------------------
# UniProtKB "Binding site" features overlapping DisProt disordered residues.
#
# Entries are fetched in batches from the UniProt REST search endpoint
# (https://www.uniprot.org/help/api_queries, accessed on 17-11-2026).
# Feature coordinates are 1-based, inclusive positions on the canonical
# UniProt sequence, i.e. the dataset `sequence`.
# ---------------------------------------------------------------------------
const UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
const UNIPROT_RELEASE_HEADER = "X-UniProt-Release"
# Accessions per request: keeps the query URL well under server length limits.
const UNIPROT_DEFAULT_BATCH_SIZE = 90
# Maximum page size accepted by the UniProt REST API.
const UNIPROT_PAGE_SIZE = 500
const UNIPROT_DEFAULT_MAX_RETRIES = 3
const UNIPROT_DEFAULT_BASE_SLEEP = 0.3
const UNIPROT_BINDING_SITE_FEATURE_TYPE = "Binding site"
const UNIPROT_PUBMED_EVIDENCE_SOURCE = "PubMed"

"""
    UniProtBindingSite

One UniProt binding-site feature: `start`/`stop` positions (inclusive) and the
number of PubMed evidences supporting it (`n_pubmed`).
"""
const UniProtBindingSite = NamedTuple{(:start, :stop, :n_pubmed), Tuple{Int, Int, Int}}

"""
    _uniprot_accession_query(accessions)::String

Build the UniProt search query matching any accession of `accessions`.
"""
function _uniprot_accession_query(accessions::AbstractVector{<:AbstractString})::String
    join(("accession:$(accession)" for accession in accessions), " OR ")
end

"""
    _uniprot_response_text(body)::String

Return the response body as text, decompressing it when it still starts with
the gzip magic bytes (`0x1f 0x8b`).
"""
function _uniprot_response_text(body::AbstractVector{UInt8})::String
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return String(transcode(GzipDecompressor, Vector{UInt8}(body)))
    end
    String(copy(body))
end

"""
    _uniprot_next_page_url(response)

Return the `rel="next"` URL of the `Link` header of `response`, or `nothing`.
"""
function _uniprot_next_page_url(response::HTTP.Response)::Union{Nothing, String}
    link = HTTP.header(response, "Link", "")
    occursin("rel=\"next\"", link) || return nothing
    matched = match(r"<([^>]+)>\s*;\s*rel=\"next\"", link)
    matched === nothing ? nothing : String(something(matched.captures[1]))
end

"""
    _default_uniprot_get_request(url, query)

Run one UniProt GET request without throwing on non-2xx statuses.
"""
function _default_uniprot_get_request(url::AbstractString, query)::HTTP.Response
    HTTP.get(url, ["Accept-Encoding" => "identity"]; query = query,
        status_exception = false)
end

"""
    fetch_uniprot_entries(accessions;
        batch_size = UNIPROT_DEFAULT_BATCH_SIZE,
        max_retries = UNIPROT_DEFAULT_MAX_RETRIES,
        base_sleep = UNIPROT_DEFAULT_BASE_SLEEP,
        get_request = _default_uniprot_get_request,
        sleeper = sleep,
        logger = println
    )

Fetch UniProtKB JSON entries for `accessions` in batches.

Returns `(entries, release)` where `entries` is a `Dict` keyed by uppercase
primary accession and `release` is the UniProt release reported by the server
(`missing` if absent). Throws if a batch still fails after `max_retries`
attempts, so a network failure is never mistaken for "no binding site".
"""
function fetch_uniprot_entries(
        accessions::AbstractVector{<:AbstractString};
        batch_size::Integer = UNIPROT_DEFAULT_BATCH_SIZE,
        max_retries::Integer = UNIPROT_DEFAULT_MAX_RETRIES,
        base_sleep::Real = UNIPROT_DEFAULT_BASE_SLEEP,
        get_request::Function = _default_uniprot_get_request,
        sleeper::Function = sleep,
        logger::Function = println
)
    batch_size > 0 || throw(ArgumentError("`batch_size` must be greater than zero."))
    entries = Dict{String, Any}()
    release = missing
    query_accessions = sort!(unique(uppercase.(strip.(String.(accessions)))))
    n_batches = cld(length(query_accessions), batch_size)

    for (batch_index, batch) in enumerate(Iterators.partition(query_accessions, batch_size))
        logger("UniProt batch $(batch_index)/$(n_batches): $(length(batch)) accession(s)")
        url = UNIPROT_SEARCH_URL
        query = Dict("query" => _uniprot_accession_query(collect(batch)),
            "format" => "json", "size" => string(UNIPROT_PAGE_SIZE))
        while true
            response = _request_with_retries(() -> get_request(url, query);
                max_retries = max_retries, base_sleep = base_sleep, sleeper = sleeper,
                is_success_response = r -> r.status == 200,
                failure_prefix = "UniProt request failed")
            header_release = HTTP.header(response, UNIPROT_RELEASE_HEADER, "")
            isempty(header_release) || (release = header_release)

            payload = JSON.parse(_uniprot_response_text(response.body))
            for entry in get(payload, "results", Any[])
                entries[uppercase(string(entry["primaryAccession"]))] = entry
            end

            next_url = _uniprot_next_page_url(response)
            next_url === nothing && break
            url = next_url
            query = Dict{String, String}()  # already encoded in `next_url`
        end
        sleeper(base_sleep)
    end

    (entries = entries, release = release)
end

"""
    _uniprot_location_value(location, side)

Return the integer `value` of `location[side]` (`"start"` or `"end"`), or
`nothing` when it is unknown.
"""
function _uniprot_location_value(location, side::AbstractString)::Union{Nothing, Int}
    endpoint = get(location, side, nothing)
    endpoint isa AbstractDict || return nothing
    value = get(endpoint, "value", nothing)
    value isa Integer ? Int(value) : nothing
end

"""
    _extract_uniprot_binding_sites(entry)::Vector{UniProtBindingSite}

Extract `"Binding site"` features from one UniProtKB JSON entry. Features with
unknown boundaries are skipped. Returns an empty vector for `nothing`.
"""
function _extract_uniprot_binding_sites(entry)::Vector{UniProtBindingSite}
    sites = UniProtBindingSite[]
    entry === nothing && return sites
    for feature in get(entry, "features", Any[])
        get(feature, "type", "") == UNIPROT_BINDING_SITE_FEATURE_TYPE || continue
        location = get(feature, "location", nothing)
        location isa AbstractDict || continue
        start = _uniprot_location_value(location, "start")
        stop = _uniprot_location_value(location, "end")
        (start === nothing || stop === nothing) && continue
        n_pubmed = count(
            evidence -> get(evidence, "source", "") ==
                        UNIPROT_PUBMED_EVIDENCE_SOURCE,
            get(feature, "evidences", Any[]))
        push!(sites, (start = start, stop = stop, n_pubmed = n_pubmed))
    end
    sites
end

"""
    _site_overlaps_disorder(site, labels)::Bool

Return `true` when `site` covers at least one residue labelled `1` or `2`.
"""
function _site_overlaps_disorder(site::UniProtBindingSite, labels::AbstractVector{<:Integer})::Bool
    low, high = max(site.start, 1), min(site.stop, length(labels))
    any(
        position -> labels[position] in
                    (_RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING), low:high)
end

"""
    _format_residue_ranges(sites)::String

Encode site ranges as one `;`-separated string, e.g. `"12-15;30"`.
"""
function _format_residue_ranges(sites::AbstractVector{UniProtBindingSite})::String
    join(
        (site.start == site.stop ? string(site.start) : "$(site.start)-$(site.stop)"
        for site in sites), ";")
end

"""
    _parse_residue_ranges(value)::Vector{Int}

Parse ranges written by `_format_residue_ranges` (`"12-15;30"`, legacy
`["12-15", "30"]` is also accepted) into unique positions `[12, 13, 14, 15, 30]`.
"""
function _parse_residue_ranges(value)::Vector{Int}
    (ismissing(value) || isempty(strip(string(value)))) && return Int[]
    positions = Int[]
    text = String(string(value))
    for matched in eachmatch(r"(\d+)\s*-\s*(\d+)|(\d+)", text)
        range_start, range_stop, single = matched.captures
        if range_start !== nothing
            first_position = parse(Int, range_start)
            last_position = parse(Int, something(range_stop))
            append!(positions, min(first_position, last_position):max(first_position, last_position))
        else
            push!(positions, parse(Int, something(single)))
        end
    end
    unique!(positions)
end

"""
    build_uniprot_binding_sites_table(dataset, entries)::DataFrame

Build one row per dataset protein from fetched UniProt `entries` (see
`fetch_uniprot_entries`). Only binding sites overlapping a DisProt disordered
residue (label `1` or `2`) are kept.

Columns:
- `protein_disprot_id`: DisProt identifier
- `protein_id`: UniProt accession
- `uniprot_entry_found`: the accession was returned by UniProt
- `binding_in_uniprot`: at least one kept binding site
- `residues_binding_in_uniprot`: kept site ranges, e.g. `"12-15;30"` (`missing` if none)
- `publication_mention_binding`: PubMed evidence count per kept site, same order (`missing` if none)
"""
function build_uniprot_binding_sites_table(
        dataset::AbstractDataFrame,
        entries::AbstractDict
)::DataFrame
    _require_columns(dataset, [:protein_id, :protein_disprot_id, :residue_labels], "dataset")
    labels_by_protein = _dataset_residue_labels_by_protein(dataset)
    disprot_ids = Dict(String(row.protein_id) => String(row.protein_disprot_id)
    for row in eachrow(dataset))
    table = DataFrame(protein_disprot_id = String[], protein_id = String[],
        uniprot_entry_found = Bool[],
        binding_in_uniprot = Bool[],
        residues_binding_in_uniprot = Union{Missing, String}[],
        publication_mention_binding = Union{Missing, String}[])

    for protein_id in sort!(collect(keys(labels_by_protein)))
        entry = get(entries, uppercase(protein_id), nothing)
        labels = labels_by_protein[protein_id]
        kept = filter(site -> _site_overlaps_disorder(site, labels),
            _extract_uniprot_binding_sites(entry))
        push!(table,
            (disprot_ids[protein_id], protein_id, entry !== nothing, !isempty(kept),
                isempty(kept) ? missing : _format_residue_ranges(kept),
                isempty(kept) ? missing : join((site.n_pubmed for site in kept), ";")))
    end
    table
end

"""
    _UNIPROT_BINDING_SITES_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_uniprot_binding_sites_table`.
"""
const _UNIPROT_BINDING_SITES_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "uniprot_entry_found" => "`true` when the accession was returned by the UniProt REST API",
    "binding_in_uniprot" => "`true` when at least one UniProt `Binding site` feature overlaps a DisProt residue labelled `1` or `2`",
    "residues_binding_in_uniprot" => "ranges of the kept binding sites on the UniProt sequence, e.g. `12-15;30` (`missing` if none)",
    "publication_mention_binding" => "number of PubMed evidences per kept site, same order (`missing` if none)"
]
