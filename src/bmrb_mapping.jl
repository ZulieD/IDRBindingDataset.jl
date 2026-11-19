using DataFrames
using HTTP
using JSON

# ---------------------------------------------------------------------------
# DisProt proteins <-> BMRB NMR entries: raw mapping with experimental metadata.
#
# Source: BMRB REST API v2 (https://api.bmrb.io/v2, accessed on 19-11-2026),
# entries in NMR-STAR JSON form. For each dataset protein (UniProt accession),
# every BMRB entry listed by the UniProt -> BMRB mapping is downloaded (with a
# disk cache) and described per assigned chemical-shift list: matched entity
# (alignment to the UniProt sequence), free / bound / homo-oligomer state,
# partner molecules, sample conditions (pH, temperature, ionic strength),
# buffer, solvent, spectrometer and shift counts. Nothing is filtered here:
# comparing free and bound states is done in a later step.
# ---------------------------------------------------------------------------
const BMRB_API_URL = "https://api.bmrb.io/v2"
const BMRB_REQUEST_HEADERS = ["Application" => "IDRBindingDataset.jl"]
const BMRB_UNIPROT_MAPPING_FILENAME = "mapping_uniprot_bmrb.json"
const BMRB_DEFAULT_MAX_RETRIES = 4
# BMRB answers HTTP 403 when the request rate is too high; wait before retrying.
const BMRB_DEFAULT_BASE_SLEEP = 12.0
# Two polymer entities with at least this identity over this fraction of the
# shorter sequence are copies of the same protein (homo-oligomer).
const BMRB_HOMO_MIN_IDENTITY = 0.95
const BMRB_HOMO_MIN_COVERAGE = 0.8
# Shift-list atoms counted per residue (amide H, N and backbone carbons).
const BMRB_BACKBONE_ATOMS = ("H", "N", "CA", "CB", "C")
# Sample components that are not buffer (solvent, reference, preservative).
const BMRB_NON_BUFFER_COMPONENTS = ("h2o", "d2o", "water", "dss", "tsp", "dss-d6",
    "sodium azide", "nan3", "azide")
const BMRB_LIST_SEPARATOR = " | "

# ---------------------------------------------------------------------------
# Download with disk cache
# ---------------------------------------------------------------------------
"""
    _default_bmrb_get_request(url)

Run one BMRB GET request without throwing on non-2xx statuses.
"""
function _default_bmrb_get_request(url::AbstractString)::HTTP.Response
    HTTP.get(url, BMRB_REQUEST_HEADERS; status_exception = false, readtimeout = 180)
end

"""
    _cached_bmrb_json(cache_path, url; get_request, max_retries, base_sleep, sleeper)

Return the parsed JSON stored at `cache_path`, downloading it from `url` first
when the file does not exist yet (atomic write).
"""
function _cached_bmrb_json(
        cache_path::AbstractString,
        url::AbstractString;
        get_request::Function = _default_bmrb_get_request,
        max_retries::Integer = BMRB_DEFAULT_MAX_RETRIES,
        base_sleep::Real = BMRB_DEFAULT_BASE_SLEEP,
        sleeper::Function = sleep
)
    if !isfile(cache_path)
        response = _request_with_retries(() -> get_request(url);
            max_retries = max_retries, base_sleep = base_sleep, sleeper = sleeper,
            is_success_response = r -> r.status == 200,
            failure_prefix = "BMRB request failed for $(url)")
        mkpath(dirname(cache_path))
        tmp_path = tempname(dirname(cache_path))
        write(tmp_path, response.body)
        mv(tmp_path, cache_path; force = true)
    end
    JSON.parse(read(cache_path, String))
end

# ---------------------------------------------------------------------------
# UniProt -> BMRB mapping
# ---------------------------------------------------------------------------
"""
    _bmrb_values(value)::Vector{String}

Return `value` as a vector of strings (scalars become one-element vectors,
`nothing` is dropped).
"""
function _bmrb_values(value)::Vector{String}
    value === nothing && return String[]
    value isa AbstractVector || return [string(value)]
    [string(item) for item in value if item !== nothing]
end

"""
    parse_bmrb_uniprot_mapping(data)::Dict{String, Vector{String}}

Parse the BMRB `mappings/uniprot/bmrb` answer into `UniProt accession =>
BMRB entry IDs`. Both a dictionary (`accession => ids`) and a list of records
(keys containing `uniprot` and `bmrb`, or `[accession, id]` pairs) are
accepted. Isoform suffixes (`P12345-2`) are removed.
"""
function parse_bmrb_uniprot_mapping(data)::Dict{String, Vector{String}}
    mapping = Dict{String, Vector{String}}()
    add!(accession, bmrb_id) = push!(
        get!(mapping, uppercase(String(first(split(strip(accession), "-")))), String[]),
        strip(bmrb_id))
    if data isa AbstractDict
        for (accession, ids) in pairs(data), bmrb_id in _bmrb_values(ids)

            add!(string(accession), bmrb_id)
        end
    elseif data isa AbstractVector
        for record in data
            if record isa AbstractDict
                record_keys = collect(keys(record))
                uniprot_key = findfirst(k -> occursin("uniprot", lowercase(string(k))),
                    record_keys)
                bmrb_key = findfirst(k -> occursin("bmrb", lowercase(string(k))),
                    record_keys)
                (uniprot_key === nothing || bmrb_key === nothing) && continue
                for accession in _bmrb_values(record[record_keys[uniprot_key]]),
                    bmrb_id in _bmrb_values(record[record_keys[bmrb_key]])

                    add!(accession, bmrb_id)
                end
            elseif record isa AbstractVector && length(record) >= 2
                add!(string(record[1]), string(record[2]))
            end
        end
    end
    for ids in values(mapping)
        unique!(ids)
    end
    mapping
end

"""
    fetch_bmrb_uniprot_mapping(; cache_dir, get_request, max_retries, base_sleep, sleeper)

Download (once, cached in `cache_dir`) and parse the BMRB UniProt -> BMRB
mapping (`mappings/uniprot/bmrb`, all match types).
"""
function fetch_bmrb_uniprot_mapping(;
        cache_dir::AbstractString,
        get_request::Function = _default_bmrb_get_request,
        max_retries::Integer = BMRB_DEFAULT_MAX_RETRIES,
        base_sleep::Real = BMRB_DEFAULT_BASE_SLEEP,
        sleeper::Function = sleep
)::Dict{String, Vector{String}}
    data = _cached_bmrb_json(joinpath(cache_dir, BMRB_UNIPROT_MAPPING_FILENAME),
        "$(BMRB_API_URL)/mappings/uniprot/bmrb?format=json&match_type=all";
        get_request = get_request, max_retries = max_retries, base_sleep = base_sleep,
        sleeper = sleeper)
    parse_bmrb_uniprot_mapping(data)
end

# ---------------------------------------------------------------------------
# NMR-STAR JSON parsing
# ---------------------------------------------------------------------------
"""
    BMRBLoop / BMRBSaveframe

Plain views of one NMR-STAR loop (`category` without leading `_`, `table` with
one string column per tag) and one saveframe (`category`, `tags`, `loops`).
"""
const BMRBLoop = NamedTuple{(:category, :table), Tuple{String, DataFrame}}
const BMRBSaveframe = NamedTuple{(:category, :tags, :loops),
    Tuple{String, Dict{String, String}, Vector{BMRBLoop}}}

"""
    _bmrb_text(value)::String

NMR-STAR value as text; JSON `null` becomes the NMR-STAR null `"."`.
"""
_bmrb_text(value)::String = value === nothing ? "." : string(value)

"""
    _bmrb_is_null(value)::Bool

`true` for NMR-STAR null values (`"."`, `"?"`, `""`).
"""
_bmrb_is_null(value)::Bool = strip(value) in (".", "?", "")

"""
    _bmrb_value(text)

`nothing` for NMR-STAR nulls, otherwise the stripped text.
"""
_bmrb_value(text::AbstractString) = _bmrb_is_null(text) ? nothing : String(strip(text))

"""
    _bmrb_float(text)::Union{Missing, Float64}

Parse a numeric NMR-STAR value; `missing` for nulls and non-numeric text.
"""
function _bmrb_float(text)::Union{Missing, Float64}
    text === nothing && return missing
    value = _bmrb_value(string(text))
    value === nothing && return missing
    parsed = tryparse(Float64, value)
    parsed === nothing ? missing : parsed
end

"""
    _parse_bmrb_loop(loop)::BMRBLoop

Parse one JSON loop (`category`, `tags`, `data`).
"""
function _parse_bmrb_loop(loop)::BMRBLoop
    category = String(lstrip(string(get(loop, "category", "")), '_'))
    tags = [String(tag) for tag in get(loop, "tags", Any[])]
    rows = get(loop, "data", Any[])
    table = DataFrame([tag => String[_bmrb_text(row[index]) for row in rows]
                       for (index, tag) in enumerate(tags)])
    (category = category, table = table)
end

"""
    parse_bmrb_entry(json)::Vector{BMRBSaveframe}

Parse one BMRB API v2 entry (`{"<id>": {"saveframes": [...]}}` or the inner
object) into saveframes.
"""
function parse_bmrb_entry(json)::Vector{BMRBSaveframe}
    root = haskey(json, "saveframes") ? json : first(values(json))
    saveframes = BMRBSaveframe[]
    for saveframe in get(root, "saveframes", Any[])
        tags = Dict{String, String}()
        for tag in get(saveframe, "tags", Any[])
            tags[String(tag[1])] = _bmrb_text(tag[2])
        end
        loops = BMRBLoop[_parse_bmrb_loop(loop) for loop in get(saveframe, "loops", Any[])]
        push!(saveframes, (category = string(get(saveframe, "category", "")),
            tags = tags, loops = loops))
    end
    saveframes
end

"""
    _bmrb_tag(saveframe, tag)

Value of `tag` in `saveframe`, or `nothing` when absent or null.
"""
function _bmrb_tag(saveframe::BMRBSaveframe, tag::AbstractString)
    _bmrb_value(get(
        saveframe.tags, tag, "."))
end

"""
    _bmrb_cell(row, column)

Value of `column` in a loop `row`, or `nothing` when absent or null.
"""
function _bmrb_cell(row, column::AbstractString)
    column in names(parent(row)) || return nothing
    _bmrb_value(row[column])
end

_bmrb_saveframes(saveframes, category) = [s for s in saveframes if s.category == category]

function _bmrb_loops(saveframes, category)
    [loop for s in saveframes for loop in s.loops
     if loop.category == category]
end

# ---------------------------------------------------------------------------
# Entities and sequence alignment
# ---------------------------------------------------------------------------
"""
    _local_alignment(query, reference; match = 2, mismatch = -1, gap = -2)

Smith-Waterman local alignment of `query` on `reference`. Returns
`(score, identity, n_aligned, indels, mapping)` where `mapping` maps aligned
`query` positions to `reference` positions (1-based) and `indels` counts
skipped residues inside the aligned region, on both sequences.
"""
function _local_alignment(
        query::AbstractString,
        reference::AbstractString;
        match::Integer = 2,
        mismatch::Integer = -1,
        gap::Integer = -2
)
    a, b = collect(query), collect(reference)
    n, m = length(a), length(b)
    scores = zeros(Int, n + 1, m + 1)
    trace = zeros(UInt8, n + 1, m + 1)   # 0 stop, 1 diagonal, 2 up, 3 left
    best = best_i = best_j = 0
    for j in 1:m, i in 1:n

        diagonal = scores[i, j] + (a[i] == b[j] ? match : mismatch)
        up = scores[i, j + 1] + gap
        left = scores[i + 1, j] + gap
        score = max(0, diagonal, up, left)
        scores[i + 1, j + 1] = score
        trace[i + 1, j + 1] = score == 0 ? 0x00 :
                              score == diagonal ? 0x01 :
                              score == up ? 0x02 : 0x03
        if score > best
            best, best_i, best_j = score, i, j
        end
    end
    mapping = Dict{Int, Int}()
    n_aligned = n_identical = 0
    i, j = best_i, best_j
    while i > 0 && j > 0 && trace[i + 1, j + 1] != 0x00
        step = trace[i + 1, j + 1]
        if step == 0x01
            mapping[i] = j
            n_aligned += 1
            n_identical += a[i] == b[j]
            i -= 1
            j -= 1
        elseif step == 0x02
            i -= 1
        else
            j -= 1
        end
    end
    indels = n_aligned == 0 ? 0 :
             (maximum(keys(mapping)) - minimum(keys(mapping)) + 1 - n_aligned) +
             (maximum(values(mapping)) - minimum(values(mapping)) + 1 - n_aligned)
    (score = best, identity = n_aligned > 0 ? n_identical / n_aligned : 0.0,
        n_aligned = n_aligned, indels = indels, mapping = mapping)
end

"""
    _clean_bmrb_sequence(text)::String

One-letter polymer sequence from NMR-STAR (`(XYZ)` modified residues become
`X`; whitespace and other symbols are removed).
"""
function _clean_bmrb_sequence(text)::String
    text === nothing && return ""
    uppercase(replace(replace(string(text), r"\([^)]*\)" => "X"), r"[^A-Za-z]" => ""))
end

"""
    _bmrb_entities(saveframes)::DataFrame

One row per `entity` saveframe: `entity_id`, `name`, `type` (`polymer`,
`non-polymer`, ...), `polymer_type`, `sequence` and `uniprot_ids` (UniProt
cross-references of the entity, `;`-separated, empty when none).
"""
function _bmrb_entities(saveframes)::DataFrame
    entities = DataFrame(entity_id = String[], name = String[], type = String[],
        polymer_type = String[], sequence = String[], uniprot_ids = String[])
    for saveframe in _bmrb_saveframes(saveframes, "entity")
        uniprot_ids = String[]
        for loop in saveframe.loops
            loop.category == "Entity_db_link" || continue
            for row in eachrow(loop.table)
                database = lowercase(something(_bmrb_cell(row, "Database_code"), ""))
                accession = _bmrb_cell(row, "Accession_code")
                (accession !== nothing && database in ("unp", "sp", "swiss-prot",
                    "uniprot", "uniprotkb", "trembl")) || continue
                push!(uniprot_ids, accession)
            end
        end
        push!(entities,
            (something(_bmrb_tag(saveframe, "ID"), "?"),
                something(_bmrb_tag(saveframe, "Name"), ""),
                lowercase(something(_bmrb_tag(saveframe, "Type"), "")),
                something(_bmrb_tag(saveframe, "Polymer_type"), ""),
                _clean_bmrb_sequence(_bmrb_tag(saveframe, "Polymer_seq_one_letter_code")),
                join(unique(uniprot_ids), ";")))
    end
    entities
end

"""
    _is_bmrb_polypeptide(entity)::Bool

`true` for polypeptide polymer entities with at least 5 residues.
"""
function _is_bmrb_polypeptide(entity)::Bool
    entity.type == "polymer" && length(entity.sequence) >= 5 &&
        (isempty(entity.polymer_type) ||
         occursin("polypeptide", lowercase(entity.polymer_type)))
end

"""
    _bmrb_target_entity(entities, protein_sequence)

Polypeptide entity with the best local-alignment score on `protein_sequence`.
Returns `(row_index, alignment)` or `nothing` when there is no polypeptide.
"""
function _bmrb_target_entity(entities::AbstractDataFrame, protein_sequence::AbstractString)
    best = nothing
    for (index, entity) in enumerate(eachrow(entities))
        _is_bmrb_polypeptide(entity) || continue
        alignment = _local_alignment(entity.sequence, protein_sequence)
        if best === nothing || alignment.score > last(best).score
            best = (index, alignment)
        end
    end
    best
end

"""
    _bmrb_assembly(saveframes, entities, target_index)

Molecules of the sample assembly and the state of the target entity:
`(state, n_target_copies, partners)` where `state` is `free` (only the target),
`bound` (other molecules present), or `homo_oligomer` (several copies of the
target protein). `partners` is a `DataFrame` of the other entities.
"""
function _bmrb_assembly(saveframes, entities::AbstractDataFrame, target_index::Integer)
    target = entities[target_index, :]
    assembly_loops = _bmrb_loops(saveframes, "Entity_assembly")
    assembly_ids = if !isempty(assembly_loops) &&
                      "Entity_ID" in names(first(assembly_loops).table)
        String.(first(assembly_loops).table.Entity_ID)
    else
        String.(entities.entity_id)
    end

    n_copies = count(==(target.entity_id), assembly_ids)
    partners = entities[1:0, :]
    for entity in eachrow(entities)
        entity.entity_id == target.entity_id && continue
        entity.entity_id in assembly_ids || continue
        if _is_bmrb_polypeptide(entity)
            alignment = _local_alignment(entity.sequence, target.sequence)
            if alignment.identity >= BMRB_HOMO_MIN_IDENTITY &&
               alignment.n_aligned >=
               BMRB_HOMO_MIN_COVERAGE *
               min(length(entity.sequence), length(target.sequence))
                n_copies += count(==(entity.entity_id), assembly_ids)
                continue
            end
        end
        push!(partners, entity)
    end
    state = n_copies > 1 ? "homo_oligomer" : isempty(partners) ? "free" : "bound"
    (state = state, n_target_copies = max(n_copies, 1), partners = partners)
end

# ---------------------------------------------------------------------------
# Experimental metadata
# ---------------------------------------------------------------------------
"""
    _bmrb_temperature_kelvin(value, units)

Temperature in K: values in °C (units containing `c`, or values below 150) are
converted.
"""
function _bmrb_temperature_kelvin(value::Union{Missing, Float64}, units::AbstractString)
    ismissing(value) && return missing
    celsius = (occursin("c", units) && !occursin("k", units)) || value < 150
    celsius ? value + 273.15 : value
end

"""
    _bmrb_sample_conditions(saveframes)

`Dict(sample_condition_list_id => (pH, temperature_K, ionic_strength_mM))`.
"""
function _bmrb_sample_conditions(saveframes)
    conditions = Dict{String, NamedTuple}()
    for saveframe in _bmrb_saveframes(saveframes, "sample_conditions")
        pH = temperature = ionic = missing
        for loop in saveframe.loops
            loop.category == "Sample_condition_variable" || continue
            for row in eachrow(loop.table)
                kind = lowercase(something(_bmrb_cell(row, "Type"), ""))
                value = _bmrb_float(_bmrb_cell(row, "Val"))
                units = lowercase(something(_bmrb_cell(row, "Val_units"), ""))
                if startswith(kind, "ph")
                    pH = value
                elseif startswith(kind, "temp")
                    temperature = _bmrb_temperature_kelvin(value, units)
                elseif occursin("ionic", kind) && !ismissing(value)
                    ionic = units == "m" ? 1000 * value : value
                end
            end
        end
        conditions[something(_bmrb_tag(saveframe, "ID"), "?")] = (pH = pH,
            temperature_K = temperature, ionic_strength_mM = ionic)
    end
    conditions
end

"""
    _bmrb_entry_metadata(saveframes, entities)

Spectrometers, experiments, samples (solvent, buffer components) and spectral
widths of one entry, indexed by their NMR-STAR IDs.
"""
function _bmrb_entry_metadata(saveframes, entities::AbstractDataFrame)
    spectrometers = Dict{String, NamedTuple}()
    for saveframe in _bmrb_saveframes(saveframes, "NMR_spectrometer")
        model = strip(join(
            filter(!isnothing, [_bmrb_tag(saveframe, "Manufacturer"),
                _bmrb_tag(saveframe, "Model")]),
            " "))
        spectrometers[something(_bmrb_tag(saveframe, "ID"), "?")] = (model = model,
            field_MHz = _bmrb_float(_bmrb_tag(saveframe, "Field_strength")))
    end

    experiments = Dict{String, NamedTuple}()
    for loop in _bmrb_loops(saveframes, "Experiment"), row in eachrow(loop.table)

        id = _bmrb_cell(row, "ID")
        id === nothing && continue
        experiments[id] = (name = something(_bmrb_cell(row, "Name"), ""),
            sample = _bmrb_cell(row, "Sample_ID"),
            spectrometer = _bmrb_cell(row, "NMR_spectrometer_ID"))
    end

    entity_names = filter(!isempty, lowercase.(entities.name))
    samples = Dict{String, NamedTuple}()
    for saveframe in _bmrb_saveframes(saveframes, "sample")
        buffer = String[]
        for loop in saveframe.loops
            loop.category == "Sample_component" || continue
            for row in eachrow(loop.table)
                # Studied molecules (protein, partner) are not buffer components.
                (_bmrb_cell(row, "Entity_ID") === nothing &&
                 _bmrb_cell(row, "Assembly_ID") === nothing) || continue
                name = _bmrb_cell(row, "Mol_common_name")
                name === nothing && continue
                lower = lowercase(name)
                lower in BMRB_NON_BUFFER_COMPONENTS && continue
                any(e -> occursin(e, lower) || occursin(lower, e), entity_names) && continue
                concentration = something(_bmrb_cell(row, "Concentration_val"), "?")
                units = something(_bmrb_cell(row, "Concentration_val_units"), "")
                push!(buffer, strip("$(name) $(concentration) $(units)"))
            end
        end
        samples[something(_bmrb_tag(saveframe, "ID"), "?")] = (
            solvent = something(_bmrb_tag(saveframe, "Solvent_system"), ""),
            buffer = sort!(buffer))
    end

    n15_widths = Dict{String, Float64}()   # Experiment_ID => 15N sweep width (ppm)
    for saveframe in _bmrb_saveframes(saveframes, "spectral_peak_list")
        experiment = something(_bmrb_tag(saveframe, "Experiment_ID"), "?")
        for loop in saveframe.loops
            loop.category == "Spectral_dim" || continue
            for row in eachrow(loop.table)
                uppercase(something(_bmrb_cell(row, "Atom_type"), "")) == "N" || continue
                width = _bmrb_float(_bmrb_cell(row, "Sweep_width"))
                ismissing(width) && continue
                if lowercase(something(_bmrb_cell(row, "Sweep_width_units"), "ppm")) == "hz"
                    frequency = _bmrb_float(_bmrb_cell(row, "Spectrometer_frequency"))
                    ismissing(frequency) && continue
                    width /= frequency   # Hz / MHz = ppm
                end
                n15_widths[experiment] = width
            end
        end
    end
    (spectrometers = spectrometers, experiments = experiments, samples = samples,
        n15_widths = n15_widths)
end

"""
    _bmrb_shift_counts(table, target_entity_id, single_polymer)

Count the shifts of the target entity in one `Atom_chem_shift` loop: all
atoms, residues with at least one shift, and backbone atoms
(`BMRB_BACKBONE_ATOMS`, `HN` counted as `H`).
"""
function _bmrb_shift_counts(table::AbstractDataFrame, target_entity_id::AbstractString,
        single_polymer::Bool)
    has_entity = "Entity_ID" in names(table)
    residue_column = "Seq_ID" in names(table) ? "Seq_ID" :
                     "Comp_index_ID" in names(table) ? "Comp_index_ID" : nothing
    atom_counts = Dict(atom => 0 for atom in BMRB_BACKBONE_ATOMS)
    residues = Set{String}()
    n_shifts = 0
    for row in eachrow(table)
        entity = has_entity ? _bmrb_value(row.Entity_ID) : nothing
        keep = !has_entity || entity == target_entity_id ||
               (single_polymer && entity === nothing)
        keep || continue
        n_shifts += 1
        residue_column === nothing || push!(residues, row[residue_column])
        atom = uppercase(something(_bmrb_cell(row, "Atom_ID"), ""))
        atom = atom == "HN" ? "H" : atom
        haskey(atom_counts, atom) && (atom_counts[atom] += 1)
    end
    (n_shifts = n_shifts, n_residues = length(residues), atoms = atom_counts)
end

# ---------------------------------------------------------------------------
# Rows of the mapping table
# ---------------------------------------------------------------------------
"""
    _empty_bmrb_mapping_table()

Return the empty mapping table (schema), see `_BMRB_MAPPING_COLUMN_DEFINITIONS`.
"""
function _empty_bmrb_mapping_table()::DataFrame
    optional_string = Union{Missing, String}
    optional_int = Union{Missing, Int}
    optional_float = Union{Missing, Float64}
    DataFrame(protein_disprot_id = String[], protein_id = String[], bmrb_id = String[],
        status = String[], entry_title = optional_string[],
        target_entity_id = optional_string[], target_entity_name = optional_string[],
        entity_identity = optional_float[], entity_n_aligned = optional_int[],
        entity_indels = optional_int[], cover_start = optional_int[],
        cover_end = optional_int[], state = optional_string[],
        n_target_copies = optional_int[], n_partners = optional_int[],
        partner_names = optional_string[], partner_types = optional_string[],
        partner_polymer_types = optional_string[], partner_lengths = optional_string[],
        partner_uniprot_ids = optional_string[], partner_sequences = optional_string[],
        shift_list_id = optional_string[], sample_condition_list_id = optional_string[],
        pH = optional_float[], temperature_K = optional_float[],
        ionic_strength_mM = optional_float[], solvent = optional_string[],
        buffer = optional_string[], experiments = optional_string[],
        spectrometers = optional_string[], field_strengths_MHz = optional_string[],
        n15_sweep_width_ppm = optional_float[], n_shifts = optional_int[],
        n_residues_with_shifts = optional_int[], n_H = optional_int[], n_N = optional_int[],
        n_CA = optional_int[], n_CB = optional_int[], n_C = optional_int[])
end

"""
    _optional_join(values)::Union{Missing, String}

Join with `BMRB_LIST_SEPARATOR`, `missing` for an empty list.
"""
_optional_join(values)::Union{Missing, String} = isempty(values) ? missing :
                                                 join(values, BMRB_LIST_SEPARATOR)

"""
    _bmrb_entry_rows(protein_disprot_id, protein_id, protein_sequence, bmrb_id, json)

Rows of the mapping table for one protein and one parsed BMRB entry: one row per
assigned chemical-shift list (one row with `status = no_shift_list` when the
entry has none, `no_polypeptide` when no entity can be the protein).
"""
function _bmrb_entry_rows(
        protein_disprot_id::AbstractString,
        protein_id::AbstractString,
        protein_sequence::AbstractString,
        bmrb_id::AbstractString,
        json
)::DataFrame
    table = _empty_bmrb_mapping_table()
    saveframes = parse_bmrb_entry(json)
    information = _bmrb_saveframes(saveframes, "entry_information")
    title = isempty(information) ? missing :
            something(_bmrb_tag(first(information), "Title"), missing)
    title = ismissing(title) ? missing : replace(title, r"\s+" => " ")
    base = Dict{Symbol, Any}(:protein_disprot_id => String(protein_disprot_id),
        :protein_id => String(protein_id), :bmrb_id => String(bmrb_id),
        :entry_title => title)

    entities = _bmrb_entities(saveframes)
    target = _bmrb_target_entity(entities, protein_sequence)
    if target === nothing
        push!(table, merge(base, Dict(:status => "no_polypeptide")); cols = :subset)
        return table
    end
    target_index, alignment = target
    target_entity = entities[target_index, :]
    covered = collect(values(alignment.mapping))
    assembly = _bmrb_assembly(saveframes, entities, target_index)
    partners = assembly.partners
    merge!(base,
        Dict(:target_entity_id => target_entity.entity_id,
            :target_entity_name => target_entity.name,
            :entity_identity => round(alignment.identity; digits = 3),
            :entity_n_aligned => alignment.n_aligned, :entity_indels => alignment.indels,
            :cover_start => isempty(covered) ? missing : minimum(covered),
            :cover_end => isempty(covered) ? missing : maximum(covered),
            :state => assembly.state, :n_target_copies => assembly.n_target_copies,
            :n_partners => nrow(partners),
            :partner_names => _optional_join(partners.name),
            :partner_types => _optional_join(partners.type),
            :partner_polymer_types => _optional_join(partners.polymer_type),
            :partner_lengths => _optional_join(string.(length.(partners.sequence))),
            :partner_uniprot_ids => _optional_join(partners.uniprot_ids),
            :partner_sequences => _optional_join(partners.sequence)))

    shift_lists = _bmrb_saveframes(saveframes, "assigned_chemical_shifts")
    if isempty(shift_lists)
        push!(table, merge(base, Dict(:status => "no_shift_list")); cols = :subset)
        return table
    end

    conditions = _bmrb_sample_conditions(saveframes)
    metadata = _bmrb_entry_metadata(saveframes, entities)
    single_polymer = count(==("polymer"), entities.type) == 1
    for shift_list in shift_lists
        condition_id = _bmrb_tag(shift_list, "Sample_condition_list_ID")
        condition = get(conditions, something(condition_id, "?"),
            (pH = missing, temperature_K = missing, ionic_strength_mM = missing))

        experiment_ids, sample_ids = String[], String[]
        for loop in shift_list.loops
            loop.category == "Chem_shift_experiment" || continue
            for row in eachrow(loop.table)
                experiment = _bmrb_cell(row, "Experiment_ID")
                experiment === nothing || push!(experiment_ids, experiment)
                sample = _bmrb_cell(row, "Sample_ID")
                sample === nothing || push!(sample_ids, sample)
            end
        end
        # Without explicit links, the whole entry describes the list.
        isempty(experiment_ids) &&
            (experiment_ids = sort!(collect(keys(metadata.experiments))))
        for experiment in experiment_ids
            sample = get(metadata.experiments, experiment, (sample = nothing,)).sample
            sample === nothing || push!(sample_ids, sample)
        end
        isempty(sample_ids) && (sample_ids = sort!(collect(keys(metadata.samples))))
        unique!(sample_ids)
        spectrometer_ids = unique([metadata.experiments[e].spectrometer
                                   for e in experiment_ids
                                   if haskey(metadata.experiments, e) &&
            metadata.experiments[e].spectrometer !== nothing])
        isempty(spectrometer_ids) &&
            (spectrometer_ids = sort!(collect(keys(metadata.spectrometers))))
        spectrometers = [metadata.spectrometers[i]
                         for i in spectrometer_ids
                         if haskey(metadata.spectrometers, i)]
        samples = [metadata.samples[i] for i in sample_ids if haskey(metadata.samples, i)]
        widths = [metadata.n15_widths[e]
                  for e in experiment_ids
                  if haskey(metadata.n15_widths, e)]

        counts = (n_shifts = 0, n_residues = 0,
            atoms = Dict(atom => 0 for atom in BMRB_BACKBONE_ATOMS))
        for loop in shift_list.loops
            loop.category == "Atom_chem_shift" || continue
            loop_counts = _bmrb_shift_counts(loop.table, target_entity.entity_id,
                single_polymer)
            counts = (n_shifts = counts.n_shifts + loop_counts.n_shifts,
                n_residues = counts.n_residues + loop_counts.n_residues,
                atoms = mergewith(+, counts.atoms, loop_counts.atoms))
        end

        push!(table,
            merge(base,
                Dict(:status => "ok",
                    :shift_list_id => something(_bmrb_tag(shift_list, "ID"), "?"),
                    :sample_condition_list_id => something(condition_id, missing),
                    :pH => condition.pH, :temperature_K => condition.temperature_K,
                    :ionic_strength_mM => condition.ionic_strength_mM,
                    :solvent => _optional_join(unique(filter(!isempty,
                        [s.solvent for s in samples]))),
                    :buffer =>
                        _optional_join(unique(reduce(vcat, [s.buffer for s in samples];
                            init = String[]))),
                    :experiments => _optional_join(unique(filter(!isempty,
                        [get(metadata.experiments, e, (name = "",)).name
                         for e in experiment_ids]))),
                    :spectrometers => _optional_join(unique(filter(!isempty,
                        [s.model for s in spectrometers]))),
                    :field_strengths_MHz => _optional_join(unique(string.(
                        [s.field_MHz for s in spectrometers if !ismissing(s.field_MHz)]))),
                    :n15_sweep_width_ppm => isempty(widths) ? missing : minimum(widths),
                    :n_shifts => counts.n_shifts, :n_residues_with_shifts =>
                        counts.n_residues,
                    :n_H => counts.atoms["H"], :n_N => counts.atoms["N"],
                    :n_CA => counts.atoms["CA"], :n_CB => counts.atoms["CB"],
                    :n_C => counts.atoms["C"]));
            cols = :subset)
    end
    table
end

"""
    build_disprot_bmrb_mapping_table(dataset, uniprot_bmrb_mapping;
        cache_dir, get_request = _default_bmrb_get_request,
        max_retries = BMRB_DEFAULT_MAX_RETRIES, base_sleep = BMRB_DEFAULT_BASE_SLEEP,
        sleeper = sleep, logger = println)

Describe every BMRB entry mapped to a dataset protein (`protein_id` =>
`uniprot_bmrb_mapping`, see `parse_bmrb_uniprot_mapping`), one row per
assigned chemical-shift list. Entries are downloaded once into `cache_dir`
(`entry_<id>.json`), so a new run only downloads missing entries.

Returns `(table, failed_entries)`; an entry that cannot be downloaded or parsed
gets one row with `status = error` and is listed in `failed_entries`.
"""
function build_disprot_bmrb_mapping_table(
        dataset::AbstractDataFrame,
        uniprot_bmrb_mapping::AbstractDict;
        cache_dir::AbstractString,
        get_request::Function = _default_bmrb_get_request,
        max_retries::Integer = BMRB_DEFAULT_MAX_RETRIES,
        base_sleep::Real = BMRB_DEFAULT_BASE_SLEEP,
        sleeper::Function = sleep,
        logger::Function = println
)
    _require_columns(dataset, [:protein_id, :protein_disprot_id, :sequence], "dataset")
    table = _empty_bmrb_mapping_table()
    failed_entries = String[]
    proteins = sort!(
        [row
         for row in eachrow(dataset)
         if haskey(uniprot_bmrb_mapping, uppercase(String(row.protein_id)))];
        by = row -> String(row.protein_id))
    n_entries = sum(
        (length(uniprot_bmrb_mapping[uppercase(String(row.protein_id))])
        for row in proteins);
        init = 0)
    logger("Dataset proteins with BMRB entries: $(length(proteins)) ($(n_entries) protein/entry pairs)")

    for (index, row) in enumerate(proteins)
        protein_id = String(row.protein_id)
        for bmrb_id in sort(uniprot_bmrb_mapping[uppercase(protein_id)])
            try
                json = _cached_bmrb_json(joinpath(cache_dir, "entry_$(bmrb_id).json"),
                    "$(BMRB_API_URL)/entry/$(bmrb_id)"; get_request = get_request,
                    max_retries = max_retries, base_sleep = base_sleep, sleeper = sleeper)
                append!(table,
                    _bmrb_entry_rows(String(row.protein_disprot_id), protein_id,
                        uppercase(string(row.sequence)), bmrb_id, json);
                    cols = :subset)
            catch err
                push!(failed_entries, "$(protein_id)/$(bmrb_id)")
                logger("BMRB entry $(bmrb_id) ($(protein_id)) failed: $(sprint(showerror, err))")
                push!(table,
                    Dict(:protein_disprot_id => String(row.protein_disprot_id),
                        :protein_id => protein_id, :bmrb_id => bmrb_id, :status => "error");
                    cols = :subset)
            end
        end
        index % 50 == 0 && logger("  $(index)/$(length(proteins)) proteins done")
    end

    sort!(table, [:protein_disprot_id, :protein_id, :bmrb_id, :shift_list_id])
    (table = table, failed_entries = failed_entries)
end

"""
    _BMRB_MAPPING_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_disprot_bmrb_mapping_table`.
Lists use `" | "` as separator, in the same order across partner columns.
"""
const _BMRB_MAPPING_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "bmrb_id" => "BMRB entry ID mapped to `protein_id` by the BMRB UniProt mapping",
    "status" => "`ok` (one row per assigned chemical-shift list), `no_shift_list`, `no_polypeptide` (no entity can be the protein) or `error` (entry not downloaded/parsed)",
    "entry_title" => "BMRB entry title",
    "target_entity_id" => "BMRB entity aligned best to the UniProt sequence (local alignment)",
    "target_entity_name" => "name of that entity",
    "entity_identity" => "identity of the entity/UniProt alignment (1.0 = identical on the aligned region)",
    "entity_n_aligned" => "number of aligned residues",
    "entity_indels" => "insertions/deletions inside the aligned region (both sequences)",
    "cover_start" => "first UniProt position covered by the entity",
    "cover_end" => "last UniProt position covered by the entity",
    "state" => "`free` (protein alone), `bound` (other molecules in the sample assembly) or `homo_oligomer` (several copies of the protein)",
    "n_target_copies" => "copies of the protein in the assembly (same entity or ≥ 95 % identical entities)",
    "n_partners" => "number of other entities in the assembly (proteins, nucleic acids, ligands, ions)",
    "partner_names" => "partner entity names",
    "partner_types" => "partner entity types (`polymer`, `non-polymer`, ...)",
    "partner_polymer_types" => "partner polymer types (e.g. `polypeptide(L)`, `polydeoxyribonucleotide`)",
    "partner_lengths" => "partner sequence lengths (`0` for non-polymers)",
    "partner_uniprot_ids" => "UniProt cross-references of each partner (`;` within one partner)",
    "partner_sequences" => "partner one-letter sequences",
    "shift_list_id" => "ID of the assigned chemical-shift list in the entry",
    "sample_condition_list_id" => "sample-condition list used by the shift list",
    "pH" => "pH of the sample conditions",
    "temperature_K" => "temperature in K (°C values converted)",
    "ionic_strength_mM" => "ionic strength in mM",
    "solvent" => "solvent system of the sample(s), e.g. `90% H2O/10% D2O`",
    "buffer" => "buffer components with concentrations and units (studied molecules, water, DSS/TSP and azide excluded)",
    "experiments" => "experiments linked to the shift list (e.g. `2D 1H-15N HSQC`)",
    "spectrometers" => "spectrometer manufacturer and model",
    "field_strengths_MHz" => "spectrometer field strengths (MHz)",
    "n15_sweep_width_ppm" => "15N spectral width (ppm) of the linked experiments, when deposited",
    "n_shifts" => "chemical shifts of the protein entity in the list",
    "n_residues_with_shifts" => "residues of the protein entity with at least one shift",
    "n_H" => "amide 1H shifts",
    "n_N" => "amide 15N shifts",
    "n_CA" => "13Cα shifts",
    "n_CB" => "13Cβ shifts",
    "n_C" => "carbonyl 13C shifts"
]
