using DataFrames

# ---------------------------------------------------------------------------
# Short linear motif (SLiM) evidence from ELM (Kumar et al. 2024, Nucleic
# Acids Res 52:D442; http://elm.eu.org), downloaded by 11_download_elm.jl.
#
# Per residue, Label_SLiM is 1 when the residue lies in
#   - an annotated ELM instance of the protein lying entirely in DisProt
#     disordered residues (labels 1/2); instances on ordered residues are
#     not binding sites of interest here, or
#   - a match of an ELM class regex lying entirely in a DisProt binding region
#     (residues labelled 2),
# and NaN elsewhere (no evidence; never a 0 vote). Only the ELM classes of the
# selected types are used (default: LIG ligand, DOC docking, DEG degron).
# There is no disorder information in this source.
# ---------------------------------------------------------------------------

"""
    SLIM_DEFAULT_ELM_TYPES

ELM class types (`ELMIdentifier` prefix) used by default: binding motifs
recognised by a partner (`LIG` ligand, `DOC` docking, `DEG` degron). `MOD`
(modification sites), `CLV` (cleavage) and `TRG` (targeting) are not used.
"""
const SLIM_DEFAULT_ELM_TYPES = ("LIG", "DOC", "DEG")

"""
    _elm_type(elm_identifier)::String

ELM class type of an identifier (`"LIG_SH3_3"` -> `"LIG"`).
"""
_elm_type(elm_identifier::AbstractString)::String = uppercase(first(split(
    elm_identifier, "_")))

"""
    _elm_class_regexes(elm_index; elm_types, max_probability, logger)

Compiled regexes of the ELM classes of `elm_types` whose chance-match
`Probability` (when the column exists) is at most `max_probability`.
Returns a vector of `identifier => Regex`; classes whose regex does not compile
are skipped (logged).
"""
function _elm_class_regexes(elm_index::AbstractDataFrame;
        elm_types = SLIM_DEFAULT_ELM_TYPES, max_probability::Real = 1.0,
        logger::Function = println)
    types = Set(uppercase.(string.(elm_types)))
    has_probability = "Probability" in names(elm_index)
    classes = Pair{String, Regex}[]
    for (identifier, pattern) in sort!(collect(_elm_identifier_to_regex(elm_index)))
        _elm_type(identifier) in types || continue
        if has_probability
            row = findfirst(
                ==(identifier), _strip_wrapping_quotes.(string.(elm_index.ELMIdentifier)))
            probability = tryparse(Float64, string(elm_index.Probability[row]))
            probability !== nothing && probability > max_probability && continue
        end
        try
            push!(classes, identifier => Regex(pattern))
        catch err
            logger("Skipping ELM class $(identifier): invalid regex $(repr(pattern)) ($(sprint(showerror, err)))")
        end
    end
    classes
end

"""
    _regex_matches_in_region(sequence, classes, in_region)

Matches (overlapping) of each class regex on the full `sequence` (so `^` / `\$`
keep their meaning of protein termini) lying entirely on residues where
`in_region` is `true`. Returns `(identifier, start, stop)` tuples.
"""
function _regex_matches_in_region(sequence::AbstractString, classes,
        in_region::AbstractVector{Bool})
    matches = Tuple{String, Int, Int}[]
    any(in_region) || return matches
    text = uppercase(String(sequence))
    for (identifier, regex) in classes
        for matched in eachmatch(regex, text; overlap = true)
            isempty(matched.match) && continue
            start = matched.offset
            stop = start + length(matched.match) - 1
            stop <= length(in_region) || continue
            all(view(in_region, start:stop)) && push!(matches, (identifier, start, stop))
        end
    end
    unique!(matches)
end

"""
    build_slim_labels_table(dataset, elm_instances, elm_index;
        elm_types = SLIM_DEFAULT_ELM_TYPES, max_probability = 1.0,
        excluded_instance_logic = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
        logger = println)

SLiM evidence per dataset protein from ELM: annotated instances of the protein
(`InstanceLogic` not excluded) lying entirely in DisProt disordered residues
(labels `1`/`2`) and matches of the ELM class regexes entirely
inside DisProt binding regions (labels `2`), restricted to the ELM class types
`elm_types`.

Returns one row per dataset protein with `n_instances`, `instances`,
`n_regex_matches`, `regex_matches` (`ELMIdentifier:start-end`, `;`-separated),
`slim_positions` (ranges) and `Label_SLiM` (`1` on motif residues, `NaN`
elsewhere; one value per UniProt residue).
"""
function build_slim_labels_table(dataset::AbstractDataFrame,
        elm_instances::AbstractDataFrame, elm_index::AbstractDataFrame;
        elm_types = SLIM_DEFAULT_ELM_TYPES, max_probability::Real = 1.0,
        excluded_instance_logic::AbstractVector = DEFAULT_ELM_EXCLUDED_INSTANCE_LOGIC,
        logger::Function = println)::DataFrame
    _require_columns(
        dataset, [:protein_id, :protein_disprot_id, :sequence, :residue_labels],
        "dataset")
    types = Set(uppercase.(string.(elm_types)))
    classes = _elm_class_regexes(elm_index; elm_types = types,
        max_probability = max_probability, logger = logger)

    instances = _expand_elm_instances_by_accessions(_filter_elm_instances_by_logic(
        elm_instances; excluded_instance_logic = excluded_instance_logic))
    filter!(row -> _elm_type(row.ELMIdentifier) in types, instances)
    instances_by_protein = Dict{String, Vector{Tuple{String, Int, Int}}}()
    for row in eachrow(instances)
        push!(get!(instances_by_protein, row.protein_id, Tuple{String, Int, Int}[]),
            (row.ELMIdentifier, row.start, row[:end]))
    end

    format_hits(hits) = isempty(hits) ? missing :
                        join(["$(id):$(start)-$(stop)" for (id, start, stop) in hits], ";")
    table = DataFrame(protein_disprot_id = String[], protein_id = String[],
        n_instances = Int[], instances = Union{Missing, String}[],
        n_regex_matches = Int[], regex_matches = Union{Missing, String}[],
        slim_positions = Union{Missing, String}[], Label_SLiM = String[])
    for row in eachrow(dataset)
        (ismissing(row.sequence) || ismissing(row.residue_labels)) && continue
        protein_id = String(row.protein_id)
        sequence = uppercase(string(row.sequence))
        residue_labels = _decode_residue_labels(string(row.residue_labels))
        n_positions = length(sequence)
        length(residue_labels) == n_positions || throw(DimensionMismatch(
            "$(protein_id): $(length(residue_labels)) residue labels for $(n_positions) residues."))

        # Instances count only when entirely in DisProt disordered residues.
        in_disorder = [label in (_RESIDUE_LABEL_DISORDERED, _RESIDUE_LABEL_BINDING)
                       for label in residue_labels]
        protein_instances = sort!([hit
                                   for hit in get(instances_by_protein, protein_id, [])
                                   if 1 <= hit[2] <= hit[3] <= n_positions &&
            all(view(in_disorder, hit[2]:hit[3]))])
        in_binding = [label == _RESIDUE_LABEL_BINDING for label in residue_labels]
        regex_hits = sort!(_regex_matches_in_region(sequence, classes, in_binding);
            by = hit -> (hit[2], hit[3], hit[1]))

        label = fill(NaN, n_positions)
        for (_, start, stop) in vcat(protein_instances, regex_hits)
            label[start:stop] .= 1.0
        end
        push!(table,
            (String(row.protein_disprot_id), protein_id,
                length(protein_instances), format_hits(protein_instances),
                length(regex_hits), format_hits(regex_hits),
                _format_position_ranges(findall(==(1.0), label)), _format_label_vector(label)))
    end
    sort!(table, [:protein_disprot_id, :protein_id])
end

"""
    _SLIM_LABELS_COLUMN_DEFINITIONS

Data dictionary of the table returned by `build_slim_labels_table`.
"""
const _SLIM_LABELS_COLUMN_DEFINITIONS = [
    "protein_disprot_id" => "DisProt identifier of the dataset protein (key shared with the other derived CSVs)",
    "protein_id" => "UniProt accession of the dataset protein",
    "n_instances" => "annotated ELM instances of the protein lying entirely in DisProt disordered residues (labels `1`/`2`; selected ELM types, excluded `InstanceLogic` removed)",
    "instances" => "annotated instances (`ELMIdentifier:start-end`, `;`-separated)",
    "n_regex_matches" => "matches of ELM class regexes lying entirely in a DisProt binding region (labels `2`)",
    "regex_matches" => "regex matches (`ELMIdentifier:start-end`, `;`-separated)",
    "slim_positions" => "residues in an instance or a regex match (ranges)",
    "Label_SLiM" => "`1` on residues of an annotated ELM instance in DisProt disordered residues or of an ELM regex match in a DisProt binding region, `NaN` elsewhere (no evidence; no disorder information)"
]
