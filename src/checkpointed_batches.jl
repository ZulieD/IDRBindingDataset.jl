using CSV
using DataFrames
using Dates

# ---------------------------------------------------------------------------
# Checkpointed, multi-threaded batch processing.
#
# Long-running per-item computations (one PDB entry, one protein, ...) are run
# in batches. Items inside a batch run in parallel on `Threads.nthreads()`
# threads. When a checkpoint directory is given, the results of every finished
# batch are appended to one CSV per output table, the identifiers of the
# finished items are appended to `done_items.txt`, and failures are written to
# `errors.log`. Restarting with the same checkpoint directory skips the items
# that are already listed in `done_items.txt`; failed items are retried.
# ---------------------------------------------------------------------------
const CHECKPOINT_DONE_FILENAME = "done_items.txt"
const CHECKPOINT_ERRORS_FILENAME = "errors.log"
const DEFAULT_CHECKPOINT_BATCH_SIZE = 15

"""
    _checkpoint_table_path(checkpoint_dir, table_name)

Return the CSV path used to accumulate rows of table `table_name` inside
`checkpoint_dir`.
"""
function _checkpoint_table_path(
        checkpoint_dir::AbstractString,
        table_name::Symbol
)::String
    joinpath(checkpoint_dir, "$(table_name).csv")
end

"""
    _read_checkpoint_done_items(checkpoint_dir)

Return the set of item identifiers already recorded as finished in
`checkpoint_dir`. Returns an empty set when no checkpoint exists yet.
"""
function _read_checkpoint_done_items(checkpoint_dir::AbstractString)::Set{String}
    done_path = joinpath(checkpoint_dir, CHECKPOINT_DONE_FILENAME)
    done_items = Set{String}()
    isfile(done_path) || return done_items
    for line in eachline(done_path)
        item = strip(line)
        isempty(item) || push!(done_items, String(item))
    end
    done_items
end

"""
    _append_checkpoint_error(checkpoint_dir, item, message; logger = println)

Append one failure entry to `errors.log`. This function never throws: when the
log cannot be written, the message is sent to `logger` instead so that a
logging problem cannot stop a long computation.
"""
function _append_checkpoint_error(
        checkpoint_dir::Union{Nothing, AbstractString},
        item::AbstractString,
        message::AbstractString;
        logger::Function = println
)::Nothing
    header = "item = $(item) | $(Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SS"))Z"
    if checkpoint_dir === nothing
        logger("ERROR $(header)\n$(message)")
        return nothing
    end
    try
        mkpath(checkpoint_dir)
        open(joinpath(checkpoint_dir, CHECKPOINT_ERRORS_FILENAME), "a") do io
            println(io, "="^80)
            println(io, header)
            println(io, message)
        end
    catch err
        logger("Could not write the checkpoint error log ($(sprint(showerror, err))).\nERROR $(header)\n$(message)")
    end
    nothing
end

"""
    _append_checkpoint_tables(checkpoint_dir, tables)

Append every non-empty `DataFrame` in the `tables` named tuple to its checkpoint
CSV, writing the header only when the CSV does not exist yet.
"""
function _append_checkpoint_tables(
        checkpoint_dir::AbstractString,
        tables::NamedTuple
)::Nothing
    for (table_name, table) in pairs(tables)
        nrow(table) == 0 && continue
        path = _checkpoint_table_path(checkpoint_dir, table_name)
        has_header = isfile(path) && filesize(path) > 0
        CSV.write(path, table; append = has_header, writeheader = !has_header)
    end
    nothing
end

"""
    _read_checkpoint_table(checkpoint_dir, table_name, template)

Read the accumulated checkpoint CSV of `table_name`, using the column types of
the empty `template` table. Returns a copy of `template` when nothing was
written yet.
"""
function _read_checkpoint_table(
        checkpoint_dir::AbstractString,
        table_name::Symbol,
        template::DataFrame
)::DataFrame
    path = _checkpoint_table_path(checkpoint_dir, table_name)
    (isfile(path) && filesize(path) > 0) || return copy(template)
    column_types = Dict(Symbol(name) => eltype(template[!, name])
    for name in names(template))
    table = CSV.read(path, DataFrame; types = column_types, stringtype = String)
    vcat(copy(template), select(table, names(template)))
end

"""
    _run_one_checkpoint_item(process, item, templates)

Run `process(item)` and check that it returns one `DataFrame` per table of
`templates`. Returns `(tables, error_message)` where exactly one is `nothing`.
"""
function _run_one_checkpoint_item(
        process::Function,
        item::AbstractString,
        templates::NamedTuple
)
    try
        tables = process(item)
        tables isa NamedTuple ||
            error("process(item) must return a NamedTuple of DataFrames, got $(typeof(tables)).")
        keys(tables) == keys(templates) ||
            error("process(item) returned tables $(keys(tables)), expected $(keys(templates)).")
        return (tables, nothing)
    catch err
        message = sprint(showerror, err) * "\n" *
                  sprint(Base.show_backtrace, catch_backtrace())
        return (nothing, message)
    end
end

"""
    run_checkpointed_batches(process, items, templates;
        checkpoint_dir = nothing,
        batch_size = DEFAULT_CHECKPOINT_BATCH_SIZE,
        logger = println
    )

Run `process(item)` for every identifier in `items`, `batch_size` items at a
time, using all available threads inside each batch.

`process(item)` must return a `NamedTuple` of `DataFrame`s with the same names
and columns as the empty `templates` tables. A failing item is logged and
skipped; it does not stop the other items.

When `checkpoint_dir` is a path, results are saved after every batch and a new
run with the same directory resumes after the last saved batch (failed items
are retried). When it is `nothing`, everything stays in memory.

Returns `(tables = NamedTuple of DataFrames, failed_items = Vector{String})`.
Row order follows processing order; callers should sort the tables when a
deterministic order is needed.
"""
function run_checkpointed_batches(
        process::Function,
        items::AbstractVector{<:AbstractString},
        templates::NamedTuple;
        checkpoint_dir::Union{Nothing, AbstractString} = nothing,
        batch_size::Integer = DEFAULT_CHECKPOINT_BATCH_SIZE,
        logger::Function = println
)
    batch_size > 0 || throw(ArgumentError("`batch_size` must be greater than zero."))
    all(t -> t isa DataFrame, values(templates)) ||
        throw(ArgumentError("`templates` must be a NamedTuple of DataFrames."))

    unique_items = unique(String.(items))
    done_items = Set{String}()
    if checkpoint_dir !== nothing
        mkpath(checkpoint_dir)
        done_items = _read_checkpoint_done_items(checkpoint_dir)
        isempty(done_items) ||
            logger("Resuming from $(checkpoint_dir): $(length(done_items)) item(s) already done.")
    end
    pending = [item for item in unique_items if !(item in done_items)]
    logger("Items: $(length(unique_items)) total, $(length(pending)) to process with $(Threads.nthreads()) thread(s).")

    in_memory = Dict(name => DataFrame[] for name in keys(templates))
    failed_items = String[]

    for batch_start in 1:batch_size:length(pending)
        batch = pending[batch_start:min(batch_start + batch_size - 1, length(pending))]
        results = Vector{Any}(undef, length(batch))

        Threads.@threads for index in eachindex(batch)
            results[index] = _run_one_checkpoint_item(process, batch[index], templates)
        end

        # Writing is sequential and happens in item order.
        for (item, (tables, message)) in zip(batch, results)
            if message !== nothing
                push!(failed_items, item)
                _append_checkpoint_error(checkpoint_dir, item, message; logger = logger)
                continue
            end
            if checkpoint_dir === nothing
                for (name, table) in pairs(tables)
                    push!(in_memory[name], table)
                end
                continue
            end
            try
                _append_checkpoint_tables(checkpoint_dir, tables)
                open(joinpath(checkpoint_dir, CHECKPOINT_DONE_FILENAME), "a") do io
                    println(io, item)
                end
            catch err
                push!(failed_items, item)
                _append_checkpoint_error(checkpoint_dir, item,
                    "Could not save results: $(sprint(showerror, err))"; logger = logger)
            end
        end

        logger("Batch $(batch_start)-$(batch_start + length(batch) - 1)/$(length(pending)) done.")
    end

    isempty(failed_items) ||
        logger("$(length(failed_items)) item(s) failed; they will be retried on the next run with the same checkpoint directory.")

    tables = map(keys(templates)) do name
        template = templates[name]
        if checkpoint_dir === nothing
            vcat(copy(template), in_memory[name]...; cols = :setequal)
        else
            _read_checkpoint_table(checkpoint_dir, name, template)
        end
    end

    (tables = NamedTuple{keys(templates)}(tables), failed_items = failed_items)
end
