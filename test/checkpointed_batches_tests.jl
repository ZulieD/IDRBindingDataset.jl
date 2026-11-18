using TestItems

@testitem "Checkpointed batches run in memory and keep failures apart" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: run_checkpointed_batches

    template = (values = DataFrame(item = String[], value = Union{Missing, Int}[]),)
    process(item) = item == "bad" ? error("boom") :
                    (values = DataFrame(item = [item], value = [length(item)]),)

    messages = String[]
    result = run_checkpointed_batches(process, ["a", "bb", "bad", "a"], template;
        batch_size = 2, logger = message -> push!(messages, message))

    @test sort(result.tables.values.item) == ["a", "bb"]
    @test result.failed_items == ["bad"]
    @test eltype(result.tables.values.value) == Union{Missing, Int}
    @test any(occursin("boom", message) for message in messages)
    @test_throws ArgumentError run_checkpointed_batches(process, ["a"], template;
        batch_size = 0)
end

@testitem "Checkpointed batches resume from the checkpoint directory" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: run_checkpointed_batches, CHECKPOINT_DONE_FILENAME,
                             CHECKPOINT_ERRORS_FILENAME

    checkpoint_dir = mktempdir()
    template = (values = DataFrame(item = String[], value = Int[]),)
    calls = String[]
    fail = Ref(true)
    function process(item)
        push!(calls, item)
        (item == "c" && fail[]) && error("transient failure")
        (values = DataFrame(item = [item], value = [Int(item[1])]),)
    end

    first_run = run_checkpointed_batches(process, ["a", "b", "c"], template;
        checkpoint_dir = checkpoint_dir, batch_size = 1, logger = _ -> nothing)
    @test first_run.failed_items == ["c"]
    @test sort(first_run.tables.values.item) == ["a", "b"]
    @test isfile(joinpath(checkpoint_dir, CHECKPOINT_ERRORS_FILENAME))
    @test readlines(joinpath(checkpoint_dir, CHECKPOINT_DONE_FILENAME)) == ["a", "b"]

    fail[] = false
    empty!(calls)
    second_run = run_checkpointed_batches(process, ["a", "b", "c"], template;
        checkpoint_dir = checkpoint_dir, batch_size = 1, logger = _ -> nothing)
    @test calls == ["c"]
    @test isempty(second_run.failed_items)
    @test sort(second_run.tables.values.item) == ["a", "b", "c"]
    @test eltype(second_run.tables.values.value) == Int
end
