using TestItems

@testitem "PDB quality filter for aggregation" begin
    using Test
    using IDRBindingDataset: _passes_pdb_quality

    @test _passes_pdb_quality("resolved", 2.0)
    @test !_passes_pdb_quality("resolved", 6.0)
    @test !_passes_pdb_quality("resolved", 5.0)
    @test _passes_pdb_quality("resolved", 6.0; max_resolution = 8.0)
    @test _passes_pdb_quality("no_resolution", missing)
    @test !_passes_pdb_quality("no_resolution", missing; keep_no_resolution = false)
    @test _passes_pdb_quality("absent", missing)
    @test !_passes_pdb_quality("absent", missing; keep_no_info = false)
end

@testitem "Final rates and per-source labels" begin
    using Test
    using IDRBindingDataset: _pdb_vote, _final_binding_rates, _source_labels,
                             _format_label_vector, _parse_label_vector

    @test _pdb_vote(0.25, :fraction) == 0.25
    @test _pdb_vote(0.25, :any) == 1.0
    @test _pdb_vote(0.25, :majority) == 0.0
    @test_throws ArgumentError _pdb_vote(0.5, :unknown)

    disprot = [0.0, 1.0, 2.0, 1.0]
    uniprot = [0.0, 2.0, 0.0, 0.0]
    pdb_n = [1, 2, 0, 1]
    pdb_p = [0.0, 0.5, NaN, 0.0]
    rates = _final_binding_rates(disprot, uniprot, true, pdb_n, pdb_p)
    @test rates.n_methods == [3, 3, 2, 3]
    @test rates.x_methods == [0.0, 1.5, 1.0, 0.0]
    @test rates.rate ≈ [0.0, 0.5, 0.5, 0.0]
    @test rates.final_label ≈ [0.0, 1.5, 1.5, 1.0]

    in_idr = [false, true, true, true]
    pdb = _source_labels(disprot, in_idr, pdb_n .> 0, [0.0, 0.5, 0.0, 0.0])
    @test isequal(pdb.label, [0.0, 1.5, NaN, 1.0])
    @test isequal(pdb.disorder, [0.0, 1.0, NaN, 1.0])
    @test isequal(pdb.binding, [NaN, 0.5, NaN, 0.0])
    @test all(isnan, _source_labels(disprot, in_idr, fill(false, 4), zeros(4)).label)

    @test _format_label_vector([0.0, 1.23456, NaN]) == "[0.0,1.235,NaN]"
    @test isequal(_parse_label_vector("[0.0,1.235,NaN]"), [0.0, 1.235, NaN])
    @test _parse_label_vector("[]") == Float64[]
end

@testitem "Aggregated labels combine DisProt, UniProt and PDB contacts" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_aggregated_labels_table, _parse_label_vector

    # P1: residues 1-2 ordered, 3-6 disordered, 5 is a DisProt binding residue.
    dataset = DataFrame(protein_id = ["P1", "P2"], protein_disprot_id = ["DP1", "DP2"],
        residue_labels = ["001121", "012"])
    uniprot = DataFrame(protein_id = ["P1", "P2"],
        residues_binding_in_uniprot = ["3", missing])
    contacts = DataFrame(
        protein_disprot_id = ["DP1", "DP1", "DP1", "DP1", "DP2"],
        protein_id = ["P1", "P1", "P1", "P1", "P2"],
        PDB = ["1aaa", "2bbb", "3ccc", "4ddd", "5eee"],
        status = ["ok", "ok", "no_heteromeric_partner", "ok", "no_mmcif"],
        quality_category = ["resolved", "no_resolution", "resolved", "resolved", "absent"],
        resolution = [2.0, missing, 2.0, 8.0, missing],
        # 1aaa: contacts on 1 (ordered, ignored) and 4; 2bbb: no contact.
        resolved_positions = ["1-4", "3-6", "1-6", "1-6", missing],
        contact_positions = ["1;4", missing, missing, "5", missing])

    result = build_aggregated_labels_table(dataset, uniprot, contacts)
    table = result.table

    @test table.protein_disprot_id == ["DP1"]
    @test only(table.n_pdb) == 2
    @test only(table.pdb_ids) == "1AAA;2BBB"
    @test only(table.has_uniprot)
    @test _parse_label_vector(only(table.pdb_N)) == [1, 1, 2, 2, 1, 1]
    @test isequal(_parse_label_vector(only(table.pdb_proba)),
        [0.0, 0.0, 0.0, 0.5, 0.0, 0.0])
    @test isequal(_parse_label_vector(only(table.Label_PDB)),
        [0.0, 0.0, 1.0, 1.5, 1.0, 1.0])
    @test isequal(_parse_label_vector(only(table.Label_Uniprot)),
        [0.0, 0.0, 2.0, 1.0, 1.0, 1.0])
    @test isequal(_parse_label_vector(only(table.Label_Uniprot_binding)),
        [NaN, NaN, 1.0, 0.0, 0.0, 0.0])
    # Residue 5: DisProt 2, UniProt 0, PDB 0 -> rate 1/3.
    @test _parse_label_vector(only(table.rate))[5] ≈ 0.333

    @test result.pdb_log.reason ==
          [
        "used", "used", "status:no_heteromeric_partner", "low_quality", "status:no_mmcif"]

    with_homomers = build_aggregated_labels_table(dataset, uniprot, contacts;
        keep_pdb_without_partner = true, max_resolution = 10.0).table
    @test only(with_homomers.n_pdb) == 4
    @test_throws ArgumentError build_aggregated_labels_table(dataset, uniprot, contacts;
        vote = :unknown)
end
