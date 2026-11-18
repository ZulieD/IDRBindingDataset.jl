using TestItems

@testitem "Dataset PDB segments merge chains and drop unused SIFTS columns" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_dataset_pdb_segments

    dataset = DataFrame(protein_id = ["P1", "P2"], protein_disprot_id = ["DP1", "DP2"])
    sifts = DataFrame(
        PDB = ["1ABC", "1abc", "1abc", "2xyz", "3aaa", "1abc"],
        CHAIN = ["B", "A", "C", "A", "A", "D"],
        SP_PRIMARY = ["P1", "P1", "P1", "P2", "P9", "P1"],
        RES_BEG = ["1", "1", "1", "1", "1", ""],
        RES_END = ["9", "9", "9", "9", "9", ""],
        PDB_BEG = ["1", "1", "1", "1", "1", ""],
        PDB_END = ["9", "9", "9", "9", "9", ""],
        SP_BEG = ["1", "1", "20", "3", "1", "x"],
        SP_END = ["10", "10", "30", "8", "5", "10"]
    )

    segments = build_dataset_pdb_segments(dataset, sifts)

    @test names(segments) ==
          ["protein_disprot_id", "protein_id", "PDB", "CHAIN", "SP_BEG", "SP_END"]
    @test segments.protein_disprot_id == ["DP1", "DP1", "DP2"]
    @test segments.protein_id == ["P1", "P1", "P2"]
    @test segments.PDB == ["1abc", "1abc", "2xyz"]
    @test segments.CHAIN == ["A;B", "C", "A"]
    @test segments.SP_BEG == [1, 20, 3]
    @test segments.SP_END == [10, 30, 8]
end

@testitem "Missing residues are counted on UniProt-mapped residues only" begin
    using Test
    using IDRBindingDataset: SIFTSResidueRecord, _count_missing_residues,
                             _select_best_chains

    labels = [0, 1, 1, 2, 2, 0]
    record(position, is_missing; uniprot_id = "P1") = (chain = "A",
        pdb_number = string(position), uniprot_id = uniprot_id,
        uniprot_position = position, is_missing = is_missing)
    records = SIFTSResidueRecord[
        record(1, true), record(2, true), record(3, false), record(4, true),
        record(5, true), record(6, true),
        record(0, true; uniprot_id = ""),          # expression tag
        record(2, true; uniprot_id = "Q9")]        # other protein

    counts = _count_missing_residues(records, labels, "p1", 2, 5)
    @test counts == (n_missing = 3, n_missing_disorder = 1, n_missing_disorder_binding = 2)

    full = _count_missing_residues(records, labels, "P1", 1, 6)
    @test full.n_missing == 5

    better = (n_missing = 1, n_missing_disorder = 0, n_missing_disorder_binding = 1)
    worse = (n_missing = 4, n_missing_disorder = 2, n_missing_disorder_binding = 0)
    @test _select_best_chains(["A" => worse, "B" => better, "C" => better]) ==
          (["B", "C"], better)
    @test _select_best_chains(Pair{String, Any}[]) == (String[], nothing)
end

@testitem "Missing-residue table keeps segments without SIFTS XML" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_pdb_missing_residues_table

    dataset = DataFrame(protein_id = ["P1"], protein_disprot_id = ["DP1"],
        sequence = ["MKTAYI"],
        residue_labels = ["011220"])
    segments = DataFrame(protein_disprot_id = ["DP1", "DP1"],
        protein_id = ["P1", "P1"], PDB = ["9zzz", "9zzz"],
        CHAIN = ["A;B", "C"], SP_BEG = [1, 2], SP_END = [6, 4])

    result = build_pdb_missing_residues_table(dataset, segments;
        sifts_xml_dirs = [mktempdir()], logger = _ -> nothing)
    table = result.table

    @test isempty(result.failed_pdbs)
    @test nrow(table) == 2
    @test table.protein_disprot_id == ["DP1", "DP1"]
    @test table.COMPLETE == [true, false]
    @test table.CHAIN_KEPT == ["A;B", "C"]
    @test all(ismissing, table.N_MISSING)
end
