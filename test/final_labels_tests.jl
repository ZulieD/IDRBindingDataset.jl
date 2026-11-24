using TestItems

@testitem "Final residue labels average the covering sources" begin
    using Test
    using IDRBindingDataset: _final_residue_labels

    # Residue 1 ordered, 2-5 disordered, 5 is a DisProt binding residue.
    residue_labels = [0, 1, 1, 1, 2]
    uniprot = [0.0, 1.0, 0.0, 0.0, 0.0]
    pdb = [NaN, 0.5, NaN, 0.0, 1.0]
    nmr = [1.0, NaN, NaN, NaN, 1.0]

    # Without NMR disorder, disorder is DisProt only.
    # 2: DisProt 0, UniProt 1, PDB 0.5 -> 1.5 / 3; 3: DisProt 0, UniProt 0 -> 0;
    # 4: DisProt 0, UniProt 0, PDB 0 -> 0; 5: all four sources -> 3 / 4.
    labels = _final_residue_labels(residue_labels, [uniprot, pdb, nmr])
    @test isequal(labels.binding, [NaN, 0.5, 0.0, 0.0, 0.75])
    @test labels.label == [0.0, 1.5, 1.0, 1.0, 1.75]
    @test labels.disorder == [0.0, 1.0, 1.0, 1.0, 1.0]
    @test labels.n_sources == [0, 3, 2, 3, 4]

    # DisProt alone: binding residues get 2, disordered ones 1.
    alone = _final_residue_labels(residue_labels, [nothing, nothing, nothing])
    @test alone.label == [0.0, 1.0, 1.0, 1.0, 2.0]
    @test isequal(alone.binding, [NaN, 0.0, 0.0, 0.0, 1.0])

    # NMR disorder: residue 1 ordered in DisProt but disordered for NMR (d = 0.5,
    # binding read from UniProt and NMR); residue 3 disordered in DisProt but
    # ordered for NMR (d = 0.5); residue 4 not evaluable by NMR (d = 1).
    nmr_disorder = [1.0, 1.0, 0.0, NaN, 1.0]
    with_nmr = _final_residue_labels(residue_labels, [uniprot, pdb, nmr];
        nmr_disorder = nmr_disorder)
    @test with_nmr.disorder == [0.5, 1.0, 0.5, 1.0, 1.0]
    @test with_nmr.binding ≈ [1 / 3, 0.5, 0.0, 0.0, 0.75]
    @test with_nmr.label ≈ [0.5 + 1 / 3, 1.5, 0.5, 1.0, 1.75]

    # Ordered for every disorder source: label 0, no binding value.
    ordered = _final_residue_labels([0], [[1.0]]; nmr_disorder = [0.0])
    @test ordered.label == [0.0]
    @test isnan(only(ordered.binding))

    @test_throws DimensionMismatch _final_residue_labels(residue_labels, [zeros(3)])
    @test_throws DimensionMismatch _final_residue_labels(residue_labels, [nothing];
        nmr_disorder = zeros(3))
end

@testitem "Final labels table combines DisProt, UniProt, PDB and chemical shifts" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_final_labels_table, _parse_label_vector

    dataset = DataFrame(protein_id = ["P2", "P1"], protein_disprot_id = ["DP2", "DP1"],
        residue_labels = ["012", "0112"])
    uniprot = DataFrame(protein_id = ["P1", "P2"],
        residues_binding_in_uniprot = ["2", missing])
    aggregated = DataFrame(protein_id = ["P1"],
        Label_PDB_binding = ["[NaN,0.5,0.0,NaN]"])
    nmr = DataFrame(protein_id = ["P1", "P2"],
        Label_chemical_shift_disorder = ["[NaN,0.0,1.0,1.0]", "[NaN,NaN,NaN]"],
        Label_chemical_shift_binding = ["[NaN,NaN,1.0,1.0]", "[NaN,NaN,NaN]"])

    # SLiM on P1 residue 3 only (no vote elsewhere).
    slim = DataFrame(protein_id = ["P1"], Label_SLiM = ["[NaN,NaN,1.0,NaN]"])

    result = build_final_labels_table(dataset, uniprot, aggregated, nmr, slim)
    table = result.table

    @test names(table) == ["protein_disprot_id", "protein_id", "Label_aggregate_all",
        "Label_aggregate_all_disorder", "Label_aggregate_all_binding"]
    @test table.protein_disprot_id == ["DP1", "DP2"]

    # P1 residue 2: disordered in DisProt, ordered for NMR (d = 0.5); binding
    # DisProt 0, UniProt 1, PDB 0.5 -> 0.5. Residue 3: d = 1; DisProt 0,
    # UniProt 0, PDB 0, NMR 1, SLiM 1 -> 2/5. Residue 4: d = 1; DisProt 1,
    # UniProt 0, NMR 1 -> 2/3.
    @test _parse_label_vector(table.Label_aggregate_all_disorder[1]) ==
          [0.0, 0.5, 1.0, 1.0]
    @test isequal(_parse_label_vector(table.Label_aggregate_all_binding[1]),
        [NaN, 0.5, 0.4, 0.667])
    @test _parse_label_vector(table.Label_aggregate_all[1]) == [0.0, 1.0, 1.4, 1.667]

    # P2: DisProt only (no UniProt site, no PDB, NMR without covered residue).
    @test _parse_label_vector(table.Label_aggregate_all[2]) == [0.0, 1.0, 2.0]

    sources = result.sources
    @test sources.sources == ["disprot;uniprot;pdb;chemical_shift;slim", "disprot"]
    @test sources.n_disordered == [3, 2]
    @test sources.n_partial_disorder == [1, 0]
    @test sources.n_nmr_only_disorder == [0, 0]
    @test sources.n_3_sources == [2, 0]
    @test sources.n_4_sources == [0, 0]
    @test sources.n_5_sources == [1, 0]
    @test sources.n_1_source == [0, 2]
end
