using TestItems

@testitem "ELM regex matches inside DisProt binding regions" begin
    using Test
    using IDRBindingDataset: _regex_matches_in_region, _elm_type

    @test _elm_type("LIG_SH3_3") == "LIG"
    @test _elm_type("DOC_MAPK_1") == "DOC"

    sequence = "MAPPLPPKAAPPLPPK"
    classes = ["LIG_SH3_test" => r"P..P", "LIG_Nterm" => r"^MA"]
    in_region = [i in 2:8 for i in 1:length(sequence)]
    # P..P at 3-6 and 4-7 (overlapping) are inside 2-8; 11-14 and 12-15 are
    # outside; ^MA (1-2) starts before the region.
    @test _regex_matches_in_region(sequence, classes, in_region) ==
          [("LIG_SH3_test", 3, 6), ("LIG_SH3_test", 4, 7)]
    @test isempty(_regex_matches_in_region(sequence, classes, falses(length(sequence))))
end

@testitem "SLiM labels from ELM instances and regexes" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_slim_labels_table, _parse_label_vector

    # P1: binding region 9-14 contains "RPLPPL" (LIG regex R.LP.L at 9-14).
    dataset = DataFrame(protein_id = ["P1", "P2"], protein_disprot_id = ["DP1", "DP2"],
        sequence = ["MKTAYIAKRPLPPLAAAA", "MAAA"],
        residue_labels = ["000011112222220000", "0011"])
    elm_index = DataFrame(ELMIdentifier = ["LIG_test", "MOD_test", "DOC_rare"],
        Regex = ["R.LP.L", "AA", "WWWW"], Probability = ["0.001", "0.1", "0.2"])
    elm_instances = DataFrame(
        ELMIdentifier = ["DOC_rare", "LIG_test", "MOD_test", "LIG_test"],
        Primary_Acc = ["P1", "P1", "P1", "P2"], Accessions = ["P1", "P1", "P1", "P2 Q9"],
        Start = ["2", "16", "15", "1"], End = ["4", "17", "18", "9"],
        InstanceLogic = ["true positive", "false positive", "true positive",
            "true positive"])

    table = build_slim_labels_table(dataset, elm_instances, elm_index;
        logger = _ -> nothing)
    @test names(table) == ["protein_disprot_id", "protein_id", "n_instances", "instances",
        "n_regex_matches", "regex_matches", "slim_positions", "Label_SLiM"]
    p1 = table[1, :]
    # Instances: DOC_rare 2-4 kept; LIG_test false positive and MOD_test (type
    # not selected) removed. Regex: LIG_test at 9-14; MOD "AA" not used.
    @test p1.instances == "DOC_rare:2-4"
    @test p1.regex_matches == "LIG_test:9-14"
    @test p1.slim_positions == "2-4;9-14"
    label = _parse_label_vector(p1.Label_SLiM)
    @test all(==(1.0), label[[2:4; 9:14]])
    @test all(isnan, label[[1; 5:8; 15:18]])

    # P2: instance outside the sequence (1-9 on 4 residues) is ignored.
    p2 = table[2, :]
    @test p2.n_instances == 0 && ismissing(p2.slim_positions)
    @test all(isnan, _parse_label_vector(p2.Label_SLiM))

    # max_probability skips permissive classes; elm_types selects class types.
    strict = build_slim_labels_table(dataset, elm_instances, elm_index;
        max_probability = 0.0001, logger = _ -> nothing)
    @test strict.n_regex_matches[1] == 0
    with_mod = build_slim_labels_table(dataset, elm_instances, elm_index;
        elm_types = ("MOD",), logger = _ -> nothing)
    @test with_mod.instances[1] == "MOD_test:15-18"
    @test with_mod.n_regex_matches[1] == 0          # no "AA" in the binding region
end
