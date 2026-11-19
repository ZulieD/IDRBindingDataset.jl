using TestItems

@testmodule BMRBFixtures begin
    # Synthetic BMRB API v2 entries (NMR-STAR JSON), built like the real ones.
    const PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
    const CONSTRUCT = "GS" * PROTEIN[5:60]

    saveframe(category, tags; loops = Any[]) = Dict("category" => category,
        "name" => category, "tags" => [[k, v] for (k, v) in tags], "loops" => loops)
    loop(category, tags, rows) = Dict("category" => "_" * category, "tags" => tags,
        "data" => rows)

    function entity(id, name, type, sequence; uniprot = nothing)
        loops = uniprot === nothing ? Any[] :
                Any[loop("Entity_db_link", ["Database_code", "Accession_code"],
            [["UNP", uniprot]])]
        saveframe("entity",
            [("ID", id), ("Name", name), ("Type", type),
                ("Polymer_type", type == "polymer" ? "polypeptide(L)" : "."),
                ("Polymer_seq_one_letter_code", sequence)]; loops = loops)
    end

    "Entry with the protein (entity 1), extra `entities` and the `assembly` entity IDs."
    function entry(bmrb_id; entities = Any[], assembly = ["1"], with_shifts = true,
            celsius = false)
        shifts = Any[]
        for (index, residue) in enumerate(CONSTRUCT)
            index <= 2 && continue
            for atom in (residue == 'P' ? ("N", "CA") : ("H", "N", "CA", "CB"))
                push!(shifts, [
                    string(length(shifts) + 1), "1", string(index), "XXX", atom, "1.0"])
            end
        end
        push!(shifts, ["999", "2", "1", "TRP", "H", "8.1"])   # partner atom: not counted
        frames = Any[
            saveframe("entry_information", [
                ("ID", bmrb_id), ("Title", "Test   entry\n$bmrb_id")]),
            entity("1", "Test protein", "polymer", CONSTRUCT),
            entities...,
            saveframe("assembly", [("ID", "1")];
                loops = Any[loop("Entity_assembly", ["ID", "Entity_ID"],
                    [[string(i), e] for (i, e) in enumerate(assembly)])]),
            saveframe("sample_conditions", [("ID", "1")];
                loops = Any[loop("Sample_condition_variable", ["Type", "Val", "Val_units"],
                    [["pH", "6.5", "pH"],
                        ["temperature", celsius ? "25" : "298", celsius ? "C" : "K"],
                        ["ionic strength", "0.15", "M"]])]),
            saveframe("sample", [("ID", "1"), ("Solvent_system", "90% H2O/10% D2O")];
                loops = Any[loop("Sample_component",
                    ["Mol_common_name", "Entity_ID", "Concentration_val",
                        "Concentration_val_units"],
                    [["Test protein", "1", "0.5", "mM"],
                        ["sodium phosphate", ".", "20", "mM"],
                        ["NaCl", ".", "50", "mM"], ["D2O", ".", "10", "%"], [
                            "DSS", ".", "0.1", "mM"]])]),
            saveframe("NMR_spectrometer",
                [("ID", "1"), ("Manufacturer", "Bruker"), ("Model", "Avance III"),
                    ("Field_strength", "600")]),
            saveframe("experiment_list", [("ID", "1")];
                loops = Any[loop(
                    "Experiment", ["ID", "Name", "Sample_ID", "NMR_spectrometer_ID"],
                    [["1", "2D 1H-15N HSQC", "1", "1"]])]),
            saveframe("spectral_peak_list", [("ID", "1"), ("Experiment_ID", "1")];
                loops = Any[loop("Spectral_dim",
                    ["Atom_type", "Sweep_width",
                        "Sweep_width_units", "Spectrometer_frequency"],
                    [["H", "14", "ppm", "600"], ["N", "2100", "Hz", "60"]])])
        ]
        if with_shifts
            push!(frames,
                saveframe("assigned_chemical_shifts",
                    [("ID", "1"), ("Sample_condition_list_ID", "1")];
                    loops = Any[
                        loop("Chem_shift_experiment", ["Experiment_ID", "Sample_ID"],
                            [["1", "1"]]),
                        loop("Atom_chem_shift",
                            ["ID", "Entity_ID", "Seq_ID", "Comp_ID", "Atom_ID", "Val"], shifts)]))
        end
        Dict(bmrb_id => Dict("bmrb_id" => bmrb_id, "saveframes" => frames))
    end
end

@testitem "BMRB UniProt mapping parser accepts dictionaries and records" begin
    using Test
    using IDRBindingDataset: parse_bmrb_uniprot_mapping

    @test parse_bmrb_uniprot_mapping(Dict("p12345" => ["1", "2"], "Q9-2" => "3")) ==
          Dict("P12345" => ["1", "2"], "Q9" => ["3"])
    records = [Dict("uniprot_id" => "P1", "bmrb_ids" => ["10", "11", "10"]),
        Dict("UniProt" => "P2-1", "BMRB_ID" => 12), ["P3", "13"], Dict("other" => 1)]
    @test parse_bmrb_uniprot_mapping(records) ==
          Dict("P1" => ["10", "11"], "P2" => ["12"], "P3" => ["13"])
end

@testitem "Local alignment maps an expression construct onto UniProt" begin
    using Test
    using IDRBindingDataset: _local_alignment, _clean_bmrb_sequence

    protein = "MKTAYIAKQRQISFVKSHFSRQ"
    alignment = _local_alignment("GS" * protein[5:15], protein)
    @test alignment.identity == 1.0
    @test alignment.n_aligned == 11
    @test alignment.indels == 0
    @test alignment.mapping[3] == 5
    @test _clean_bmrb_sequence("GS(MSE)KT\nAY") == "GSXKTAY"
    @test _clean_bmrb_sequence(nothing) == ""
end

@testitem "BMRB entry rows describe free, bound and homo-oligomer samples" setup = [BMRBFixtures] begin
    using Test
    using IDRBindingDataset: _bmrb_entry_rows
    F = BMRBFixtures

    free = _bmrb_entry_rows("DP1", "P1", F.PROTEIN, "90001", F.entry("90001"))
    row = only(eachrow(free))
    @test row.status == "ok"
    @test row.state == "free"
    @test row.entry_title == "Test entry 90001"
    @test row.entity_identity == 1.0
    @test (row.cover_start, row.cover_end) == (5, 60)
    @test row.n_partners == 0
    @test ismissing(row.partner_names)
    @test row.pH == 6.5
    @test row.temperature_K == 298.0
    @test row.ionic_strength_mM == 150.0
    @test row.solvent == "90% H2O/10% D2O"
    @test row.buffer == "NaCl 50 mM | sodium phosphate 20 mM"
    @test row.experiments == "2D 1H-15N HSQC"
    @test row.spectrometers == "Bruker Avance III"
    @test row.field_strengths_MHz == "600.0"
    @test row.n15_sweep_width_ppm == 35.0
    @test row.n_residues_with_shifts == 56
    @test row.n_H + row.n_N + row.n_CA + row.n_CB == row.n_shifts

    peptide = F.entity("2", "Partner peptide", "polymer", "WWPEEFLKRA"; uniprot = "Q99999")
    zinc = F.entity("3", "ZINC ION", "non-polymer", ".")
    bound = _bmrb_entry_rows("DP1", "P1", F.PROTEIN, "90002",
        F.entry("90002"; entities = [peptide, zinc], assembly = ["1", "2", "3"],
            celsius = true))
    row = only(eachrow(bound))
    @test row.state == "bound"
    @test row.n_partners == 2
    @test row.partner_names == "Partner peptide | ZINC ION"
    @test row.partner_types == "polymer | non-polymer"
    @test row.partner_lengths == "10 | 0"
    @test row.partner_uniprot_ids == "Q99999 | "
    @test row.temperature_K ≈ 298.15

    homo = _bmrb_entry_rows("DP1", "P1", F.PROTEIN, "90003",
        F.entry("90003"; assembly = ["1", "1"]))
    @test only(homo.state) == "homo_oligomer"
    @test only(homo.n_target_copies) == 2

    no_shifts = _bmrb_entry_rows("DP1", "P1", F.PROTEIN, "90004",
        F.entry("90004"; with_shifts = false))
    @test only(no_shifts.status) == "no_shift_list"
    @test only(no_shifts.state) == "free"
end

@testitem "BMRB mapping table downloads entries once and reports failures" setup = [BMRBFixtures] begin
    using Test
    using HTTP
    using JSON
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_disprot_bmrb_mapping_table, fetch_bmrb_uniprot_mapping,
                             BMRB_API_URL
    F = BMRBFixtures

    requested = String[]
    function mock_get(url)
        push!(requested, url)
        occursin("mappings", url) &&
            return HTTP.Response(200, JSON.json(Dict("P1" => ["90001", "90009"])))
        occursin("90009", url) && return HTTP.Response(404, "not found")
        HTTP.Response(200, JSON.json(F.entry("90001")))
    end
    cache_dir = mktempdir()
    options = (cache_dir = cache_dir, get_request = mock_get, max_retries = 2,
        base_sleep = 0.0, sleeper = _ -> nothing)

    mapping = fetch_bmrb_uniprot_mapping(; options...)
    @test mapping == Dict("P1" => ["90001", "90009"])

    dataset = DataFrame(protein_id = ["P1", "P2"], protein_disprot_id = ["DP1", "DP2"],
        sequence = [F.PROTEIN, "MKT"])
    result = build_disprot_bmrb_mapping_table(dataset, mapping; options...,
        logger = _ -> nothing)
    @test result.table.bmrb_id == ["90001", "90009"]
    @test result.table.status == ["ok", "error"]
    @test result.failed_entries == ["P1/90009"]
    @test isfile(joinpath(cache_dir, "entry_90001.json"))

    # Second run: the mapping and the cached entry are not downloaded again.
    empty!(requested)
    fetch_bmrb_uniprot_mapping(; options...)
    build_disprot_bmrb_mapping_table(dataset, mapping; options..., logger = _ -> nothing)
    @test requested == ["$(BMRB_API_URL)/entry/90009",
        "$(BMRB_API_URL)/entry/90009"]
end
