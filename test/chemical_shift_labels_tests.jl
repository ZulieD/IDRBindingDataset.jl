using TestItems

@testmodule ChemicalShiftFixtures begin
    using IDRBindingDataset: _random_coil

    # Same synthetic system as the original pipeline self-test:
    # protein of 80 residues, NMR construct = "GS" + residues 5-74;
    # folded (helix) 5-34, disordered 35-74, binding site 50-56, peak 58 lost.
    const PROTEIN = "MKTAYIAKQRDEFGHLNSWVPKTAYIAKQRDEFGSEEKPSLEQNGSDKTEAPSQEGKNSDATEGSPKQEDLTGSASEPKA"
    const CONSTRUCT = "GS" * PROTEIN[5:74]
    const THREE = Dict('A' => "ALA", 'R' => "ARG", 'N' => "ASN", 'D' => "ASP",
        'C' => "CYS", 'Q' => "GLN", 'E' => "GLU", 'G' => "GLY", 'H' => "HIS",
        'I' => "ILE", 'L' => "LEU", 'K' => "LYS", 'M' => "MET", 'F' => "PHE",
        'P' => "PRO", 'S' => "SER", 'T' => "THR", 'W' => "TRP", 'Y' => "TYR",
        'V' => "VAL")

    "Deterministic pseudo-random value in [-1, 1]."
    noise(i, salt) = sin(12.9898 * i + 78.233 * salt) * 0.5 + cos(4.1414 * i * salt) * 0.5

    "True shifts (CA, CB, H, N) of protein position `p`."
    function true_shifts(p)
        residue = PROTEIN[p]
        coil = _random_coil(PROTEIN, p, residue, :ncidp)
        folded = p <= 34
        ca = coil[1] + (folded ? 3.0 : 0.15) + 0.1 * noise(p, 1)
        cb = isnan(coil[2]) ? NaN : coil[2] + (folded ? -1.0 : 0.0) + 0.1 * noise(p, 2)
        h = folded ? 8.4 + 1.4 * noise(p, 3) : 8.2 + 0.2 * noise(p, 4)
        n = folded ? 118 + 10 * noise(p, 5) : 122 + 4 * noise(p, 6)
        (ca, cb, h, n)
    end

    "Atom_chem_shift rows for one experiment."
    function shift_rows(; bound = false, binding = 50:56, salt = 0, construct = CONSTRUCT)
        rows = Vector{Vector{String}}()
        for (index, residue) in enumerate(construct)
            index <= 2 && continue                 # GS tag not assigned
            p = index - 2 + 4                      # protein position
            ca, cb, h, n = true_shifts(p)
            ca += 0.03 * noise(p, 7 + salt)
            cb += 0.03 * noise(p, 8 + salt)
            h += 0.003 * noise(p, 9 + salt)
            n += 0.02 * noise(p, 10 + salt)
            if bound
                h += 0.02                          # global referencing offset
                n += 0.1
                p in binding && (h += 0.25; n += 1.5)
                p == 58 && continue                # lost peak
            end
            for (atom, value) in (("H", h), ("N", n), ("CA", ca), ("CB", cb))
                (residue == 'P' && atom == "H") && continue
                isnan(value) && continue
                push!(rows,
                    [string(length(rows) + 1), "1", string(index), string(index),
                        THREE[residue], atom, string(round(value; digits = 3))])
            end
        end
        rows
    end

    sf(category, tags; loops = Any[]) = Dict("category" => category,
        "tags" => [[k, v] for (k, v) in tags], "loops" => loops)
    lp(category, tags, rows) = Dict("category" => "_" * category, "tags" => tags,
        "data" => rows)

    "BMRB entry JSON with the protein and optional partner entities."
    function entry(bmrb_id, rows; partners = Any[], assembly = ["1"],
            buffer = [("sodium phosphate", "20"), ("NaCl", "50")], n_width = "35",
            construct = CONSTRUCT)
        entities = Any[sf("entity",
            [("ID", "1"), ("Name", "Test protein"), ("Type", "polymer"),
                ("Polymer_type", "polypeptide(L)"),
                ("Polymer_seq_one_letter_code", construct)])]
        for (i, (name, sequence)) in enumerate(partners)
            push!(entities,
                sf("entity",
                    [("ID", string(i + 1)), ("Name", name), ("Type", "polymer"),
                        ("Polymer_type", "polypeptide(L)"),
                        ("Polymer_seq_one_letter_code", sequence)]))
        end
        components = [[name, ".", conc, "mM"] for (name, conc) in buffer]
        push!(components, ["Test protein", "1", "0.5", "mM"])
        frames = vcat(Any[sf("entry_information", [("Title", "Synthetic $bmrb_id")])],
            entities,
            Any[
                sf("assembly", [("ID", "1")];
                    loops = Any[lp("Entity_assembly", ["ID", "Entity_ID"],
                        [[string(i), e] for (i, e) in enumerate(assembly)])]),
                sf("sample_conditions", [("ID", "1")];
                    loops = Any[lp(
                        "Sample_condition_variable", ["Type", "Val", "Val_units"],
                        [["pH", "6.5", "pH"], ["temperature", "298", "K"],
                            ["ionic strength", "70", "mM"]])]),
                sf("sample", [("ID", "1"), ("Solvent_system", "90% H2O/10% D2O")];
                    loops = Any[lp("Sample_component",
                        ["Mol_common_name", "Entity_ID", "Concentration_val",
                            "Concentration_val_units"], components)]),
                sf("NMR_spectrometer",
                    [("ID", "1"), ("Manufacturer", "Bruker"), ("Model", "Avance"),
                        ("Field_strength", "600")]),
                sf("experiment_list", [("ID", "1")];
                    loops = Any[lp("Experiment",
                        ["ID", "Name", "Sample_ID", "NMR_spectrometer_ID"],
                        [["1", "2D 1H-15N HSQC", "1", "1"]])]),
                sf("spectral_peak_list", [("ID", "1"), ("Experiment_ID", "1")];
                    loops = Any[lp(
                        "Spectral_dim", ["Atom_type", "Sweep_width",
                            "Sweep_width_units"],
                        [["N", n_width, "ppm"]])]),
                sf("assigned_chemical_shifts",
                    [("ID", "1"), ("Sample_condition_list_ID", "1")];
                    loops = Any[
                        lp("Chem_shift_experiment", ["Experiment_ID", "Sample_ID"],
                            [["1", "1"]]),
                        lp("Atom_chem_shift",
                            ["ID", "Entity_ID", "Seq_ID", "Comp_index_ID", "Comp_ID",
                                "Atom_ID", "Val"], rows)])])
        Dict(bmrb_id => Dict("saveframes" => frames))
    end

    "All synthetic entries (as in the original self-test)."
    function entries()
        Dict(
            "90001" => entry("90001", shift_rows()),
            "90002" => entry("90002", shift_rows(bound = true);
                partners = [("Partner peptide", "WWPEEFLKRA")], assembly = ["1", "2"]),
            "90003" => entry("90003", shift_rows(bound = true); assembly = ["1", "1"]),
            "90004" => entry("90004", shift_rows(salt = 3)),
            "90005" => entry("90005", shift_rows(bound = true, binding = 50:53, salt = 5);
                partners = [("Other partner", "DDEEWKLPRS")], assembly = ["1", "2"]),
            "90006" => entry("90006", shift_rows(bound = true, binding = 40:45);
                partners = [("Partner in Tris", "KKRRWWLLPE")], assembly = ["1", "2"],
                buffer = [("Tris-HCl", "50"), ("NaCl", "150")]))
    end

    "Construct with protein position `p` mutated to `residue`."
    function mutated_construct(p, residue)
        characters = collect(CONSTRUCT)
        characters[p - 4 + 2] = residue
        String(characters)
    end
end

@testitem "Random coil values and run filtering" begin
    using Test
    using IDRBindingDataset: _random_coil, _keep_runs

    # ncIDP: Ala alone 52.599/19.102; Pro at i+1 lowers Cα(i) by 2.423 ppm.
    @test _random_coil("A", 1, 'A', :ncidp) == (52.599, 19.102)
    ca, _ = _random_coil("AP", 1, 'A', :ncidp)
    @test ca ≈ 52.599 - 2.423
    @test isnan(_random_coil("G", 1, 'G', :ncidp)[2])
    @test _random_coil("A", 1, 'A', :wishart) == (52.5, 19.1)
    @test _random_coil("X", 1, 'X', :ncidp) === nothing

    flags = Bool[1, 1, 1, 0, 1, 1, 1, 1, 1]
    positions = [1, 2, 3, 4, 10, 11, 13, 14, 15]
    # 10-15 spans 6 positions (gap of one at 12); 1-3 spans 3.
    @test _keep_runs(flags, positions, 5) == Bool[0, 0, 0, 0, 1, 1, 1, 1, 1]
    @test _keep_runs(flags, positions, 3) == Bool[1, 1, 1, 0, 1, 1, 1, 1, 1]
end

@testitem "Experiment comparability criteria" begin
    using Test
    using IDRBindingDataset: _compare_experiments, CHEMICAL_SHIFT_DEFAULT_SETTINGS

    base = (pH = 6.5, temperature_K = 298.0, ionic_strength_mM = 70.0,
        fields_MHz = [600.0], spectrometers = ["Bruker Avance"],
        buffer_mM = Dict{String, Union{Missing, Float64}}(
            "sodium phosphate" => 20.0, "nacl" => 50.0),
        solvents = ["90% H2O/10% D2O"], n15_sweep_width_ppm = 35.0)
    settings = CHEMICAL_SHIFT_DEFAULT_SETTINGS

    same = _compare_experiments(base, base, settings)
    @test same.accepted
    @test all(==(:same), values(same.criteria))

    other_buffer = merge(base,
        (buffer_mM = Dict{String, Union{Missing, Float64}}(
            "tris-hcl" => 20.0, "nacl" => 50.0),))
    @test !_compare_experiments(base, other_buffer, settings).accepted
    @test _compare_experiments(base, other_buffer, settings).criteria[:buffer] == :different

    # Same buffering species: concentrations, salts and additives may differ.
    same_species = merge(base,
        (buffer_mM = Dict{String, Union{Missing, Float64}}(
            "napi" => 50.0, "nacl" => 150.0, "dtt" => 1.0),))
    @test _compare_experiments(base, same_species, settings).criteria[:buffer] == :same
    no_species = merge(base,
        (buffer_mM = Dict{String, Union{Missing, Float64}}("nacl" => 50.0),))
    @test _compare_experiments(base, no_species, settings).criteria[:buffer] == :unknown

    # Temperature: up to 10 K apart.
    warm = merge(base, (temperature_K = 308.0,))
    @test _compare_experiments(base, warm, settings).criteria[:temperature] == :same

    hot = merge(base, (temperature_K = 310.0, fields_MHz = [800.0]))
    comparison = _compare_experiments(base, hot, settings)
    @test !comparison.accepted
    @test comparison.criteria[:temperature] == :different
    @test comparison.criteria[:field] == :different        # reported, not required

    unknown = merge(base, (pH = missing,))
    @test _compare_experiments(base, unknown, settings).accepted
    @test !_compare_experiments(base, unknown, merge(settings, (strict_metadata = true,))).accepted

    narrow = merge(base, (n15_sweep_width_ppm = 20.0,))
    # A 15N window narrower than 30 ppm differs, but the window is not required.
    @test _compare_experiments(base, narrow, settings).criteria[:window] == :different
    @test _compare_experiments(base, narrow, settings).accepted
    @test _compare_experiments(narrow, narrow, settings).criteria[:window] == :same
end

@testitem "Buffering species of sample components" begin
    using Test
    using IDRBindingDataset: _buffer_species, _buffer_species_mM, _format_buffer_species

    @test _buffer_species("Sodium phosphate") == "phosphate"
    @test _buffer_species("NaPi") == "phosphate"
    @test _buffer_species("K2HPO4") == "phosphate"
    @test _buffer_species("Tris-HCl") == "tris"
    @test _buffer_species("Tris-d11") == "tris"
    @test _buffer_species("Bis-Tris") == "bis-tris"
    @test _buffer_species("2-(N-morpholino)ethanesulfonic acid (MES)") == "mes"
    @test _buffer_species("HEPES-d18") == "hepes"
    @test _buffer_species("sodium acetate") == "acetate"
    @test _buffer_species("TCEP") === nothing
    @test _buffer_species("tris(2-carboxyethyl)phosphine") === nothing
    for additive in ("DTT", "EDTA", "PMSF", "NaCl", "KCl", "protease inhibitor cocktail")
        @test _buffer_species(additive) === nothing
    end

    species = _buffer_species_mM(Dict{String, Union{Missing, Float64}}(
        "na2hpo4" => 10.0, "nah2po4" => 10.0, "nacl" => 50.0, "tris" => missing))
    @test species["phosphate"] == 20.0
    @test ismissing(species["tris"])
    @test _format_buffer_species(species) == "phosphate 20.0 mM; tris ? mM"
    @test ismissing(_format_buffer_species(Dict{String, Union{Missing, Float64}}()))
end

@testitem "Mutated NMR constructs" setup = [ChemicalShiftFixtures] begin
    using Test
    using JSON
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_disprot_bmrb_mapping_table, build_chemical_shift_labels,
                             _parse_label_vector
    F = ChemicalShiftFixtures

    # 90001: wild type; 90010: mutation on 20 (outside the binding region);
    # 90011: mutation on 48, within 3 residues of the DisProt binding region 50-58.
    outside = F.mutated_construct(20, F.PROTEIN[20] == 'W' ? 'F' : 'W')
    inside = F.mutated_construct(48, F.PROTEIN[48] == 'W' ? 'F' : 'W')
    entries = Dict("90001" => F.entry("90001", F.shift_rows()),
        "90010" => F.entry("90010", F.shift_rows(construct = outside);
            construct = outside),
        "90011" => F.entry("90011", F.shift_rows(construct = inside); construct = inside))
    cache_dir = mktempdir()
    for (id, json) in entries
        write(joinpath(cache_dir, "entry_$(id).json"), JSON.json(json))
    end
    labels = join([p in 35:74 ? (p in 50:58 ? '2' : '1') : '0' for p in 1:80])
    dataset = DataFrame(protein_id = ["PSYN1"], protein_disprot_id = ["DP99999"],
        sequence = [F.PROTEIN], residue_labels = [labels])
    mapping = build_disprot_bmrb_mapping_table(dataset,
        Dict("PSYN1" => sort(collect(keys(entries)))); cache_dir = cache_dir,
        logger = _ -> nothing).table
    @test all(>=(0.95), parse.(Float64, string.(mapping.entity_identity)))

    result = build_chemical_shift_labels(dataset, mapping; cache_dir = cache_dir,
        logger = _ -> nothing)
    @test result.lists.reason == ["used", "used", "mutation_in_binding_region"]

    # The mutated residue 20 is not used: one free list instead of two.
    n_free = _parse_label_vector(only(result.labels.cs_n_free))
    @test n_free[20] == 1
    @test n_free[21] == 2
end

@testitem "Chemical shift labels on the synthetic free/bound system" setup = [ChemicalShiftFixtures] begin
    using Test
    using JSON
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_disprot_bmrb_mapping_table, build_chemical_shift_labels,
                             _parse_label_vector
    F = ChemicalShiftFixtures

    cache_dir = mktempdir()
    entries = F.entries()
    for (id, json) in entries
        write(joinpath(cache_dir, "entry_$(id).json"), JSON.json(json))
    end
    labels = join([p in 35:74 ? (p in 50:58 ? '2' : '1') : '0' for p in 1:80])
    dataset = DataFrame(protein_id = ["PSYN1"], protein_disprot_id = ["DP99999"],
        sequence = [F.PROTEIN], residue_labels = [labels])
    mapping = build_disprot_bmrb_mapping_table(dataset,
        Dict("PSYN1" => sort(collect(keys(entries)))); cache_dir = cache_dir,
        logger = _ -> nothing).table
    @test mapping.state == ["free", "bound", "homo_oligomer", "free", "bound", "bound"]

    result = build_chemical_shift_labels(dataset, mapping; cache_dir = cache_dir,
        logger = _ -> nothing)
    @test isempty(result.failed_proteins)

    lists = result.lists
    @test lists.reason == ["used", "used", "homo_oligomer", "used", "used",
        "no_comparable_free"]

    comparisons = result.comparisons
    @test nrow(comparisons) == 6                     # 3 bound lists × 2 free lists
    rejected = comparisons[comparisons.bound_bmrb_id .== "90006", :]
    @test all(rejected.buffer .== "different")
    @test !any(rejected.accepted)
    @test count(comparisons.selected) == 2

    row = only(eachrow(result.labels))
    @test row.n_free_lists == 2
    @test row.n_pairs == 2
    disorder = _parse_label_vector(row.Label_chemical_shift_disorder)
    binding = _parse_label_vector(row.Label_chemical_shift_binding)
    label = _parse_label_vector(row.Label_chemical_shift)

    # Construct ends (5 first / 5 last residues) are not analysed: data on 8-69.
    @test all(isnan, disorder[1:7]) && all(isnan, disorder[70:80])
    @test all(==(0.0), disorder[10:32])              # folded region
    @test all(==(1.0), disorder[37:67])              # disordered region
    # Prolines have no amide H, hence no CSP: they keep the disorder label only.
    observed(range) = [p for p in range if F.PROTEIN[p] != 'P']
    @test F.PROTEIN[52] == 'P' && isnan(binding[52]) && label[52] == 1.0
    # Binding 50-56 in one pair out of two (a proline gap does not break the run);
    # the isolated lost peak 58 and the 4-residue run 50-53 of 90005 are ignored.
    @test all(==(0.5), binding[observed(50:56)])
    @test binding[58] == 0.0
    @test all(==(0.0), binding[observed(37:48)])
    @test all(==(1.5), label[observed(50:56)])
    @test label[40] == 1.0
    @test label[20] == 0.0
end
