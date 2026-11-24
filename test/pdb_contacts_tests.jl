using TestItems

@testitem "PDB quality table parsing" begin
    using Test
    using IDRBindingDataset: _parse_resolution, _pdb_quality_category,
                             read_pdb_quality_table

    @test _parse_resolution("2.10") == 2.1
    @test _parse_resolution("1.9, 2.3") == 1.9
    @test _parse_resolution(3) == 3.0
    @test _parse_resolution("NOT APPLICABLE") === missing
    @test _parse_resolution(missing) === missing

    @test _pdb_quality_category(nothing) === :absent
    @test _pdb_quality_category((method = "SOLUTION NMR", resolution = missing)) ===
          :no_resolution
    @test _pdb_quality_category((method = "X-RAY DIFFRACTION", resolution = 2.0)) ===
          :resolved

    path = joinpath(mktempdir(), "pdb_info.csv")
    write(path, "IDCODE,METHOD,RESOLUTION\n1ABC,x-ray diffraction,2.5\n2XYZ,,NOT APPLICABLE\n")
    quality = read_pdb_quality_table(path)
    @test quality["1abc"] == (method = "X-RAY DIFFRACTION", resolution = 2.5)
    @test quality["2xyz"].method == "missing"
    @test quality["2xyz"].resolution === missing
end

@testitem "Homologous chains are detected by k-mers" begin
    using Test
    using IDRBindingDataset: _kmer_fraction, _homologous_chains

    protein = "MKTAYIAKQRQISFVKSHFSRQ"
    @test _kmer_fraction("MKTAYIAKQR", protein; k = 5) == 1.0
    @test _kmer_fraction("WWWWWWWW", protein; k = 5) == 0.0
    @test _kmer_fraction("MKT", protein; k = 5) == 1.0
    @test _kmer_fraction("", protein) == 0.0

    chains = Dict("A" => "MKTAYIAKQR", "B" => "QISFVKSHFS", "C" => "GGGGGGGGGG")
    @test _homologous_chains(chains, ["A", "Z"], protein) == Set(["A", "B"])
end

@testitem "Interface residues are mapped to UniProt positions" begin
    using Test
    using IDRBindingDataset: _interface_residue_indices, _chain_contact_positions
    using IDRBindingDataset.MIToS.PDB: PDBResidue, PDBResidueIdentifier, PDBAtom,
                                       Coordinates

    atom(x) = PDBAtom(Coordinates(x, 0.0, 0.0), "CA", "C", 1.0, "0.0", "", "")
    residue(chain, number, name, x) = PDBResidue(
        PDBResidueIdentifier(number, number, name, "ATOM", "1", chain), [atom(x)])

    query = [residue("A", "1", "MET", 0.0), residue("A", "2", "LYS", 20.0),
        residue("A", "3", "GLY", 40.0), residue("A", "4", "THR", 60.0)]
    partner = [residue("B", "1", "GLY", 4.0), residue("B", "2", "GLY", 45.0),
        residue("B", "3", "GLY", 63.0)]

    hits = _interface_residue_indices(query, partner, 6.0)
    @test hits == Set([1, 3, 4])

    # PDB residue 4 has no UniProt mapping; residue 3 is GLY in the PDB but T in UniProt.
    mapping = Dict("1" => 1, "2" => 2, "3" => 3)
    result = _chain_contact_positions(query, mapping, hits, "MKT")
    @test result.resolved == Set([1, 2, 3])
    @test result.contacts == Set([1, 3])
    @test result.n_mapped == 3
    @test result.n_identical == 2
    @test result.mismatches == Set([3])
end

@testitem "Positions are encoded as compact ranges" begin
    using Test
    using IDRBindingDataset: _format_position_ranges, _parse_residue_ranges

    @test _format_position_ranges([7, 1, 2, 3, 3, 9, 10]) == "1-3;7;9-10"
    @test _format_position_ranges([5]) == "5"
    @test _format_position_ranges(Int[]) === missing
    @test _parse_residue_ranges(_format_position_ranges([1, 2, 3, 7, 9, 10])) ==
          [1, 2, 3, 7, 9, 10]
end

@testitem "Contacts table selects binding proteins and reports missing files" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_pdb_contacts_table

    dataset = DataFrame(protein_id = ["P1", "P2"], protein_disprot_id = ["DP1", "DP2"],
        sequence = ["MKT", "MKT"], residue_labels = ["011", "012"])
    mapping = DataFrame(protein_disprot_id = ["DP1", "DP2", "DP2"],
        protein_id = ["P1", "P2", "P2"], PDB = ["1abc", "9zzz", "9zzz"],
        CHAIN = ["A", "A", "B"])

    result = build_pdb_contacts_table(dataset, mapping;
        mmcif_dirs = [mktempdir()], sifts_xml_dirs = [mktempdir()],
        logger = _ -> nothing)
    table = result.table

    # P1 has no DisProt binding residue (label 2); P2 has two rows for one PDB.
    @test table.protein_id == ["P2"]
    @test table.protein_disprot_id == ["DP2"]
    @test table.status == ["no_mmcif"]
    @test table.quality_category == ["absent"]
    @test ismissing(only(table.contact_positions))
    @test isempty(result.failed_items)

    all_proteins = build_pdb_contacts_table(dataset, mapping;
        mmcif_dirs = [mktempdir()], sifts_xml_dirs = [mktempdir()],
        only_binding_proteins = false, logger = _ -> nothing)
    @test all_proteins.table.protein_id == ["P1", "P2"]
end

@testitem "Contacts table on a real mmCIF and SIFTS file" begin
    using Test
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: build_pdb_contacts_table, _read_sifts_residue_records
    import IDRBindingDataset.MIToS

    # MIToS test files: PDB 2VQC (UniProt P20220), a homo-oligomer.
    data_dir = normpath(joinpath(dirname(pathof(MIToS)), "..", "test", "data"))
    records = _read_sifts_residue_records(joinpath(data_dir, "2vqc.xml.gz"))
    mapped = [r for r in records["A"] if r.uniprot_id == "P20220"]
    sequence_length = maximum(r.uniprot_position for r in mapped)
    sequence = fill('X', sequence_length)
    for r in mapped
        r.is_missing && continue
        sequence[r.uniprot_position] = 'A'
    end

    dataset = DataFrame(protein_id = ["P20220"], protein_disprot_id = ["DP00847"],
        sequence = [String(sequence)], residue_labels = [repeat("2", sequence_length)])
    mapping = DataFrame(protein_disprot_id = ["DP00847"], protein_id = ["P20220"],
        PDB = ["2vqc"], CHAIN = ["A"])

    table = build_pdb_contacts_table(dataset, mapping;
        mmcif_dirs = [data_dir], sifts_xml_dirs = [data_dir],
        logger = _ -> nothing).table

    @test nrow(table) == 1
    @test only(table.status) in ("ok", "no_heteromeric_partner")
    @test only(table.chains) == "A"
    @test !("n_mapped" in names(table))
    @test !("n_resolved" in names(table))
    @test !("n_contacts_unmapped" in names(table))
    @test only(table.n_contacts) == 0 || only(table.status) == "ok"
    @test !ismissing(only(table.resolved_positions))
end

@testitem "mmCIF structures are read with author chain identifiers" begin
    using Test
    using IDRBindingDataset: _load_first_model_heavy_atoms
    import IDRBindingDataset.MIToS

    # In 2VQC, waters and ligands have their own label chain (B) but the author
    # chain is A, as in SIFTS.
    path = normpath(joinpath(dirname(pathof(MIToS)), "..", "test", "data", "2vqc.cif"))
    label_chains = unique(r.id.chain
    for r in MIToS.PDB.read_file(path, MIToS.PDB.MMCIFFile))
    @test "B" in label_chains
    residues = _load_first_model_heavy_atoms(path)
    @test unique(r.id.chain for r in residues) == ["A"]
    @test all(atom -> !(uppercase(strip(atom.element)) in ("H", "D")),
        atom for residue in residues for atom in residue.atoms)
end
