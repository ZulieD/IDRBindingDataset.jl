using TestItems

@testitem "PDB entry identifiers and file lookup" begin
    using Test
    using IDRBindingDataset: _pdb_entry_id, _find_pdb_entry_file, _sifts_xml_path

    @test _pdb_entry_id("1ABC") == "1abc"
    @test _pdb_entry_id(" 1abc_A ") == "1abc"
    @test _pdb_entry_id("not_pdb") == "not_pdb"

    flat_dir, nested_dir = mktempdir(), mktempdir()
    mkpath(joinpath(nested_dir, "bc"))
    nested_path = joinpath(nested_dir, "bc", "1bcd.xml.gz")
    flat_path = joinpath(flat_dir, "2xyz.xml")
    write(nested_path, "x")
    write(flat_path, "x")

    extensions = (".xml.gz", ".xml")
    @test _find_pdb_entry_file([flat_dir, nested_dir], "1BCD", extensions) == nested_path
    @test _find_pdb_entry_file([flat_dir, nested_dir], "2xyz", extensions) == flat_path
    @test _find_pdb_entry_file([flat_dir], "9zzz", extensions) === nothing
    @test _sifts_xml_path("9zzz"; sifts_xml_dirs = [flat_dir]) === nothing
end

@testitem "SIFTS residue records map PDB numbers to UniProt positions" begin
    using Test
    using IDRBindingDataset: SIFTSResidueRecord, _pdb_to_uniprot_positions

    records = SIFTSResidueRecord[
        (chain = "A", pdb_number = "10", uniprot_id = "P12345",
            uniprot_position = 5, is_missing = false),
        (chain = "A", pdb_number = "11", uniprot_id = "P12345",
            uniprot_position = 6, is_missing = true),
        (chain = "A", pdb_number = "12A", uniprot_id = "P12345",
            uniprot_position = 7, is_missing = false),
        (chain = "A", pdb_number = "13", uniprot_id = "",
            uniprot_position = 0, is_missing = false),
        (chain = "A", pdb_number = "14", uniprot_id = "Q99999",
            uniprot_position = 8, is_missing = false)
    ]

    @test _pdb_to_uniprot_positions(records, "p12345") == Dict("10" => 5, "12A" => 7)
    @test _pdb_to_uniprot_positions(records, "O00000") == Dict{String, Int}()
end

@testitem "Chain lists are encoded with semicolons" begin
    using Test
    using IDRBindingDataset: _format_chain_list, _parse_chain_list

    @test _format_chain_list(["A", "B"]) == "A;B"
    @test _parse_chain_list("A;B") == ["A", "B"]
    @test _parse_chain_list("A") == ["A"]
    @test _parse_chain_list("[\"A\", \"BB\"]") == ["A", "BB"]
    @test _parse_chain_list(missing) == String[]
end
