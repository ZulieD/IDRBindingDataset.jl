using TestItems

@testitem "UniProt binding sites are extracted from JSON entries" begin
    using Test
    using JSON
    using IDRBindingDataset: _extract_uniprot_binding_sites

    entry = JSON.parse("""
    {"primaryAccession": "P1", "features": [
      {"type": "Binding site", "location": {"start": {"value": 2}, "end": {"value": 4}},
       "evidences": [{"source": "PubMed"}, {"source": "PDB"}, {"source": "PubMed"}]},
      {"type": "Binding site", "location": {"start": {"value": 6}, "end": {"value": 6}}},
      {"type": "Binding site", "location": {"start": {"value": null}, "end": {"value": 9}}},
      {"type": "Region", "location": {"start": {"value": 1}, "end": {"value": 3}}}
    ]}
    """)

    sites = _extract_uniprot_binding_sites(entry)
    @test sites ==
          [(start = 2, stop = 4, n_pubmed = 2), (start = 6, stop = 6, n_pubmed = 0)]
    @test isempty(_extract_uniprot_binding_sites(nothing))
end

@testitem "UniProt residue ranges round-trip" begin
    using Test
    using IDRBindingDataset: _format_residue_ranges, _parse_residue_ranges

    sites = [(start = 12, stop = 15, n_pubmed = 0), (start = 30, stop = 30, n_pubmed = 1)]
    @test _format_residue_ranges(sites) == "12-15;30"
    @test _parse_residue_ranges("12-15;30") == [12, 13, 14, 15, 30]
    @test _parse_residue_ranges("[\"12-15\", \"30\"]") == [12, 13, 14, 15, 30]
    @test _parse_residue_ranges(missing) == Int[]
end

@testitem "UniProt entries are fetched in batches and summarized per protein" begin
    using Test
    using HTTP
    using JSON
    using IDRBindingDataset.DataFrames
    using IDRBindingDataset: fetch_uniprot_entries, build_uniprot_binding_sites_table

    requested_queries = String[]
    function mock_get(url, query)
        push!(requested_queries, query["query"])
        payload = Dict("results" => [
            Dict("primaryAccession" => "P1",
            "features" => [
                Dict("type" => "Binding site",
                    "location" => Dict("start" => Dict("value" => 2),
                        "end" => Dict("value" => 3)),
                    "evidences" => [Dict("source" => "PubMed")]),
                Dict("type" => "Binding site",
                    "location" => Dict("start" => Dict("value" => 5),
                        "end" => Dict("value" => 5)))
            ])
        ])
        HTTP.Response(200, ["X-UniProt-Release" => "2026_04"], JSON.json(payload))
    end

    fetched = fetch_uniprot_entries(["p2", "P1"]; batch_size = 1, base_sleep = 0.0,
        get_request = mock_get, sleeper = _ -> nothing, logger = _ -> nothing)
    @test requested_queries == ["accession:P1", "accession:P2"]
    @test collect(keys(fetched.entries)) == ["P1"]
    @test fetched.release == "2026_04"

    dataset = DataFrame(protein_id = ["P2", "P1"], protein_disprot_id = ["DP2", "DP1"],
        residue_labels = ["0000", "01200"])
    table = build_uniprot_binding_sites_table(dataset, fetched.entries)

    @test table.protein_disprot_id == ["DP1", "DP2"]
    @test table.protein_id == ["P1", "P2"]
    @test table.uniprot_entry_found == [true, false]
    @test table.binding_in_uniprot == [true, false]
    # The site at position 5 is on an ordered residue (label 0) and is dropped.
    @test isequal(table.residues_binding_in_uniprot, ["2-3", missing])
    @test isequal(table.publication_mention_binding, ["1", missing])
end

@testitem "UniProt fetch fails loudly after retries" begin
    using Test
    using HTTP
    using IDRBindingDataset: fetch_uniprot_entries

    failing_get(url, query) = HTTP.Response(500, "error")
    @test_throws ErrorException fetch_uniprot_entries(["P1"]; max_retries = 2,
        base_sleep = 0.0, get_request = failing_get, sleeper = _ -> nothing,
        logger = _ -> nothing)
end
