#!/home/julie.daniel/.julia/juliaup/julia-1.11.7+0.x64.linux.gnu/bin/julia

#=
#SBATCH --nodelist=node48
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --output=check_csv_created.jl.o%j.out
=#

import Pkg
Pkg.activate("/stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDP/scripts/env/env1_020226")

# Print the first rows of the dataset CSVs (DisProt, PDB mapping, UniProt,
# PDB contacts, aggregated labels).
#
# The rows shown are the first rows WITHOUT missing value, so that every column
# is filled; when a table has fewer complete rows, the rows with the fewest
# missing values are shown after them.

using CSV
using DataFrames

const N_ROWS = isempty(ARGS) ? 5 : parse(Int, ARGS[1])
const MAX_CELL_WIDTH = 40   # long cells (label vectors) are cut at this width
# Columns whose values are categories: the number of rows per value is printed.
const CATEGORY_COLUMNS = ["quality_category", "status", "reason"]

const CSV_FILES = [
    "DisProt (final dataset)" => "data/dataset/cluster_dataset_labels.csv.gz",
    "DisProt -> PDB mapping" => "data/derived/disprot_pdb_mapping/disprot_pdb_mapping.csv.gz",
    "UniProt binding sites" => "data/derived/uniprot_binding_sites/uniprot_binding_sites.csv.gz",
    "PDB contacts" => "data/derived/pdb_contacts/pdb_contacts.csv.gz",
    "Aggregated labels" => "data/derived/aggregated_labels/aggregated_labels.csv.gz",
    "Aggregated labels PDB log" => "data/derived/aggregated_labels/aggregated_labels_pdb_log.csv.gz"
]

"Print the number of rows per value of the category columns present in `table`."
function show_category_counts(table::DataFrame)
    for column in CATEGORY_COLUMNS
        column in names(table) || continue
        counts = sort!(combine(groupby(table, column), nrow => :n), :n; rev = true)
        println("Rows per `", column, "`:")
        for row in eachrow(counts)
            value = coalesce(row[column], "missing")
            println("  ", rpad(value, 32), lpad(row.n, 8), "  (",
                round(100 * row.n / nrow(table); digits = 1), " %)")
        end
    end
end

"Print how many rows have at least one missing value, and missing counts per column."
function show_missing_summary(table::DataFrame)
    n_missing_per_row = [count(ismissing, row) for row in eachrow(table)]
    n_rows_with_missing = count(>(0), n_missing_per_row)
    println("Rows with at least one missing value: ", n_rows_with_missing, " / ",
        nrow(table), " (complete rows: ", nrow(table) - n_rows_with_missing, ")")
    columns_with_missing = [(name, count(ismissing, table[!, name]))
                            for name in names(table)]
    filter!(pair -> last(pair) > 0, columns_with_missing)
    if isempty(columns_with_missing)
        println("No missing value in any column.")
    else
        println("Missing values per column:")
        for (name, n) in columns_with_missing
            println("  ", rpad(name, 32), n)
        end
    end
    n_missing_per_row
end

for (name, path) in CSV_FILES
    println("\n", "="^80, "\n", name, "  (", path, ")\n", "="^80)
    if !isfile(path)
        println("File not found.")
        continue
    end
    # Every column read as text: keeps label strings such as "0000" intact.
    table = DataFrame(CSV.File(path; types = String))
    println(nrow(table), " rows × ", ncol(table), " columns")
    show_category_counts(table)
    n_missing_per_row = show_missing_summary(table)

    # Complete rows first (stable sort keeps file order among equals).
    order = sortperm(n_missing_per_row; alg = MergeSort)
    examples = table[order[1:min(N_ROWS, nrow(table))], :]
    n_complete = count(i -> n_missing_per_row[i] == 0, order[1:nrow(examples)])
    println("\nExample rows (", n_complete, " without missing value",
        n_complete < nrow(examples) ? ", then rows with the fewest missing values" : "",
        "):")
    show(stdout, examples; allcols = true, truncate = MAX_CELL_WIDTH)
    println()
end