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

# Print example rows of the dataset CSVs (DisProt, PDB mapping, UniProt,
# PDB contacts, aggregated labels, BMRB mapping, chemical shift labels), count
# rows per category (quality_category, status, reason, state, used), count
# missing values, and summarise free / bound NMR shift lists (kept, compared,
# rejected and why).
#
# The rows shown are the first rows WITHOUT missing value, so that every column
# is filled; when a table has fewer complete rows, the rows with the fewest
# missing values are shown after them.


using CSV
using DataFrames

const N_ROWS = isempty(ARGS) ? 5 : parse(Int, ARGS[1])
const MAX_CELL_WIDTH = 40   # long cells (label vectors) are cut at this width
# Columns whose values are categories: the number of rows per value is printed.
const CATEGORY_COLUMNS = ["quality_category", "status", "reason", "state", "used", "sources"]

const CSV_FILES = [
    "DisProt (final dataset)" => "data/dataset/cluster_dataset_labels.csv.gz",
    "DisProt -> PDB mapping" => "data/derived/disprot_pdb_mapping/disprot_pdb_mapping.csv.gz",
    "UniProt binding sites" => "data/derived/uniprot_binding_sites/uniprot_binding_sites.csv.gz",
    "PDB contacts" => "data/derived/pdb_contacts/pdb_contacts.csv.gz",
    "Aggregated labels" => "data/derived/aggregated_labels/aggregated_labels.csv.gz",
    "Aggregated labels PDB log" => "data/derived/aggregated_labels/aggregated_labels_pdb_log.csv.gz",
    "DisProt -> BMRB mapping" => "data/derived/disprot_bmrb_mapping/disprot_bmrb_mapping.csv.gz",
    "Chemical shift labels" => "data/derived/chemical_shift_labels/chemical_shift_labels.csv.gz",
    "Chemical shift lists" => "data/derived/chemical_shift_labels/chemical_shift_lists.csv.gz",
    "Chemical shift free/bound comparisons" => "data/derived/chemical_shift_labels/chemical_shift_comparisons.csv.gz",
    "Final labels" => "data/derived/final_labels/final_labels.csv.gz",
    "Final label sources" => "data/derived/final_labels/final_labels_sources.csv.gz"
]

const BMRB_MAPPING_PATH = "data/derived/disprot_bmrb_mapping/disprot_bmrb_mapping.csv.gz"
const CS_LISTS_PATH = "data/derived/chemical_shift_labels/chemical_shift_lists.csv.gz"
const CS_COMPARISONS_PATH = "data/derived/chemical_shift_labels/chemical_shift_comparisons.csv.gz"
const CS_CRITERIA = ["pH", "temperature", "ionic", "buffer", "solvent", "window",
    "field", "spectrometer"]

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

"Print `label: n` lines with percentages of `total`."
function print_counts(pairs, total)
    for (label, n) in pairs
        println("  ", rpad(label, 46), lpad(n, 8),
            total > 0 ? "  ($(round(100 * n / total; digits = 1)) %)" : "")
    end
end

"Free / bound shift lists from the BMRB mapping to the chemical-shift labels."
function show_nmr_summary()
    println("\n", "="^80, "\nNMR free / bound summary\n", "="^80)
    if isfile(BMRB_MAPPING_PATH)
        mapping = DataFrame(CSV.File(BMRB_MAPPING_PATH; types = String))
        listed = mapping[coalesce.(mapping.status, "") .== "ok", :]
        println("Shift lists at the start (BMRB mapping, status = ok): ", nrow(listed))
        print_counts([string(coalesce(state, "missing")) => count(isequal(state), listed.state)
                      for state in unique(listed.state)], nrow(listed))
    else
        println("BMRB mapping not found: ", BMRB_MAPPING_PATH)
    end
    if !isfile(CS_LISTS_PATH)
        println("Chemical shift lists not found: ", CS_LISTS_PATH)
        return
    end
    lists = DataFrame(CSV.File(CS_LISTS_PATH; types = String))
    for state in ("free", "bound")
        rows = lists[coalesce.(lists.state, "") .== state, :]
        n_used = count(==("true"), rows.used)
        println("\n", uppercase(state), " shift lists: ", nrow(rows), " at the start, ",
            n_used, state == "free" ? " used for disorder" : " compared with a free state (pair used)")
        reasons = sort!(combine(groupby(rows[rows.used .!= "true", :], :reason), nrow => :n),
            :n; rev = true)
        nrow(reasons) == 0 || println("  set aside:")
        print_counts([row.reason => row.n for row in eachrow(reasons)], nrow(rows))
    end
    other = lists[.!in.(coalesce.(lists.state, ""), Ref(("free", "bound"))), :]
    nrow(other) == 0 || (println("\nOther shift lists (never used):");
     print_counts([string(r.reason) => r.n
                   for r in eachrow(combine(groupby(other, :reason), nrow => :n))],
         nrow(other)))

    isfile(CS_COMPARISONS_PATH) || return
    comparisons = DataFrame(CSV.File(CS_COMPARISONS_PATH; types = String))
    rejected = comparisons[comparisons.accepted .!= "true", :]
    println("\nBound x free comparisons: ", nrow(comparisons), " (accepted: ",
        nrow(comparisons) - nrow(rejected), ", rejected: ", nrow(rejected), ")")
    println("Rejected comparisons per criterion (one comparison can differ on several):")
    print_counts(["different $(c)" => count(==("different"), rejected[!, c])
                  for c in CS_CRITERIA if c in names(rejected)], nrow(rejected))
    println("Criteria not deposited (unknown) in all comparisons:")
    print_counts(["unknown $(c)" => count(==("unknown"), comparisons[!, c])
                  for c in CS_CRITERIA if c in names(comparisons)], nrow(comparisons))
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

show_nmr_summary()
