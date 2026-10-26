# =============================================================================
# build_interaction_long_table.jl
#
# Regenerate the explicit long-format species interaction table (issue C1) from
# the semicolon-delimited qualitative matrix.
#
# Usage:
#   julia --project=. scripts/build_interaction_long_table.jl
#
#   Optional arguments:
#     --input PATH    source matrix CSV
#                     (default data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv)
#     --output PATH   long-format CSV destination
#                     (default data/BIOTIC/interaction_matrix_long.csv)
#     -h, --help      show help
#
# The output has columns source,target,mechanism,alpha,ambiguous,raw_text with
# one row per directed effect (an undirected symmetric relationship gives two
# rows).  Cells that mix an undirected clause with a directed clause are flagged
# ambiguous = true.  "No coexist" (allopatry) contributes no rows (alpha = 0).
# =============================================================================

using Pkg; Pkg.activate(".")
using CSV, DataFrames, Guadex

const DEFAULT_INPUT = "data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv"
const DEFAULT_OUTPUT = "data/BIOTIC/interaction_matrix_long.csv"

function parse_args()
    input = DEFAULT_INPUT
    output = DEFAULT_OUTPUT
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--input"
            i += 1; i <= length(ARGS) || error("--input requires a value")
            input = ARGS[i]
        elseif startswith(arg, "--input=")
            input = arg[9:end]
        elseif arg == "--output"
            i += 1; i <= length(ARGS) || error("--output requires a value")
            output = ARGS[i]
        elseif startswith(arg, "--output=")
            output = arg[10:end]
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/build_interaction_long_table.jl [--input PATH] [--output PATH]")
            exit(0)
        else
            error("Unknown argument: $arg (use --help for usage)")
        end
        i += 1
    end
    return input, output
end

function count_ambiguous_cells(input::String)
    df = CSV.read(input, DataFrame; delim=';')
    rename!(df, 1 => :Species)
    labels = string.(df.Species)
    cols = String.(names(df)[2:end])
    n = 0
    for (i, row) in enumerate(eachrow(df)), c in cols
        effects = Guadex.parse_interaction_cell(row[Symbol(c)], labels[i], c)
        if !isempty(effects) && any(e -> e.ambiguous, effects)
            n += 1
        end
    end
    return n
end

function main()
    input, output = parse_args()
    println("Reading interaction matrix: $input")
    table = Guadex.build_interaction_long_table(input)

    println("Long-format rows: $(nrow(table))")
    println("Ambiguous cells:  $(count_ambiguous_cells(input))")
    println("Ambiguous rows:   $(count(table.ambiguous))")
    println("Non-zero alphas:  $(count(x -> x != 0.0, table.alpha))")

    mkpath(dirname(output))
    CSV.write(output, table)
    println("Wrote $output")
end

main()
