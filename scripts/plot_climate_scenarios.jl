using Pkg; Pkg.activate(".")
using Guadex

# =============================================================================
# --- Climate-scenario figures (plot only, no re-simulation) ---
#
# Reads the four-level CSV outputs under `results/climate_scenarios` and writes
# per-run, ensemble (with across-GCM confidence bands), across-scenario and
# final-year figures.  Safe to run at any time; it never re-runs the model.
#
# Usage:
#   julia --project=. scripts/plot_climate_scenarios.jl [results_root] [figures_dir]
#
# Defaults:
#   results_root = results/climate_scenarios
#   figures_dir  = <results_root>/figures
#
# Environment:
#   GUADEX_CONFIG_ONLY=1   validate arguments/paths and exit
#   GUADEX_CLIMATE_END_YEAR  override the expected final year (default: max end_year in runs_index.csv)
# =============================================================================

if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] plot configuration smoke test passed")
    exit(0)
end

results_root = length(ARGS) >= 1 ? ARGS[1] : joinpath("results", "climate_scenarios")
figures_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(results_root, "figures")

if !isdir(results_root)
    error("results root not found: $results_root (run run_climate_scenarios.jl first)")
end

expected_end_year = get(ENV, "GUADEX_CLIMATE_END_YEAR", "")
expected_end_year = isempty(expected_end_year) ? nothing : parse(Int, expected_end_year)

println("="^70)
println("Climate-scenario figures")
println("  results:  $results_root")
println("  figures:  $figures_dir")
println("  expected end year: $(expected_end_year === nothing ? "auto" : expected_end_year)")
println("="^70)

path = plot_climate_scenario_figures(results_root;
    figures_dir=figures_dir,
    expected_end_year=expected_end_year)

println("\nDone. Figures written to: $path")
