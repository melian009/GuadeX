using Pkg; Pkg.activate(".")
using Guadex

# =============================================================================
# --- Climate-scenario diagnostic figures (plot only, no re-simulation) ---
#
# Explains the shape of the climate ensemble rather than describing every run:
#
#   diagnostics/forcing_and_response.png  warming forcing vs richness/biomass
#   diagnostics/thermal_niches.png        species optima vs site temperatures
#   diagnostics/community_filling.png     initial vs final per-site richness
#
# Usage:
#   julia --project=. scripts/plot_climate_diagnostics.jl [results_root] [figures_dir]
# =============================================================================

if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] climate diagnostics configuration smoke test passed")
    exit(0)
end

results_root = length(ARGS) >= 1 ? ARGS[1] : joinpath("results", "climate_scenarios")
figures_dir = length(ARGS) >= 2 ? ARGS[2] : joinpath(results_root, "figures")

if !isdir(results_root)
    error("results root not found: $results_root (run run_climate_scenarios.jl first)")
end

println("="^70)
println("Climate-scenario diagnostics")
println("  results:  $results_root")
println("  figures:  $figures_dir")
println("="^70)

written = plot_climate_diagnostics(results_root; figures_dir=figures_dir)

println("\nDone. $(length(written)) diagnostic figure(s) written.")
for path in written
    println("  $path")
end
