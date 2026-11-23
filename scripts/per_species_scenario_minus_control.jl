using Pkg; Pkg.activate(".")
using DataFrames
using CSV
using Statistics

# Per-species scenario-minus-control change in basin total density at end-of-2045.
#
# Reads only the already-exported Phase 2 tables (no simulation) and writes
# results/climate_scenarios_phase2/per_species_scenario_minus_control.csv, the
# source of Table 6 in docs/FinalReportOctober2026.tex.
#
#   julia --project=. scripts/per_species_scenario_minus_control.jl
#
# Override the root with GUADEX_CLIMATE_OUTPUT_DIR (default:
# results/climate_scenarios_phase2).

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "scripts", "plot_paths.jl"))

const CLIMATE_ROOT = report_output_root("GUADEX_CLIMATE_OUTPUT_DIR",
    joinpath("results", "climate_scenarios_phase2"); base=ROOT)

run_path(row) = normpath(joinpath(ROOT, replace(String(row.run_dir), "\\" => "/")))

function total_density_2045(run_dir)
    path = joinpath(run_dir, "export", "levels", "species_timeseries.csv")
    isfile(path) || return Dict{String,Float64}()
    df = CSV.read(path, DataFrame)
    sub = df[df.year .== maximum(Int.(df.year)), :]
    totals = Dict{String,Float64}()
    for row in eachrow(sub)
        totals[string(row.species)] = get(totals, string(row.species), 0.0) +
                                      Float64(row.density)
    end
    return totals
end

function main()
    index = CSV.read(joinpath(CLIMATE_ROOT, "runs_index.csv"), DataFrame)
    control_rows = index[index.scenario .== "control", :]
    nrow(control_rows) == 1 || error("expected exactly one control run")
    control = total_density_2045(run_path(control_rows[1, :]))
    scenarios = index[index.scenario .!= "control", :]
    run_totals = [(r, total_density_2045(run_path(r))) for r in eachrow(scenarios)]
    species = sort(collect(keys(control)))

    out = DataFrame(species=String[], control_total_density_2045=Float64[],
        mean_delta=Float64[], median_delta=Float64[], min_delta=Float64[],
        max_delta=Float64[], n_positive=Int[], n_negative=Int[], n_zero=Int[],
        mean_relative_change_pct=Float64[])
    for sp in species
        deltas = Float64[]
        rels = Float64[]
        for (_, totals) in run_totals
            delta = get(totals, sp, 0.0) - get(control, sp, 0.0)
            push!(deltas, delta)
            push!(rels, control[sp] == 0 ? NaN : 100 * delta / control[sp])
        end
        push!(out, (sp, control[sp], mean(deltas), median(deltas), minimum(deltas),
            maximum(deltas), count(>(0), deltas), count(<(0), deltas),
            count(==(0), deltas), mean(filter(isfinite, rels))))
    end

    path = joinpath(CLIMATE_ROOT, "per_species_scenario_minus_control.csv")
    CSV.write(path, out)
    println("wrote $path ($(nrow(out)) species, $(nrow(scenarios)) runs)")
end

main()
