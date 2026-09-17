# =============================================================================
# Burn-in convergence probe.
#
# Integrates the baseline (seasonal daily) burn-in with NO early stopping and
# reports, for every year block, the basin-total, q95 and max relative changes in
# site biomass.  Used to choose `spin_up_max_years` / `spin_up_tol` for the K = 1x
# ensemble and to prove the burn-in actually reaches a stationary state.
#
# Env overrides:
#   GUADEX_PROBE_YEARS      number of one-year blocks (default 100)
#   GUADEX_PROBE_K_BASE     base carrying-capacity multiplier (default 1.0)
#   GUADEX_PROBE_FORCING_FILE  wide daily forcing (default from parameters.toml)
#   GUADEX_PROBE_OUT        output directory (default results/climate_scenarios_k1x_burnin)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Guadex
using DataFrames
using CSV
using Statistics
using SparseArrays
using LinearAlgebra

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "parameters.jl"))
using .SimulationParameters
const PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "run_climate_scenarios")

function _arg(flag, default)
    idx = findfirst(==(flag), ARGS)
    idx === nothing && return default
    idx == length(ARGS) && error("$flag needs a value")
    return ARGS[idx + 1]
end

const DAYS_PER_YEAR = Int(PARAMS["general"]["days_per_year"])
const RC = PARAMS["run_climate_scenarios"]
const FORCING_FILE = _arg("--forcing", get(ENV, "GUADEX_PROBE_FORCING_FILE",
    string(get(RC, "daily_forcing_file", "guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide.csv"))))
const K_BASE = parse(Float64, _arg("--kbase", get(ENV, "GUADEX_PROBE_K_BASE", "1.0")))
const MAX_YEARS = parse(Int, _arg("--years", get(ENV, "GUADEX_PROBE_YEARS", "100")))
const TOL = parse(Float64, _arg("--tol", get(ENV, "GUADEX_PROBE_TOL", "1.0e-5")))
const OUT_DIR = _arg("--out", get(ENV, "GUADEX_PROBE_OUT", "results/climate_scenarios_k1x_burnin"))

println("="^70)
println("Burn-in probe: K_BASE=$(K_BASE)x observed, max $(MAX_YEARS) yr")
println("  forcing: $FORCING_FILE")
println("="^70)

# Honour the same input/obstacle and heat-stress overrides as
# run_climate_scenarios.jl so the probed burn-in is the same ODE as the ensemble.
const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE",
    get(PARAMS["inputs"], "cedex_var_file", nothing))
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE",
    get(PARAMS["inputs"], "cedex_esc_uts_file", nothing))
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE",
    get(PARAMS["inputs"], "obstacles_file", nothing))
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE",
    string(get(PARAMS["obstacles"], "mode", "legacy"))))
const OBSTACLE_TOLERANCE = parse(Float64, get(ENV, "GUADEX_OBSTACLE_TOLERANCE_M",
    string(get(PARAMS["obstacles"], "matching_tolerance_m", 2000.0))))
const OBSTACLE_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_PASSABILITY",
    string(get(PARAMS["obstacles"], "upstream_passability", 0.1))))
const OBSTACLE_DOWNSTREAM_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY",
    string(get(PARAMS["obstacles"], "downstream_passability", 0.5))))

const STRESS_PARAMS = get(PARAMS, "temperature_stress", Dict{String,Any}())
const HEAT_STRESS_ENABLED = begin
    cli = lowercase(_arg("--heat-stress", ""))
    cli in ("1", "true", "on", "yes") ? true :
    cli in ("0", "false", "off", "no") ? false :
    lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS",
        string(get(STRESS_PARAMS, "enabled", false)))) in ("1", "true", "yes")
end
const HEAT_STRESS_K = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_K",
    string(get(STRESS_PARAMS, "k", 0.0))))
const HEAT_STRESS_CALIBRATE = lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE",
    string(get(STRESS_PARAMS, "calibrate", false)))) in ("1", "true", "yes")
const HEAT_STRESS_MAX_LOSS = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_MAX_LOSS",
    string(get(STRESS_PARAMS, "max_annual_loss", 0.05))))

data = prepare_ode_data(
    upstream_cost = Float64(get(RC, "upstream_cost", 0.05)),
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_TOLERANCE,
    obstacle_passability = OBSTACLE_PASSABILITY,
    obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY,
    carrying_capacity_base_scaling = K_BASE
)

n_sites = data.params.n_sites
n_species = data.params.n_species
density_cols = [Symbol("$(sp)_DEN") for sp in data.species]
density_df = filter(row -> row.CODIGO in data.sites, data.density_df)
u0 = max.(replace(Matrix(density_df[:, density_cols]), NaN => 0.0), 0.0)
u0_flat = vec(u0)
println("Sites=$(n_sites), species=$(n_species); mean observed total = " *
        "$(round(sum(u0) / n_sites, digits=4))")

daily = load_daily_forcing_any(FORCING_FILE)
matrix_of(df, sc) = is_wide_daily_forcing(df) ?
    wide_forcing_matrix(df, data.sites; scenario=sc) :
    daily_forcing_matrix(df, data.sites; scenario=sc)
bdates, btemps = matrix_of(daily, "historical")
println("Baseline daily series: $(length(bdates)) days from 'historical'")

schedule = baseline_climatology_schedule(btemps, bdates, 1; days_per_year=Float64(DAYS_PER_YEAR))

k_heat = HEAT_STRESS_ENABLED ? HEAT_STRESS_K : 0.0
if HEAT_STRESS_ENABLED && HEAT_STRESS_CALIBRATE
    energies = Float64[]
    for s in 1:n_species, i in 1:n_sites
        push!(energies, exceedance_energy(@view(btemps[i, :]), data.params.thermal_upper_limits[s]))
    end
    k_heat = calibrate_heat_stress_rate(energies; max_annual_loss=HEAT_STRESS_MAX_LOSS)
    println("Heat-stress calibration: k=$k_heat (max annual loss $HEAT_STRESS_MAX_LOSS)")
end
println("Heat-stress mortality: $(k_heat > 0 ? "k=$k_heat" : "disabled")")
params = set_heat_stress_rate(data.params, k_heat)

println("\nIntegrating burn-in (tol=$TOL, criterion=:basin)...")
t0 = time()
spin = spin_up(params; initial_state=u0_flat, schedule=schedule,
    days_per_year=Float64(DAYS_PER_YEAR), max_years=MAX_YEARS, tol=TOL,
    criterion=:basin, progress_every=25)
elapsed = time() - t0
println("Integrated $(spin.years) year-blocks in $(round(elapsed, digits=1)) s " *
        "($(round(elapsed / max(spin.years, 1), digits=2)) s/yr); " *
        "converged=$(spin.converged)")

rows = NamedTuple[]
let prev = site_totals(u0_flat, n_sites, n_species)
    for (k, total) in enumerate(spin.history)
        rel = abs.(total .- prev) ./ max.(abs.(prev), 1e-12)
        push!(rows, (year=k,
            basin_change=abs(sum(total) - sum(prev)) / max(abs(sum(prev)), 1e-12),
            q95_change=quantile(rel, 0.95),
            max_change=maximum(rel),
            basin_total=sum(total)))
        prev = total
    end
end
df = DataFrame(rows)
mkpath(OUT_DIR)
out_path = joinpath(OUT_DIR, "burnin_probe.csv")
CSV.write(out_path, df)

checkpoints = sort(unique(vcat(1:min(10, MAX_YEARS), [10, 20, 30, 40, 50, 75, 100,
    150, 200, 300, 400, 500], [MAX_YEARS])))
println("\nyear   basin_change   q95_change   max_change   basin_total")
for k in checkpoints
    (1 <= k <= nrow(df)) || continue
    r = df[k, :]
    println(lpad(r.year, 4), "   ", lpad(round(r.basin_change, sigdigits=4), 12),
        "   ", lpad(round(r.q95_change, sigdigits=4), 9),
        "   ", lpad(round(r.max_change, sigdigits=4), 10),
        "   ", lpad(round(r.basin_total, digits=3), 11))
end
println("\nProbe written to $out_path")
