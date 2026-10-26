# =============================================================================
# Robust burn-in composition probe (E4 fix).
#
# Integrates the baseline (seasonal daily) burn-in with NO early stopping and
# records, for every year block, the raw all-cells composition statistic
# (95th percentile |Δ log N| over every site × species cell) alongside candidate
# *robust* statistics that restrict the percentile to cells above a density
# floor in at least one of the two compared years.  Used to choose the floor,
# the composition tolerance and min_years for the corrected config.
#
# Env overrides:
#   GUADEX_PARAMETERS_FILE   parameter file (default parameters.toml)
#   GUADEX_PROBE_YEARS       number of one-year blocks (default 120)
#   GUADEX_PROBE_K_BASE      base carrying-capacity multiplier (default 1.0)
#   GUADEX_PROBE_FORCING_FILE  wide daily forcing (default from parameters)
#   GUADEX_PROBE_OUT         output directory (default results/probe_burnin_robust)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Guadex
using DifferentialEquations
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

const DAYS_PER_YEAR = Int(PARAMS["general"]["days_per_year"])
const RC = PARAMS["run_climate_scenarios"]
const MAX_YEARS = parse(Int, get(ENV, "GUADEX_PROBE_YEARS", "120"))
const K_BASE = parse(Float64, get(ENV, "GUADEX_PROBE_K_BASE", "1.0"))
const FORCING_FILE = get(ENV, "GUADEX_PROBE_FORCING_FILE",
    string(get(RC, "daily_forcing_file", "guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide.csv")))
const OUT_DIR = get(ENV, "GUADEX_PROBE_OUT", joinpath("results", "probe_burnin_robust"))

const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE", get(PARAMS["inputs"], "cedex_var_file", nothing))
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE", get(PARAMS["inputs"], "cedex_esc_uts_file", nothing))
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE", get(PARAMS["inputs"], "obstacles_file", nothing))
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE", string(get(PARAMS["obstacles"], "mode", "legacy"))))
const OBSTACLE_TOLERANCE = parse(Float64, get(ENV, "GUADEX_OBSTACLE_TOLERANCE_M",
    string(get(PARAMS["obstacles"], "matching_tolerance_m", 2000.0))))
const OBSTACLE_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_PASSABILITY",
    string(get(PARAMS["obstacles"], "upstream_passability", 0.1))))
const OBSTACLE_DOWNSTREAM_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY",
    string(get(PARAMS["obstacles"], "downstream_passability", 0.5))))

const STRESS_PARAMS = get(PARAMS, "temperature_stress", Dict{String,Any}())
const HEAT_STRESS_ENABLED = lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS",
    string(get(STRESS_PARAMS, "enabled", false)))) in ("1", "true", "yes")
const HEAT_STRESS_K = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_K",
    string(get(STRESS_PARAMS, "k", 0.0))))
const HEAT_STRESS_CALIBRATE = lowercase(get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE",
    string(get(STRESS_PARAMS, "calibrate", false)))) in ("1", "true", "yes")
const HEAT_STRESS_MAX_LOSS = parse(Float64, get(ENV, "GUADEX_CLIMATE_HEAT_STRESS_MAX_LOSS",
    string(get(STRESS_PARAMS, "max_annual_loss", 0.05))))

const FLOORS = [1e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]

println("="^70)
println("Robust composition probe: max $(MAX_YEARS) yr, K_BASE=$(K_BASE)x")
println("  parameter file: $(SimulationParameters.default_file())")
println("  forcing: $FORCING_FILE")
println("="^70)

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
        "$(round(sum(u0) / n_sites, digits=4)); max observed cell = $(round(maximum(u0), digits=4))")

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
params = set_heat_stress_rate(data.params, k_heat)
println("Heat-stress mortality: $(k_heat > 0 ? "k=$k_heat" : "disabled")")

# ---- Manual burn-in loop (same solver/callback as Guadex.spin_up) ----------
clamp_cb = DiscreteCallback(
    (u, t, integrator) -> any(x -> !isfinite(x) || x < 0, u),
    integrator -> begin
        @inbounds for i in eachindex(integrator.u)
            x = integrator.u[i]
            integrator.u[i] = (isfinite(x) && x > 0) ? x : 0.0
        end
    end;
    save_positions=(false, false))
rhs_params = ScheduledMetacommunityParams(params, schedule)

function robust_stats(prev, cur, floor; log_floor=1e-12)
    deltas = Float64[]
    n_active = 0
    for i in eachindex(prev)
        a = Float64(prev[i]); b = Float64(cur[i])
        max(a, b) > floor || continue
        n_active += 1
        push!(deltas, abs(log(max(b, log_floor)) - log(max(a, log_floor))))
    end
    return (q95 = isempty(deltas) ? 0.0 : quantile(deltas, 0.95), n_active = n_active)
end

function run_probe(u0_flat)
    rows = NamedTuple[]
    prev_state = copy(u0_flat)
    prev_total = site_totals(u0_flat, n_sites, n_species)
    for year in 1:MAX_YEARS
        prob = ODEProblem(metacommunity_ode_scheduled!, prev_state, (0.0, Float64(DAYS_PER_YEAR)), rhs_params)
        sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6,
            save_everystep=false, save_start=false, save_end=true, callback=clamp_cb)
        cur_state = vec(Float64.(sol.u[end]))
        total = site_totals(cur_state, n_sites, n_species)
        basin = abs(sum(total) - sum(prev_total)) / max(abs(sum(prev_total)), 1e-12)
        raw = robust_stats(prev_state, cur_state, 0.0; log_floor=1e-12).q95
        entry = Dict{Symbol,Any}(:year => year, :basin_change => basin, :all_cells_q95 => raw,
            :frac_below_1e_6 => count(x -> x <= 1e-6, cur_state) / length(cur_state),
            :n_above_1e_2 => count(x -> x > 1e-2, cur_state),
            :n_above_1e_1 => count(x -> x > 1e-1, cur_state))
        for f in FLOORS
            st = robust_stats(prev_state, cur_state, f)
            entry[Symbol("q95_floor_$(f)")] = st.q95
            entry[Symbol("n_active_$(f)")] = st.n_active
        end
        push!(rows, (; entry...))
        prev_state = cur_state
        prev_total = total
        if year % 10 == 0
            println("  yr $(lpad(year,3)): basin=$(round(basin, sigdigits=4)) " *
                    "raw=$(round(raw, sigdigits=4)) " *
                    "robust[1e-3]=$(round(entry[Symbol("q95_floor_$(1e-3)")], sigdigits=4)) " *
                    "robust[1e-2]=$(round(entry[Symbol("q95_floor_$(1e-2)")], sigdigits=4)) " *
                    "n>1e-2=$(entry[:n_above_1e_2]) n>1e-1=$(entry[:n_above_1e_1])")
        end
    end
    return rows
end

t0 = time()
rows = run_probe(u0_flat)
elapsed = time() - t0
println("Integrated $(MAX_YEARS) year-blocks in $(round(elapsed, digits=1)) s " *
        "($(round(elapsed / MAX_YEARS, digits=2)) s/yr)")

df = DataFrame(rows)
mkpath(OUT_DIR)
out_path = joinpath(OUT_DIR, "burnin_robust_probe.csv")
CSV.write(out_path, df)
println("Probe written to $out_path")
