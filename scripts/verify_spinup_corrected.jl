# =============================================================================
# Verify the corrected E4 burn-in on REAL data using the production `spin_up`.
#
# Runs the exact corrected stop rule (criterion=:both, tol=5e-5,
# composition_tol=1e-2, min_years=10, active floor 0.1) under the corrected
# config and reports the year it stops and both composition diagnostics.
#
# Env overrides:
#   GUADEX_PARAMETERS_FILE   default parameters_climate_scenarios_corrected.toml
#   GUADEX_VERIFY_MAX_YEARS  default 600
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Guadex
using DataFrames
using Statistics

const ROOT = dirname(@__DIR__)
include(joinpath(ROOT, "parameters.jl"))
using .SimulationParameters
get!(ENV, "GUADEX_PARAMETERS_FILE", "parameters_climate_scenarios_corrected.toml")
const PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "run_climate_scenarios")

const DAYS_PER_YEAR = Int(PARAMS["general"]["days_per_year"])
const RC = PARAMS["run_climate_scenarios"]
const MAX_YEARS = parse(Int, get(ENV, "GUADEX_VERIFY_MAX_YEARS", "600"))
const TOL = 5.0e-5
const COMPOSITION_TOL = 1.0e-2
const MIN_YEARS = 10
const FORCING_FILE = string(get(RC, "daily_forcing_file", ""))

data = prepare_ode_data(
    upstream_cost = Float64(get(RC, "upstream_cost", 0.05)),
    cedex_var_file = get(PARAMS["inputs"], "cedex_var_file", nothing),
    cedex_esc_uts_file = get(PARAMS["inputs"], "cedex_esc_uts_file", nothing),
    obstacles_file = get(PARAMS["inputs"], "obstacles_file", nothing),
    obstacle_mode = Symbol(get(PARAMS["obstacles"], "mode", "legacy")),
    obstacle_matching_tolerance = Float64(get(PARAMS["obstacles"], "matching_tolerance_m", 2000.0)),
    obstacle_passability = Float64(get(PARAMS["obstacles"], "upstream_passability", 0.1)),
    obstacle_downstream_passability = Float64(get(PARAMS["obstacles"], "downstream_passability", 0.5)),
    connectivity_method = SimulationParameters.connectivity_method(PARAMS),
    carrying_capacity_base_scaling = Float64(get(RC, "carrying_capacity_base_scaling", 1.0))
)
n_sites, n_species = data.params.n_sites, data.params.n_species
density_cols = [Symbol("$(sp)_DEN") for sp in data.species]
density_df = filter(row -> row.CODIGO in data.sites, data.density_df)
u0_flat = vec(max.(replace(Matrix(density_df[:, density_cols]), NaN => 0.0), 0.0))

daily = load_daily_forcing_any(FORCING_FILE)
matrix_of(df, sc) = is_wide_daily_forcing(df) ?
    wide_forcing_matrix(df, data.sites; scenario=sc) :
    daily_forcing_matrix(df, data.sites; scenario=sc)
bdates, btemps = matrix_of(daily, "historical")
schedule = baseline_climatology_schedule(btemps, bdates, 1; days_per_year=Float64(DAYS_PER_YEAR))

const STRESS = get(PARAMS, "temperature_stress", Dict{String,Any}())
k_heat = 0.0
if lowercase(string(get(STRESS, "enabled", false))) in ("1", "true", "yes")
    if lowercase(string(get(STRESS, "calibrate", false))) in ("1", "true", "yes")
        energies = Float64[]
        for s in 1:n_species, i in 1:n_sites
            push!(energies, exceedance_energy(@view(btemps[i, :]), data.params.thermal_upper_limits[s]))
        end
        k_heat = calibrate_heat_stress_rate(energies;
            max_annual_loss=Float64(get(STRESS, "max_annual_loss", 0.05)))
    else
        k_heat = Float64(get(STRESS, "k", 0.0))
    end
end
params = set_heat_stress_rate(data.params, k_heat)

println("="^70)
println("Verify corrected E4 burn-in: both, tol=$TOL, composition_tol=$COMPOSITION_TOL, min_years=$MIN_YEARS")
println("  active floor = $(Guadex.spin_up_composition_active_floor()), heat k=$k_heat, cap=$MAX_YEARS")
println("="^70)
t0 = time()
spin = spin_up(params; initial_state=u0_flat, schedule=schedule,
    days_per_year=Float64(DAYS_PER_YEAR), max_years=MAX_YEARS, tol=TOL,
    composition_tol=COMPOSITION_TOL, min_years=MIN_YEARS, criterion=:both,
    progress_every=25)
elapsed = time() - t0
println("\nRESULT: converged=$(spin.converged), years=$(spin.years) " *
        "($(round(elapsed, digits=1)) s, $(round(elapsed / max(spin.years, 1), digits=2)) s/yr)")
println("  basin change            = $(spin.last_basin_change)")
println("  composition (robust)    = $(spin.last_composition_change)")
println("  composition (all cells) = $(spin.last_composition_change_all_cells)")
println("  active cells            = $(spin.composition_active_cells)")
println("  criteria                = $(join(String.(spin.criteria), "+"))")
