using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using Statistics
using Guadex

# =============================================================================
# --- Staged climate experiments (WP6) ---
#
# Implements the staged design from the modelling improvement plan, cheapest
# stage first and each gating the next:
#
#   E0 spin-up realism  : spun-up state vs observed (report, don't hide)
#   E1 seasonality      : baseline annual-mean vs baseline daily climatology
#   E2 heat stress      : heat-stress off vs on (Warm lowland sites)
#   E3 K sensitivity    : carrying capacity x {1, 3, 10}
#   E4 optimum sweep    : thermal optimum across the empirical range
#   E5 full ensemble    : delegated to run_climate_scenarios.jl
#
# Every case writes the standard four-level export under
# results/climate_experiments/<stage>/<case>/ and appends one row to
# results/climate_experiments/stage_index.csv.  Settings come from the
# [run_climate_experiments] section of parameters.toml.
# =============================================================================

const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "parameters.jl")) ? (@__DIR__) : dirname(@__DIR__)
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "scenarios", "temperature_stress", "run_climate_experiments")

const DAYS_PER_YEAR = Int(GUADEX_PARAMS["general"]["days_per_year"])
const CFG = GUADEX_PARAMS["run_climate_experiments"]

const START_YEAR = Int(get(CFG, "start_year", 2026))
const END_YEAR = Int(get(CFG, "end_year", 2045))
const SIMULATION_YEARS = END_YEAR - START_YEAR + 1
const T_END = Float64(SIMULATION_YEARS * DAYS_PER_YEAR)
const YEAR_OFFSETS = collect(0:(SIMULATION_YEARS - 1))
const YEAR_LABELS = collect(START_YEAR:END_YEAR)
const SAVE_INTERVAL_DAYS = 365.0 / 12.0
const PRESENCE_THRESHOLD = 0.1
const UPSTREAM_COST = 0.05

const STAGES = begin
    override = get(ENV, "GUADEX_EXPERIMENTS_STAGES", "")
    isempty(override) ? String.(get(CFG, "stages", ["E0", "E1", "E2", "E3", "E4"])) :
        String.(strip.(split(override, ",")))
end
const EXPERIMENT_SCENARIO = string(get(CFG, "scenario", "ssp585"))
const EXPERIMENT_GCMS = String.(get(CFG, "gcms", ["ACCESS-CM2"]))
const K_SCALING_GRID = Float64.(get(CFG, "k_scaling_grid", [1.0, 3.0, 10.0]))
const OPTIMUM_FRACTIONS = Float64.(get(CFG, "optimum_fractions", [0.0, 0.25, 0.5, 0.75, 1.0]))
const SPIN_UP_MAX_YEARS = Int(get(CFG, "spin_up_max_years", 50))
const SPIN_UP_TOL = Float64(get(CFG, "spin_up_tol", 1.0e-6))
const BASE_OUTPUT_DIR = string(get(CFG, "output_dir", "results/climate_experiments"))
const PROJECTION_FILE = string(get(CFG, "temperature_projections_file",
    "guadex_tw/outputs/tables/water_temp_future_2045.csv"))
const DAILY_FORCING_FILE = string(get(CFG, "daily_forcing_file",
    "guadex_tw/outputs/tables/water_temp_daily_guadex_sites.csv"))

const STRESS = get(GUADEX_PARAMS, "temperature_stress", Dict{String,Any}())
const HEAT_STRESS_K = parse(Float64, string(get(STRESS, "k", 0.0)))
const HEAT_STRESS_CALIBRATE = lowercase(string(get(STRESS, "calibrate", false))) in ("1", "true", "yes")
const HEAT_STRESS_MAX_LOSS = parse(Float64, string(get(STRESS, "max_annual_loss", 0.05)))

const NATIVE_SPECIES = String.(GUADEX_PARAMS["species"]["native"])
const INVASIVE_SPECIES = String.(GUADEX_PARAMS["species"]["invasive"])
const MIGRATORY_SPECIES = String.(get(GUADEX_PARAMS["species"], "migratory", String[]))

println("="^70)
println("Staged climate experiments (WP6): $(join(STAGES, ", "))")
println("  scenario: $EXPERIMENT_SCENARIO; GCMs: $(join(EXPERIMENT_GCMS, ", "))")
println("  horizon:  $START_YEAR-$END_YEAR")
println("="^70)

# --- Base data (once) ---
println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = UPSTREAM_COST,
    cedex_var_file = get(GUADEX_PARAMS["inputs"], "cedex_var_file", nothing),
    cedex_esc_uts_file = get(GUADEX_PARAMS["inputs"], "cedex_esc_uts_file", nothing),
    obstacles_file = get(GUADEX_PARAMS["inputs"], "obstacles_file", nothing),
    obstacle_mode = Symbol(get(GUADEX_PARAMS["obstacles"], "mode", "legacy")),
    obstacle_matching_tolerance = Float64(get(GUADEX_PARAMS["obstacles"], "matching_tolerance_m", 2000.0)),
    obstacle_passability = Float64(get(GUADEX_PARAMS["obstacles"], "upstream_passability", 0.1)),
    obstacle_downstream_passability = Float64(get(GUADEX_PARAMS["obstacles"], "downstream_passability", 0.5))
)
n_sites = data_base.params.n_sites
n_species = data_base.params.n_species

density_cols = [Symbol("$(sp)_DEN") for sp in data_base.species]
density_df_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
u0_obs = max.(replace(Matrix(density_df_filtered[:, density_cols]), NaN => 0.0), 0.0)
u0_flat = vec(u0_obs)

saveat = unique(vcat(collect(0.0:SAVE_INTERVAL_DAYS:T_END), T_END))
positivity_cb = DiscreteCallback(
    (u, t, integrator) -> any(x -> x < 0, u),
    integrator -> (integrator.u .= max.(integrator.u, 0.0));
    save_positions=(false, true))

mkpath(BASE_OUTPUT_DIR)
index_path = joinpath(BASE_OUTPUT_DIR, "stage_index.csv")
index_df = isfile(index_path) ? CSV.read(index_path, DataFrame) : DataFrame(
    stage=String[], case=String[], scenario=String[], gcm=String[],
    forcing_mode=String[], heat_stress_k=Float64[], k_scaling=Float64[],
    optimum_fraction=Float64[], end_basin_biomass=Float64[],
    end_native_richness=Float64[], end_native_extinction_risk=Float64[],
    run_dir=String[])

function record_case(; stage, case_name, scenario, gcm, forcing_mode, k, k_scaling,
        fraction, run_dir)
    global index_df
    basin_path = joinpath(run_dir, "export", "levels", "level_basin.csv")
    biomass = richness = risk = NaN
    if isfile(basin_path)
        basin = CSV.read(basin_path, DataFrame)
        final = basin[basin.year .== END_YEAR, :]
        if !isempty(final)
            biomass = final[1, :mean_total_biomass]
            richness = final[1, :mean_native_richness]
            risk = final[1, :mean_native_extinction_risk]
        end
    end
    filter!(row -> !(row.stage == stage && row.case == case_name), index_df)
    push!(index_df, (stage=stage, case=case_name, scenario=scenario, gcm=gcm,
        forcing_mode=forcing_mode, heat_stress_k=k, k_scaling=k_scaling,
        optimum_fraction=fraction, end_basin_biomass=biomass,
        end_native_richness=richness, end_native_extinction_risk=risk, run_dir=run_dir))
    CSV.write(index_path, index_df)
end

function solve_case(; stage, case_name, params, schedule, u0, warming,
        forcing_temps=nothing, forcing_dates=nothing, baseline_species_density,
        metadata=Dict())
    run_dir = joinpath(BASE_OUTPUT_DIR, stage, case_name)
    mkpath(run_dir)
    scheduled = ScheduledMetacommunityParams(params, schedule)
    prob = ODEProblem(metacommunity_ode_scheduled!, u0, (0.0, T_END), scheduled)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=saveat, callback=positivity_cb)
    export_run_outputs(joinpath(run_dir, "export");
        sol_t=sol.t, sol_u=sol.u, sites=data_base.sites, species=data_base.species,
        site_df=data_base.site_df,
        crosswalk_path=joinpath(_GUADEX_ROOT, "data", "site_waterbody_crosswalk.csv"),
        native_species=NATIVE_SPECIES, invasive_species=INVASIVE_SPECIES,
        migratory_species=MIGRATORY_SPECIES,
        days_per_year=DAYS_PER_YEAR, threshold=PRESENCE_THRESHOLD,
        report_year_offsets=YEAR_OFFSETS, report_year_labels=YEAR_LABELS,
        temperature_baseline=params.temperatures, warming=warming,
        dams=data_base.dams, distance_matrix=data_base.distance_matrix,
        habitat_suitability=params.habitat_suitability, upstream_cost=UPSTREAM_COST,
        require_crosswalk=true,
        baseline_species_density=baseline_species_density,
        species_upper_limits=params.thermal_upper_limits,
        daily_forcing=forcing_temps === nothing ? nothing :
            (temps=forcing_temps, dates=forcing_dates),
        run_metadata=merge(Dict("script" => "run_climate_experiments.jl",
            "stage" => stage, "case" => case_name), metadata))
    return run_dir
end

# ---------------------------------------------------------------------------
# Setup shared across stages
# ---------------------------------------------------------------------------
baseline_params = data_base.params

"""Spin up `params` under baseline forcing and return (result, n_sites × n_species state)."""
function spin_for(params)
    sched = baseline_temps === nothing ? nothing :
        baseline_climatology_schedule(baseline_temps, baseline_dates, 1;
            days_per_year=Float64(DAYS_PER_YEAR))
    s = spin_up(params; initial_state=u0_flat, schedule=sched,
        days_per_year=Float64(DAYS_PER_YEAR),
        max_years=SPIN_UP_MAX_YEARS, tol=SPIN_UP_TOL)
    return s, reshape(s.state, n_sites, n_species)
end

# Optional daily baseline climatology (WP1 product), loaded before the spin-up
# so the equilibrium is computed under the same seasonal forcing as the runs.
baseline_dates = nothing
baseline_temps = nothing
if isfile(DAILY_FORCING_FILE)
    daily = load_daily_forcing_any(DAILY_FORCING_FILE)
    if "historical" in Set(String.(daily.scenario))
        baseline_dates, baseline_temps = is_wide_daily_forcing(daily) ?
            wide_forcing_matrix(daily, data_base.sites; scenario="historical") :
            daily_forcing_matrix(daily, data_base.sites; scenario="historical")
        println("daily baseline available: $(length(baseline_dates)) days")
    end
end

spin, baseline_species_density = spin_for(baseline_params)
u0_spun = spin.state

# Heat-stress slope: calibrate from the baseline if requested and possible.
heat_k = HEAT_STRESS_K
if HEAT_STRESS_CALIBRATE && baseline_temps !== nothing
    energies = Float64[]
    for s in 1:n_species, i in 1:n_sites
        push!(energies, exceedance_energy(@view(baseline_temps[i, :]),
            baseline_params.thermal_upper_limits[s]))
    end
    heat_k = calibrate_heat_stress_rate(energies; max_annual_loss=HEAT_STRESS_MAX_LOSS)
    println("WP3 heat-stress calibration: k=$heat_k")
end

projections = isfile(PROJECTION_FILE) ? load_temperature_projections(PROJECTION_FILE) : nothing

function experiment_schedule(scenario, gcm; zero=false)
    if zero || projections === nothing
        return annual_mean_schedule(zeros(n_sites, SIMULATION_YEARS);
            days_per_year=Float64(DAYS_PER_YEAR)), zeros(n_sites, SIMULATION_YEARS)
    end
    years, curve = basin_warming_curve(projections, scenario, gcm, START_YEAR, END_YEAR)
    warming = warming_matrix(n_sites, years, curve; elevations=data_base.elevations)
    return annual_mean_schedule(warming; days_per_year=Float64(DAYS_PER_YEAR)), warming
end

const BASELINE_START = 1986
const BASELINE_END = 2005

"""
    experiment_daily_schedule(scenario, gcm)

Per-site daily schedule from the WP1 forcing (per-GCM file when present, else the
ensemble-median file).  Used by E2-E4: heat stress acts only when daily
temperatures exceed the empirical upper limits, which annual means never reach.
Returns `(schedule, annual_mean_warming, forcing_temps, forcing_dates)`.
"""
function experiment_daily_schedule(scenario, gcm)
    path = DAILY_FORCING_FILE
    if isfile(DAILY_FORCING_FILE)
        root, ext = splitext(DAILY_FORCING_FILE)
        candidate = string(root, "_", scenario, "_", gcm, ext)
        isfile(candidate) && (path = candidate)
    end
    daily = load_daily_forcing_any(path)
    m(df, sc) = is_wide_daily_forcing(df) ?
        wide_forcing_matrix(df, data_base.sites; scenario=sc) :
        daily_forcing_matrix(df, data_base.sites; scenario=sc)
    fdates, ftemps = m(daily, scenario)
    bdates, btemps = m(daily, "historical")
    sched, _ = daily_temperature_schedule(ftemps; dates=fdates,
        baseline_temps=btemps, baseline_dates=bdates,
        baseline_start=BASELINE_START, baseline_end=BASELINE_END)
    warming = annual_mean_deltas(sched; days_per_year=Float64(DAYS_PER_YEAR))
    return sched, warming, ftemps, fdates
end

# ---------------------------------------------------------------------------
# E0 - spin-up realism
# ---------------------------------------------------------------------------
if "E0" in STAGES
    observed_total = site_totals(vec(u0_obs), n_sites, n_species)
    spun_total = spin.site_biomass
    observed_native_richness = [count(>(PRESENCE_THRESHOLD), @view u0_obs[i,
        species_indices(data_base.species, NATIVE_SPECIES)]) for i in 1:n_sites]
    spun_matrix = reshape(u0_spun, n_sites, n_species)
    spun_native_richness = [count(>(PRESENCE_THRESHOLD), @view spun_matrix[i,
        species_indices(data_base.species, NATIVE_SPECIES)]) for i in 1:n_sites]
    ratio = fill(NaN, n_sites)
    positive = observed_total .> 0
    ratio[positive] = spun_total[positive] ./ observed_total[positive]
    summary = DataFrame(
        site=data_base.sites,
        observed_total_biomass=observed_total,
        spun_up_total_biomass=spun_total,
        total_biomass_ratio=ratio,
        observed_native_richness=observed_native_richness,
        spun_up_native_richness=spun_native_richness)
    mkpath(joinpath(BASE_OUTPUT_DIR, "E0"))
    CSV.write(joinpath(BASE_OUTPUT_DIR, "E0", "spinup_summary.csv"), summary)
    spin_meta = Dict("years" => spin.years, "converged" => spin.converged,
        "last_relative_change" => spin.last_relative_change)
    open(joinpath(BASE_OUTPUT_DIR, "E0", "spinup_metadata.json"), "w") do io
        print(io, spin_meta)
    end
    finite_ratio = filter(isfinite, ratio)
    println("E0: spin-up $(spin.years) yr, converged=$(spin.converged), " *
            "change=$(spin.last_relative_change); " *
            "median spun/observed biomass = $(isempty(finite_ratio) ? NaN : round(median(finite_ratio), digits=3)) " *
            "(≥ 1 means the spun-up state is above the observed snapshot)")
end

# ---------------------------------------------------------------------------
# E1 - seasonality (baseline climatology, no trend)
# ---------------------------------------------------------------------------
if "E1" in STAGES
    p_nostress = set_heat_stress_rate(baseline_params, 0.0)
    schedule_annual, warming_zero = experiment_schedule(EXPERIMENT_SCENARIO,
        first(EXPERIMENT_GCMS); zero=true)
    dir = solve_case(; stage="E1", case_name="annual_mean", params=p_nostress,
        schedule=schedule_annual, u0=u0_spun, warming=warming_zero,
        baseline_species_density=baseline_species_density,
        metadata=Dict("forcing_mode" => "annual_mean", "heat_stress_k" => 0.0,
            "k_scaling" => 1.0, "optimum_fraction" => 0.5))
    record_case(stage="E1", case_name="annual_mean", scenario="baseline",
        gcm="none", forcing_mode="annual_mean", k=0.0, k_scaling=1.0,
        fraction=0.5, run_dir=dir)

    if baseline_temps !== nothing
        clim_schedule = baseline_climatology_schedule(baseline_temps, baseline_dates,
            SIMULATION_YEARS; days_per_year=Float64(DAYS_PER_YEAR))
        warming_clim = annual_mean_deltas(clim_schedule; days_per_year=Float64(DAYS_PER_YEAR))
        dir = solve_case(; stage="E1", case_name="daily_climatology", params=p_nostress,
            schedule=clim_schedule, u0=u0_spun, warming=warming_clim,
            forcing_temps=baseline_temps, forcing_dates=baseline_dates,
            baseline_species_density=baseline_species_density,
            metadata=Dict("forcing_mode" => "daily", "heat_stress_k" => 0.0,
                "k_scaling" => 1.0, "optimum_fraction" => 0.5))
        record_case(stage="E1", case_name="daily_climatology", scenario="baseline",
            gcm="none", forcing_mode="daily", k=0.0, k_scaling=1.0,
            fraction=0.5, run_dir=dir)
    else
        @warn "E1 daily case skipped: no daily baseline series ($DAILY_FORCING_FILE)"
    end
end

# ---------------------------------------------------------------------------
# E2 - heat stress off vs on
# ---------------------------------------------------------------------------
if "E2" in STAGES
    gcm = first(EXPERIMENT_GCMS)
    schedule, warming, ftemps, fdates = experiment_daily_schedule(EXPERIMENT_SCENARIO, gcm)
    for (name, k) in (("heat_off", 0.0), ("heat_on", heat_k))
        params = set_heat_stress_rate(baseline_params, k)
        _, baseline_density = spin_for(params)
        dir = solve_case(; stage="E2", case_name=name, params=params,
            schedule=schedule, u0=vec(baseline_density), warming=warming,
            forcing_temps=ftemps, forcing_dates=fdates,
            baseline_species_density=baseline_density,
            metadata=Dict("forcing_mode" => "daily", "heat_stress_k" => k,
                "k_scaling" => 1.0, "optimum_fraction" => 0.5))
        record_case(stage="E2", case_name=name, scenario=EXPERIMENT_SCENARIO,
            gcm=gcm, forcing_mode="daily", k=k, k_scaling=1.0,
            fraction=0.5, run_dir=dir)
    end
end

# ---------------------------------------------------------------------------
# E3 - carrying-capacity sensitivity
# ---------------------------------------------------------------------------
if "E3" in STAGES
    gcm = first(EXPERIMENT_GCMS)
    schedule, warming, ftemps, fdates = experiment_daily_schedule(EXPERIMENT_SCENARIO, gcm)
    for scaling in K_SCALING_GRID
        name = "K_$(scaling)x"
        params = scale_carrying_capacity(set_heat_stress_rate(baseline_params, heat_k), scaling)
        # Re-spin under the scaled carrying capacity so the case starts from its
        # own equilibrium rather than importing the transient of the 1x K state.
        _, baseline_density = spin_for(params)
        dir = solve_case(; stage="E3", case_name=name, params=params,
            schedule=schedule, u0=vec(baseline_density), warming=warming,
            forcing_temps=ftemps, forcing_dates=fdates,
            baseline_species_density=baseline_density,
            metadata=Dict("forcing_mode" => "daily", "heat_stress_k" => heat_k,
                "k_scaling" => scaling, "optimum_fraction" => 0.5))
        record_case(stage="E3", case_name=name, scenario=EXPERIMENT_SCENARIO,
            gcm=gcm, forcing_mode="daily", k=heat_k, k_scaling=scaling,
            fraction=0.5, run_dir=dir)
    end
end

# ---------------------------------------------------------------------------
# E4 - thermal-optimum sweep
# ---------------------------------------------------------------------------
if "E4" in STAGES
    gcm = first(EXPERIMENT_GCMS)
    schedule, warming, ftemps, fdates = experiment_daily_schedule(EXPERIMENT_SCENARIO, gcm)
    for fraction in OPTIMUM_FRACTIONS
        name = "optimum_$(fraction)"
        optima = optimum_sweep_optima(data_base.species_chars_df, data_base.species;
            fraction=fraction)
        params = set_heat_stress_rate(set_thermal_optima(baseline_params, optima), heat_k)
        _, baseline_density = spin_for(params)
        dir = solve_case(; stage="E4", case_name=name, params=params,
            schedule=schedule, u0=vec(baseline_density), warming=warming,
            forcing_temps=ftemps, forcing_dates=fdates,
            baseline_species_density=baseline_density,
            metadata=Dict("forcing_mode" => "daily", "heat_stress_k" => heat_k,
                "k_scaling" => 1.0, "optimum_fraction" => fraction))
        record_case(stage="E4", case_name=name, scenario=EXPERIMENT_SCENARIO,
            gcm=gcm, forcing_mode="daily", k=heat_k, k_scaling=1.0,
            fraction=fraction, run_dir=dir)
    end
end

# ---------------------------------------------------------------------------
# E5 - full ensemble (delegated)
# ---------------------------------------------------------------------------
if "E5" in STAGES
    println("\nE5 (full GCM x SSP ensemble) is delegated to run_climate_scenarios.jl:")
    println("  GUADEX_CLIMATE_FORCING_MODE=daily GUADEX_CLIMATE_SPIN_UP=1 " *
            "GUADEX_CLIMATE_HEAT_STRESS=1 julia --project=. run_climate_scenarios.jl")
end

if nrow(index_df) > 0
    println("\n" * "="^70)
    println("Staged experiments complete: $(nrow(index_df)) case(s)")
    println("Index: $index_path")
    println("="^70)
end
