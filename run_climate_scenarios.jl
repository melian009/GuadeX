using Pkg; Pkg.activate(".");
using DifferentialEquations
using DataFrames
using CSV
using SparseArrays
using LinearAlgebra
using JLD2
using Statistics
using Guadex

# =============================================================================
# --- Climate-driven simulations (all temperature projection scenarios) ---
#
# Runs the metacommunity model from the current year to 2045 with a year-by-year
# warming trajectory derived from the guadex_tw water-temperature projections,
# for every GCM x SSP combination (default 11 x 4 = 44 runs).  Each run writes
# the four reporting levels (sampling point, sub-basin, water body, whole basin)
# and the viz/ viewer files.
#
# Time step is daily (1 model time unit = 1 day); states are saved monthly so
# the outputs stay small while annual/seasonal metrics remain available.
#
# Settings come from the [run_climate_scenarios] section of parameters.toml.
# Individual settings can be overridden with the GUADEX_CLIMATE_* environment
# variables (see docs/climate_scenarios.md).
# =============================================================================

# --- Load run parameters (single source of truth: parameters.toml) ---
const _GUADEX_ROOT = isfile(joinpath(@__DIR__, "parameters.jl")) ? (@__DIR__) : dirname(@__DIR__)
include(joinpath(_GUADEX_ROOT, "parameters.jl"))
using .SimulationParameters
const GUADEX_PARAMS = SimulationParameters.load()
SimulationParameters.require_sections(GUADEX_PARAMS, "general", "inputs", "obstacles",
    "species", "subcatchments", "scenarios", "run_climate_scenarios")

const DAYS_PER_YEAR = Int(GUADEX_PARAMS["general"]["days_per_year"])

# Updated 2045 inputs, overridable per run.
const CEDEX_VAR_FILE = get(ENV, "GUADEX_CEDEX_VAR_FILE", GUADEX_PARAMS["inputs"]["cedex_var_file"])
const CEDEX_UTS_FILE = get(ENV, "GUADEX_CEDEX_UTS_FILE", GUADEX_PARAMS["inputs"]["cedex_esc_uts_file"])
const OBSTACLES_FILE = get(ENV, "GUADEX_OBSTACLES_FILE", GUADEX_PARAMS["inputs"]["obstacles_file"])
const OBSTACLE_MODE = Symbol(get(ENV, "GUADEX_OBSTACLE_MODE", GUADEX_PARAMS["obstacles"]["mode"]))
const OBSTACLE_MATCHING_TOLERANCE = parse(Float64, get(ENV, "GUADEX_OBSTACLE_TOLERANCE_M", string(GUADEX_PARAMS["obstacles"]["matching_tolerance_m"])))
const OBSTACLE_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["upstream_passability"])))
const OBSTACLE_DOWNSTREAM_PASSABILITY = parse(Float64, get(ENV, "GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY", string(GUADEX_PARAMS["obstacles"]["downstream_passability"])))

const NATIVE_SPECIES = String.(GUADEX_PARAMS["species"]["native"])
const INVASIVE_SPECIES = String.(GUADEX_PARAMS["species"]["invasive"])
const UPSTREAM_COST = Float64(get(GUADEX_PARAMS["run_climate_scenarios"], "upstream_cost", 0.05))

# --- Climate run settings ---
function _climate_setting(key, default)
    return get(GUADEX_PARAMS["run_climate_scenarios"], key, default)
end

const PROJECTION_FILE = get(ENV, "GUADEX_CLIMATE_PROJECTIONS_FILE",
    _climate_setting("temperature_projections_file", "guadex_tw/outputs/tables/water_temp_future_2045.csv"))
const START_YEAR = parse(Int, get(ENV, "GUADEX_CLIMATE_START_YEAR", string(_climate_setting("start_year", 2026))))
const END_YEAR = parse(Int, get(ENV, "GUADEX_CLIMATE_END_YEAR", string(_climate_setting("end_year", 2045))))
const SAVE_INTERVAL_DAYS = parse(Float64, get(ENV, "GUADEX_CLIMATE_SAVE_DAYS", string(_climate_setting("save_interval_days", 365.0 / 12.0))))
const ELEVATION_SCALING = lowercase(get(ENV, "GUADEX_CLIMATE_ELEVATION_SCALING",
    string(_climate_setting("elevation_scaling", false)))) in ("1", "true", "yes")
const PRESENCE_THRESHOLD = Float64(_climate_setting("presence_threshold", 0.1))
const MAX_RUNS = parse(Int, get(ENV, "GUADEX_CLIMATE_MAX_RUNS", "0"))  # 0 = no cap (testing aid)
const FORCE_RERUN = lowercase(get(ENV, "GUADEX_CLIMATE_FORCE", "0")) in ("1", "true", "yes")
const MAKE_FIGURES = lowercase(get(ENV, "GUADEX_CLIMATE_PLOT",
    string(_climate_setting("make_figures", true)))) in ("1", "true", "yes")

const CLIMATE_SCENARIOS = String.(get(GUADEX_PARAMS["run_climate_scenarios"], "scenarios",
    ["ssp126", "ssp245", "ssp370", "ssp585"]))
const DEFAULT_GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg",
    "IITM-ESM", "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0",
    "NorESM2-MM", "UKESM1-0-LL"]
const CLIMATE_GCMS = begin
    override = get(ENV, "GUADEX_CLIMATE_GCMS", "")
    isempty(override) ? String.(get(GUADEX_PARAMS["run_climate_scenarios"], "gcms", DEFAULT_GCMS)) :
        String.(split(override, ","))
end

const SIMULATION_YEARS = END_YEAR - START_YEAR + 1
const T_END = Float64(SIMULATION_YEARS * DAYS_PER_YEAR)
const YEAR_OFFSETS = collect(0:(SIMULATION_YEARS - 1))
const YEAR_LABELS = collect(START_YEAR:END_YEAR)

if get(ENV, "GUADEX_CONFIG_ONLY", "0") == "1"
    println("[parameters] climate-scenario configuration smoke test passed")
    println("  scenarios: $(join(CLIMATE_SCENARIOS, ", "))")
    println("  gcms: $(length(CLIMATE_GCMS))")
    println("  years: $START_YEAR-$END_YEAR (daily step, $(SAVE_INTERVAL_DAYS)-day saves)")
    println("  projections: $PROJECTION_FILE")
    exit(0)
end

# =============================================================================
# --- Main ---
# =============================================================================

println("="^70)
println("Climate-driven simulations: $START_YEAR-$END_YEAR")
println("  scenarios: $(join(CLIMATE_SCENARIOS, ", "))")
println("  GCMs:      $(length(CLIMATE_GCMS))")
println("  save:      every $(SAVE_INTERVAL_DAYS) days (monthly by default)")
println("  elevation scaling: $ELEVATION_SCALING")
println("="^70)

projections = load_temperature_projections(PROJECTION_FILE)
available_scenarios = Set(String.(projections.scenario))
available_gcms = Set(String.(projections.gcm))

println("\nLoading base data (once)...")
data_base = prepare_ode_data(
    upstream_cost = UPSTREAM_COST,
    cedex_var_file = CEDEX_VAR_FILE,
    cedex_esc_uts_file = CEDEX_UTS_FILE,
    obstacles_file = OBSTACLES_FILE,
    obstacle_mode = OBSTACLE_MODE,
    obstacle_matching_tolerance = OBSTACLE_MATCHING_TOLERANCE,
    obstacle_passability = OBSTACLE_PASSABILITY,
    obstacle_downstream_passability = OBSTACLE_DOWNSTREAM_PASSABILITY
)

n_sites = data_base.params.n_sites

# Initial conditions from the observed density matrix.
density_cols = [Symbol("$(sp)_DEN") for sp in data_base.species]
density_df_filtered = filter(row -> row.CODIGO in data_base.sites, data_base.density_df)
u0 = Matrix(density_df_filtered[:, density_cols])
replace!(u0, NaN => 0.0)
u0 = max.(u0, 0.0)
u0_flat = vec(u0)

saveat = unique(vcat(collect(0.0:SAVE_INTERVAL_DAYS:T_END), T_END))

function positivity_condition(u, t, integrator)
    return any(x -> x < 0, u)
end
function positivity_affect!(integrator)
    integrator.u .= max.(integrator.u, 0.0)
end
positivity_cb = DiscreteCallback(positivity_condition, positivity_affect!; save_positions=(false, true))

# Restrict path components to a safe character set so projection labels can
# never escape the results directory.
function safe_component(value)
    safe = replace(string(value), r"[^A-Za-z0-9._-]" => "_")
    safe = replace(safe, ".." => "__")
    return isempty(safe) ? "unnamed" : safe
end

function basin_row_from_csv(path, year)
    isfile(path) || return (NaN, NaN, NaN)
    df = CSV.read(path, DataFrame)
    ("year" in names(df)) || return (NaN, NaN, NaN)
    sub = df[df.year .== year, :]
    isempty(sub) && return (NaN, NaN, NaN)
    r = sub[1, :]
    return (r.mean_native_richness, r.mean_native_extinction_risk, r.mean_total_biomass)
end

base_output_dir = joinpath("results", "climate_scenarios")
mkpath(base_output_dir)
index_path = joinpath(base_output_dir, "runs_index.csv")
index_df = isfile(index_path) ? CSV.read(index_path, DataFrame) : DataFrame(
    scenario=String[], gcm=String[], start_year=Int[], end_year=Int[],
    warming_end_degc=Float64[], basin_native_richness=Float64[],
    basin_native_extinction_risk=Float64[], basin_total_biomass=Float64[],
    run_dir=String[])

run_number = 0
total_runs = length(CLIMATE_SCENARIOS) * length(CLIMATE_GCMS)

for scenario in CLIMATE_SCENARIOS
    scenario in available_scenarios || (@warn "scenario $scenario absent from projections; skipping"; continue)
    for gcm in CLIMATE_GCMS
        gcm in available_gcms || (@warn "GCM $gcm absent from projections; skipping"; continue)
        global run_number += 1
        if MAX_RUNS > 0 && run_number > MAX_RUNS
            println("[cap] GUADEX_CLIMATE_MAX_RUNS=$MAX_RUNS reached; stopping after $MAX_RUNS run(s)")
            break
        end

        run_dir = joinpath(base_output_dir, safe_component(scenario), safe_component(gcm))
        curve_years, curve = basin_warming_curve(projections, scenario, gcm, START_YEAR, END_YEAR)

        # Resume support: a completed run is skipped unless explicitly forced,
        # and its summary row is preserved in the incremental run index.
        if !FORCE_RERUN && isfile(joinpath(run_dir, "simulation_output.jld2"))
            println("\n[$run_number/$total_runs] $scenario / $gcm -> skipped (already exists)")
            if !any((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
                richness, risk, biomass = basin_row_from_csv(
                    joinpath(run_dir, "export", "levels", "level_basin.csv"), END_YEAR)
                push!(index_df, (scenario=scenario, gcm=gcm, start_year=START_YEAR,
                    end_year=END_YEAR, warming_end_degc=curve[end],
                    basin_native_richness=richness, basin_native_extinction_risk=risk,
                    basin_total_biomass=biomass, run_dir=run_dir))
                CSV.write(index_path, index_df)
            end
            continue
        end

        mkpath(run_dir)
        println("\n[$run_number/$total_runs] $scenario / $gcm -> $run_dir")

        warming = warming_matrix(n_sites, curve_years, curve;
            elevations=data_base.elevations, elevation_scaling=ELEVATION_SCALING)

        schedule = TemperatureSchedule(warming, Float64(DAYS_PER_YEAR))
        scheduled_params = ScheduledMetacommunityParams(data_base.params, schedule)

        prob = ODEProblem(metacommunity_ode_scheduled!, u0_flat, (0.0, T_END), scheduled_params)
        sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=saveat, callback=positivity_cb)
        println("  solved: $(length(sol.t)) snapshots, final basin delta = $(round(curve[end], digits=3)) degC")

        jldsave(joinpath(run_dir, "simulation_output.jld2");
            sol_t=sol.t, sol_u=sol.u,
            sites=data_base.sites, species=data_base.species,
            scenario=scenario, gcm=gcm,
            start_year=START_YEAR, end_year=END_YEAR,
            years=YEAR_LABELS, warming=warming,
            temperature_baseline=data_base.params.temperatures,
            habitat_suitability=data_base.params.habitat_suitability,
            save_interval_days=SAVE_INTERVAL_DAYS,
            elevation_scaling=ELEVATION_SCALING,
            upstream_cost=UPSTREAM_COST,
            cedex_var_file=CEDEX_VAR_FILE,
            cedex_uts_file=CEDEX_UTS_FILE,
            obstacles_file=OBSTACLES_FILE,
            obstacle_mode=string(OBSTACLE_MODE),
            obstacle_matched_count=count(data_base.obstacle_mapping_diagnostics.matched),
            obstacle_total_count=nrow(data_base.obstacle_mapping_diagnostics))

        export_result = export_run_outputs(joinpath(run_dir, "export");
            sol_t=sol.t, sol_u=sol.u,
            sites=data_base.sites, species=data_base.species, site_df=data_base.site_df,
            crosswalk_path=joinpath(_GUADEX_ROOT, "data", "site_waterbody_crosswalk.csv"),
            native_species=NATIVE_SPECIES, invasive_species=INVASIVE_SPECIES,
            days_per_year=DAYS_PER_YEAR, threshold=PRESENCE_THRESHOLD,
            report_year_offsets=YEAR_OFFSETS, report_year_labels=YEAR_LABELS,
            temperature_baseline=data_base.params.temperatures,
            warming=warming,
            dams=data_base.dams, distance_matrix=data_base.distance_matrix,
            habitat_suitability=data_base.params.habitat_suitability,
            upstream_cost=UPSTREAM_COST,
            require_crosswalk=true,
            run_metadata=Dict(
                "script" => "run_climate_scenarios.jl",
                "scenario" => scenario,
                "gcm" => gcm,
                "start_year" => START_YEAR,
                "end_year" => END_YEAR,
                "elevation_scaling" => ELEVATION_SCALING,
                "save_interval_days" => SAVE_INTERVAL_DAYS,
                "warming_end_degc" => curve[end],
                "projections_file" => PROJECTION_FILE,
                "obstacle_mode" => string(OBSTACLE_MODE)))

        basin = export_result.level_tables[:basin]
        final_basin = basin[basin.year .== END_YEAR, :]
        duplicate = findall((index_df.scenario .== scenario) .& (index_df.gcm .== gcm))
        isempty(duplicate) || deleteat!(index_df, duplicate)
        push!(index_df, (
            scenario=scenario, gcm=gcm,
            start_year=START_YEAR, end_year=END_YEAR,
            warming_end_degc=curve[end],
            basin_native_richness=isempty(final_basin) ? NaN : final_basin[1, :mean_native_richness],
            basin_native_extinction_risk=isempty(final_basin) ? NaN : final_basin[1, :mean_native_extinction_risk],
            basin_total_biomass=isempty(final_basin) ? NaN : final_basin[1, :mean_total_biomass],
            run_dir=run_dir))
        CSV.write(index_path, index_df)
    end
    if MAX_RUNS > 0 && run_number >= MAX_RUNS
        break
    end
end

if nrow(index_df) > 0
    println("\n" * "="^70)
    println("Climate scenarios complete: $(nrow(index_df)) run(s) written under '$base_output_dir'")
    println("Run index: $index_path")
    println("="^70)
else
    println("\nNo climate runs were produced.")
end

if MAKE_FIGURES
    println("\nGenerating climate figures...")
    try
        plot_climate_scenario_figures(base_output_dir;
            figures_dir=joinpath(base_output_dir, "figures"))
    catch err
        @warn "climate figure generation failed" exception=err
    end
else
    println("\nSkipping figures (GUADEX_CLIMATE_PLOT=0).")
end
