# =============================================================================
# Shared climate-consistent core for the obstacle / upstream-cost sensitivity
# sweeps (`run_sensitivity_report.jl`, `run_alt_interactions.jl`).
#
# The sensitivity runs are made directly comparable to the climate-scenario
# ensemble (`run_climate_scenarios.jl`): per-site daily water temperature, a
# seasonal baseline burn-in, calibrated heat stress, the WP0 species
# classification and the effective K = `carrying_capacity_base_scaling`
# convention.  The only new axes are the passability scenario and the upstream
# dispersal cost.
#
# This file is `include`d after `Guadex` and `parameters.jl`; it defines plain
# functions in the including module.
# =============================================================================

using Statistics

const SENSITIVITY_BASELINE_START = 1986
const SENSITIVITY_BASELINE_END = 2005

# ---------------------------------------------------------------------------
# Scenario / passability helpers
# ---------------------------------------------------------------------------

"""
    build_dam_passability_vector(site_df, sites; passability_per_subcatchment)

Per-site multiplier applied to the obstacle/dam passability matrix, taken from
the per-subcatchment passability scenario (sites in unlisted subcatchments get
1.0).  Mirrors the legacy sensitivity scripts.
"""
function build_dam_passability_vector(site_df, sites;
        passability_per_subcatchment::Dict{Float64,Float64}=Dict{Float64,Float64}())
    n_sites = length(sites)
    passability_vector = ones(n_sites)
    site_to_sc = Dict(row.CODIGO => row.CODIGO_S for row in eachrow(site_df))
    for (i, site) in enumerate(sites)
        if haskey(site_to_sc, site) && haskey(passability_per_subcatchment, site_to_sc[site])
            passability_vector[i] = passability_per_subcatchment[site_to_sc[site]]
        end
    end
    return passability_vector
end

"""
    classify_species_indices(all_species, target_group)

Indices of `target_group` species inside `all_species` (order preserved, codes
absent from the model are dropped).
"""
function classify_species_indices(all_species, target_group)
    return [findfirst(==(sp), all_species) for sp in target_group if sp in all_species]
end

"""
    modified_dams_for_scenario(base_dams, passability_vector)

Element-wise scenario multiplier applied to every restricted link of the base
obstacle/dam matrix (links already at 1.0 are unrestricted and stay 1.0).  A
factor above 1.0 improves passage, capped at full passability.
"""
function modified_dams_for_scenario(base_dams, passability_vector)
    modified = copy(base_dams)
    n = size(modified, 1)
    for j in 1:n, i in 1:n
        if i != j && modified[i, j] < 1.0
            modified[i, j] = min(1.0, modified[i, j] * passability_vector[j])
        end
    end
    return modified
end

"""
    per_gcm_forcing_path(base, scenario, gcm)

Path of the per-GCM wide daily forcing file written by
`guadex_tw/scripts/13_project_guadex_sites.py`:
`<base>_<scenario>_<gcm>.csv`.
"""
function per_gcm_forcing_path(base::AbstractString, scenario::AbstractString, gcm::AbstractString)
    root, ext = splitext(base)
    return string(root, "_", safe_component(scenario), "_", safe_component(gcm), ext)
end

"""
    safe_component(value)

Filesystem-safe scenario/GCM label.
"""
function safe_component(value)
    safe = replace(string(value), r"[^A-Za-z0-9._-]" => "_")
    safe = replace(safe, ".." => "__")
    return isempty(safe) ? "unnamed" : safe
end

# ---------------------------------------------------------------------------
# Shared settings
# ---------------------------------------------------------------------------

function _bool_setting(value)
    return lowercase(string(value)) in ("1", "true", "yes", "on")
end

"""
    load_sensitivity_settings(GUADEX_PARAMS)

Collect the climate-consistent settings shared by the sensitivity sweeps.  The
`[obstacle_sensitivity]` section is the explicit source of truth (falling back
to `[run_climate_scenarios]` and `[temperature_stress]`); its shipped values
reproduce the `climate_scenarios_k1x_burnin` configuration: per-GCM daily
forcing, a seasonal baseline burn-in, calibrated heat stress and K = 1x
observed.  Shorten the horizon for smoke tests with
`GUADEX_SENSITIVITY_END_YEAR`.
"""
function load_sensitivity_settings(GUADEX_PARAMS)
    climate = GUADEX_PARAMS["run_climate_scenarios"]
    os = get(GUADEX_PARAMS, "obstacle_sensitivity", Dict{String,Any}())
    stress = get(GUADEX_PARAMS, "temperature_stress", Dict{String,Any}())
    shared(key, default) = get(os, key, get(climate, key, default))

    days_per_year = Int(GUADEX_PARAMS["general"]["days_per_year"])
    start_year = Int(shared("start_year", 2026))
    end_year = parse(Int, get(ENV, "GUADEX_SENSITIVITY_END_YEAR",
        string(shared("end_year", 2045))))
    end_year >= start_year || error("end_year ($end_year) must be >= start_year ($start_year)")
    simulation_years = end_year - start_year + 1

    smoke = _bool_setting(get(ENV, "GUADEX_SENSITIVITY_SMOKE", "0"))

    return (
        days_per_year = days_per_year,
        start_year = start_year,
        end_year = end_year,
        simulation_years = simulation_years,
        year_offsets = collect(0:(simulation_years - 1)),
        year_labels = collect(start_year:end_year),
        save_interval_days = Float64(shared("save_interval_days", 30.4167)),
        presence_threshold = Float64(shared("presence_threshold", 0.1)),
        daily_forcing_file = get(ENV, "GUADEX_SENSITIVITY_DAILY_FILE",
            string(shared("daily_forcing_file",
                "guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide.csv"))),
        baseline_period_start = Int(shared("baseline_period_start",
            SENSITIVITY_BASELINE_START)),
        baseline_period_end = Int(shared("baseline_period_end",
            SENSITIVITY_BASELINE_END)),
        carrying_capacity_base_scaling = parse(Float64, get(ENV, "GUADEX_SENSITIVITY_K_BASE",
            string(shared("carrying_capacity_base_scaling", 1.0)))),
        carrying_capacity_scaling = Float64(shared("carrying_capacity_scaling", 1.0)),
        thermal_optima_fraction = Float64(shared("thermal_optima_fraction", 0.5)),
        spin_up = _bool_setting(get(ENV, "GUADEX_SENSITIVITY_SPIN_UP",
            string(shared("spin_up", true)))),
        spin_up_max_years = parse(Int, get(ENV, "GUADEX_SENSITIVITY_SPINUP_MAX_YEARS",
            string(smoke ? 3 : shared("spin_up_max_years", 600)))),
        spin_up_tol = parse(Float64, get(ENV, "GUADEX_SENSITIVITY_SPINUP_TOL",
            string(smoke ? 1.0 : shared("spin_up_tol", 5.0e-5)))),
        spin_up_criterion = Symbol(lowercase(get(ENV, "GUADEX_SENSITIVITY_SPINUP_CRITERION",
            string(shared("spin_up_criterion", "basin"))))),
        spin_up_progress_every = Int(shared("spin_up_progress_every", 50)),
        heat_stress_enabled = _bool_setting(get(ENV, "GUADEX_SENSITIVITY_HEAT_STRESS",
            string(get(os, "heat_stress_enabled",
                get(stress, "enabled", true))))),
        heat_stress_k = Float64(get(stress, "k", 0.0)),
        heat_stress_calibrate = _bool_setting(get(ENV, "GUADEX_SENSITIVITY_HEAT_STRESS_CALIBRATE",
            string(get(os, "heat_stress_calibrate",
                get(stress, "calibrate", true))))),
        heat_stress_max_loss = Float64(get(os, "heat_stress_max_loss",
            get(stress, "max_annual_loss", 0.05))),
        quasi_extinction_q = Float64(shared("quasi_extinction_q", 0.1)),
        quasi_extinction_persistence = Int(shared("quasi_extinction_persistence", 3)),
        upstream_cost_default = Float64(shared("upstream_cost", 0.05)),
        output_dir = get(ENV, "GUADEX_SENSITIVITY_OUTPUT_DIR",
            string(shared("output_dir", "results/sensitivity_obstacles"))),
        smoke = smoke,
    )
end

"""
    parse_climate_models(spec)

Parse a list of `{ scenario = "...", gcm = "..." }` entries (from
`parameters.toml` or the `GUADEX_SENSITIVITY_CLIMATE_MODELS` override
`"ssp126:IITM-ESM,ssp245:UKESM1-0-LL"`) into named tuples.  Errors loudly on a
malformed entry so a typo cannot silently drop a climate model.
"""
function parse_climate_models(spec)
    models = NamedTuple{(:scenario, :gcm),Tuple{String,String}}[]
    if spec isa AbstractString
        entries = split(spec, ",")
        for entry in entries
            isempty(strip(entry)) && continue
            parts = split(entry, ":")
            length(parts) == 2 || error("climate model override entry '$entry' must be scenario:gcm")
            push!(models, (scenario=strip(parts[1]), gcm=strip(parts[2])))
        end
    else
        for entry in spec
            haskey(entry, "scenario") && haskey(entry, "gcm") ||
                error("each climate model needs a 'scenario' and a 'gcm': $entry")
            push!(models, (scenario=String(entry["scenario"]), gcm=String(entry["gcm"])))
        end
    end
    isempty(models) && error("at least one climate model is required")
    return models
end

function climate_model_tag(model)
    return "$(safe_component(model.scenario))__$(safe_component(model.gcm))"
end

# ---------------------------------------------------------------------------
# Daily forcing
# ---------------------------------------------------------------------------

"""
    load_baseline_forcing(path, sites; scenario="historical")

Per-site daily baseline (historical) series from a wide or long forcing file.
"""
function load_baseline_forcing(path::AbstractString, sites; scenario::AbstractString="historical")
    isfile(path) || error("daily temperature forcing file not found: $path")
    df = load_daily_forcing_any(path)
    scenario in Set(string.(df.scenario)) ||
        error("forcing file $path has no '$scenario' scenario")
    if is_wide_daily_forcing(df)
        return wide_forcing_matrix(df, sites; scenario=scenario)
    end
    return daily_forcing_matrix(df, sites; scenario=scenario)
end

"""
    load_climate_model_forcing(forcing_file, sites, scenario, gcm; kwargs...)

Load the per-GCM daily series and build the daily [`TemperatureSchedule`](@ref)
plus the annual-mean warming matrix and the exposure forcing.  The anomaly is
taken against the same file's `historical` scenario (1986-2005 by default), so
the run reproduces `run_climate_scenarios.jl` in per-GCM daily mode.
"""
function load_climate_model_forcing(forcing_file::AbstractString, sites,
        scenario::AbstractString, gcm::AbstractString;
        baseline_start::Int=SENSITIVITY_BASELINE_START,
        baseline_end::Int=SENSITIVITY_BASELINE_END,
        days_per_year::Real=365.0,
        year_labels::AbstractVector=Int[])
    path = per_gcm_forcing_path(forcing_file, scenario, gcm)
    isfile(path) || error("per-GCM daily forcing not found: $path")
    src = load_daily_forcing_any(path)
    wide = is_wide_daily_forcing(src)
    pivot(df, sc) = wide ? wide_forcing_matrix(df, sites; scenario=sc) :
        daily_forcing_matrix(df, sites; scenario=sc)
    fdates, ftemps = pivot(src, scenario)
    bdates, btemps = pivot(src, "historical")
    schedule, _ = daily_temperature_schedule(ftemps; dates=fdates,
        baseline_temps=btemps, baseline_dates=bdates,
        baseline_start=baseline_start, baseline_end=baseline_end)
    years_out, warming_out = annual_mean_deltas_by_year(schedule, fdates)
    columns = [findfirst(==(y), years_out) for y in year_labels]
    any(isnothing, columns) &&
        error("daily forcing ($scenario/$gcm) does not cover every year in " *
              "$(first(year_labels))-$(last(year_labels))")
    warming = warming_out[:, [something(c) for c in columns]]
    return (
        schedule = schedule,
        warming = warming,
        forcing_temps = ftemps,
        forcing_dates = fdates,
        baseline_means = [mean(@view btemps[i, :]) for i in 1:size(btemps, 1)],
    )
end

# ---------------------------------------------------------------------------
# Heat stress
# ---------------------------------------------------------------------------

"""
    resolve_heat_stress_rate(settings, thermal_upper_limits, baseline_temps)

Return the shared heat-stress slope `k`: calibrated from the baseline
exceedance energy when requested, the configured value otherwise, or `0.0` when
heat stress is disabled.
"""
function resolve_heat_stress_rate(settings, thermal_upper_limits, baseline_temps)
    settings.heat_stress_enabled || return 0.0
    if settings.heat_stress_calibrate
        baseline_temps === nothing &&
            error("heat-stress calibration needs a daily baseline series")
        energies = Float64[]
        for s in eachindex(thermal_upper_limits), i in 1:size(baseline_temps, 1)
            push!(energies, exceedance_energy(@view(baseline_temps[i, :]),
                thermal_upper_limits[s]))
        end
        return calibrate_heat_stress_rate(energies; max_annual_loss=settings.heat_stress_max_loss)
    end
    return settings.heat_stress_k
end

# ---------------------------------------------------------------------------
# Burn-in
# ---------------------------------------------------------------------------

"""
    _stable_digest(s)

Deterministic 16-hex-character FNV-1a digest, used to keep cache filenames short
enough for Windows path limits.  Unlike `Base.hash`, it is stable across
sessions and Julia versions.
"""
function _stable_digest(s::AbstractString)
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ UInt64(b)) * 0x00000100000001b3
    end
    return string(h, base=16, pad=16)
end

function _spinup_cache_key(settings; tag="default", extra::AbstractString="",
        initial_state=nothing, baseline_temps=nothing)
    initial_fp = initial_state === nothing ? "na" :
        "n$(length(initial_state))_s$(round(sum(initial_state), digits=6))"
    baseline_fp = baseline_temps === nothing ? "na" :
        "d$(size(baseline_temps, 1))x$(size(baseline_temps, 2))_s$(round(sum(baseline_temps), digits=4))"
    forcing_path = String(settings.daily_forcing_file)
    forcing_fp = isfile(forcing_path) ?
        "sz$(stat(forcing_path).size)_mt$(round(Int, stat(forcing_path).mtime))" : "missing"
    pieces = [
        tag,
        "kbase=$(settings.carrying_capacity_base_scaling)",
        "kcal=$(settings.carrying_capacity_scaling)",
        "opt=$(settings.thermal_optima_fraction)",
        "heaton=$(settings.heat_stress_enabled)",
        "heatk=$(settings.heat_stress_k)",
        "heatcal=$(settings.heat_stress_calibrate)",
        "maxloss=$(settings.heat_stress_max_loss)",
        "years=$(settings.spin_up_max_years)",
        "tol=$(settings.spin_up_tol)",
        "crit=$(settings.spin_up_criterion)",
        "base=$(settings.baseline_period_start)-$(settings.baseline_period_end)",
        "file=$(basename(forcing_path))",
        "forc=$forcing_fp",
        "uc=$(settings.upstream_cost_default)",
        "extra=$extra",
        "init=$initial_fp",
        "baseline=$baseline_fp",
    ]
    # Descriptive prefix for humans, plus a stable digest of the full key so the
    # filename stays short (Windows path limits) without losing any settings.
    prefix = safe_component(tag)
    return "$(first(prefix, min(length(prefix), 24)))_$(_stable_digest(join(pieces, "|")))"
end

"""
    get_or_compute_spinup(settings, params, initial_state, baseline_temps, baseline_dates; tag, extra)

Integrate to the baseline equilibrium once and cache the state under
`results/sensitivity_spinup_cache/` so repeated scripts (and a resumed sweep)
reuse it.  The cache key covers the burn-in controls, the effective upstream
cost, the caller-supplied `extra` fingerprint (obstacle configuration and, for
the alt sweep, the interaction matrix + sigma), the forcing file size/mtime and
the initial state, so a changed run configuration cannot silently reuse a stale
equilibrium.  Disable reuse with `GUADEX_SPINUP_REUSE=0`, or clear the cache with
`GUADEX_SPINUP_FORCE=1`.
"""
function get_or_compute_spinup(settings, params::MetacommunityParams, initial_state,
        baseline_temps, baseline_dates; tag::AbstractString="default",
        extra::AbstractString="")
    cache_dir = get(ENV, "GUADEX_SPINUP_CACHE_DIR", joinpath("results", "sensitivity_spinup_cache"))
    reuse = _bool_setting(get(ENV, "GUADEX_SPINUP_REUSE", "1"))
    force = _bool_setting(get(ENV, "GUADEX_SPINUP_FORCE", "0"))
    path = joinpath(cache_dir, "spinup_" *
        _spinup_cache_key(settings; tag=tag, extra=extra,
            initial_state=initial_state, baseline_temps=baseline_temps) * ".jld2")

    if reuse && !force && isfile(path)
        cached = jldopen(path, "r") do f
            (state=collect(Float64, f["state"]), converged=Bool(f["converged"]),
             years=Int(f["years"]), last_basin_change=Float64(f["last_basin_change"]),
             last_q95_change=Float64(f["last_q95_change"]),
             last_relative_change=Float64(f["last_relative_change"]))
        end
        println("  burn-in: reusing cached state ($(cached.years) yr, " *
                "converged=$(cached.converged)) from $path")
        return cached
    end

    schedule = baseline_temps === nothing ? nothing :
        baseline_climatology_schedule(baseline_temps, baseline_dates, 1;
            days_per_year=Float64(settings.days_per_year))
    println("  burn-in: integrating under $(baseline_temps === nothing ? "static" : "seasonal") " *
            "baseline (max $(settings.spin_up_max_years) yr, tol=$(settings.spin_up_tol), " *
            "criterion=:$(settings.spin_up_criterion))...")
    spin = spin_up(params; initial_state=initial_state, schedule=schedule,
        days_per_year=Float64(settings.days_per_year),
        max_years=settings.spin_up_max_years, tol=settings.spin_up_tol,
        criterion=settings.spin_up_criterion,
        progress_every=settings.spin_up_progress_every)
    mkpath(cache_dir)
    jldsave(path; state=spin.state, converged=spin.converged, years=spin.years,
        last_basin_change=spin.last_basin_change, last_q95_change=spin.last_q95_change,
        last_relative_change=spin.last_relative_change)
    println("  burn-in: years=$(spin.years), converged=$(spin.converged), " *
            "basin change=$(spin.last_basin_change), q95=$(spin.last_q95_change), " *
            "max=$(spin.last_relative_change)")
    return spin
end

# ---------------------------------------------------------------------------
# Scenario integration and export
# ---------------------------------------------------------------------------

"""
    run_sensitivity_scenario(data_base, params, schedule, u0_flat, upstream_cost,
                             passability_vector; settings, saveat, positivity_cb)

Rebuild the dispersal network for one (upstream cost, passability) scenario,
solve the daily scheduled ODE and return `(sol, modified_dams, scenario_params)`.
"""
function run_sensitivity_scenario(data_base, params::MetacommunityParams,
        schedule::TemperatureSchedule, u0_flat, upstream_cost::Real,
        passability_vector; settings, saveat, positivity_cb)
    modified_dams = modified_dams_for_scenario(data_base.dams, passability_vector)
    new_dispersal_matrix = precompute_dispersal_matrix(
        params.n_sites,
        Matrix(data_base.distance_matrix),
        data_base.elevations,
        Float64(upstream_cost),
        modified_dams,
        data_base.species)
    scenario_params = set_dispersal_matrix(params, new_dispersal_matrix)
    t_end = Float64(settings.simulation_years * settings.days_per_year)
    prob = ODEProblem(metacommunity_ode_scheduled!, u0_flat, (0.0, t_end),
        ScheduledMetacommunityParams(scenario_params, schedule))
    sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, saveat=saveat,
        callback=positivity_cb)
    return sol, modified_dams, scenario_params
end

function sensitivity_saveat(settings)
    t_end = Float64(settings.simulation_years * settings.days_per_year)
    return unique(vcat(collect(0.0:settings.save_interval_days:t_end), t_end))
end

function sensitivity_positivity_callback()
    return DiscreteCallback(
        (u, t, integrator) -> any(x -> x < 0, u),
        integrator -> (integrator.u .= max.(integrator.u, 0.0));
        save_positions=(false, false))
end

"""
    export_sensitivity_run(run_dir; data_base, sol, params, forcing, settings, dams,
                           upstream_cost, baseline_species_density, native_species,
                           invasive_species, migratory_species, run_metadata)

Four-level + viewer + quasi-extinction export for one sensitivity run, using the
same reporting pipeline as the climate scenarios.
"""
function export_sensitivity_run(run_dir::AbstractString;
        data_base, sol, params::MetacommunityParams, forcing, settings, dams,
        upstream_cost::Real, baseline_species_density, native_species,
        invasive_species, migratory_species, run_metadata::AbstractDict)
    crosswalk_path = joinpath(@__DIR__, "data", "site_waterbody_crosswalk.csv")
    isfile(crosswalk_path) || (crosswalk_path = joinpath("data", "site_waterbody_crosswalk.csv"))
    isfile(crosswalk_path) || error("site water-body crosswalk not found: $crosswalk_path")
    return export_run_outputs(joinpath(run_dir, "export");
        sol_t=sol.t, sol_u=sol.u,
        sites=data_base.sites, species=data_base.species, site_df=data_base.site_df,
        crosswalk_path=crosswalk_path,
        native_species=native_species, invasive_species=invasive_species,
        migratory_species=migratory_species,
        days_per_year=settings.days_per_year, threshold=settings.presence_threshold,
        report_year_offsets=settings.year_offsets, report_year_labels=settings.year_labels,
        temperature_baseline=params.temperatures, warming=forcing.warming,
        dams=dams, distance_matrix=data_base.distance_matrix,
        habitat_suitability=params.habitat_suitability, upstream_cost=upstream_cost,
        require_crosswalk=true,
        baseline_species_density=baseline_species_density,
        quasi_extinction_q=settings.quasi_extinction_q,
        quasi_extinction_persistence=settings.quasi_extinction_persistence,
        species_upper_limits=params.thermal_upper_limits,
        daily_forcing=(temps=forcing.forcing_temps, dates=forcing.forcing_dates),
        run_metadata=run_metadata)
end

# ---------------------------------------------------------------------------
# Run identity, completion marker and provenance
# ---------------------------------------------------------------------------

const SENSITIVITY_COMPLETE_MARKER = ".sensitivity_complete.jld2"

# Artifacts a completed export must contain.  `level_basin.csv` alone is written
# before the species/viewer/metadata outputs, so it is not a completeness marker.
const SENSITIVITY_REQUIRED_EXPORTS = [
    "export/levels/level_basin.csv",
    "export/levels/species_timeseries.csv",
    "export/levels/quasi_extinction_summary.csv",
    "export/run_metadata.json",
]

"""
    sensitivity_run_fingerprint(settings, model; upstream_cost, passability_scenario,
                                interaction_matrix=nothing, sigma=nothing)

Science-affecting identity of one run.  Written with the completion marker and
compared against both the marker and, for legacy runs, the stored JLD2, so a run
produced under different settings is never reused or relabelled.
"""
function sensitivity_run_fingerprint(settings, model; upstream_cost, passability_scenario,
        interaction_matrix::Union{Nothing,AbstractString}=nothing,
        sigma::Union{Nothing,Real}=nothing)
    fingerprint = Dict{String,Any}(
        "start_year" => settings.start_year, "end_year" => settings.end_year,
        "climate_scenario" => String(model.scenario), "gcm" => String(model.gcm),
        "upstream_cost" => Float64(upstream_cost),
        "passability_scenario" => String(passability_scenario),
        "carrying_capacity_base_scaling" => settings.carrying_capacity_base_scaling,
        "carrying_capacity_scaling" => settings.carrying_capacity_scaling,
        "thermal_optima_fraction" => settings.thermal_optima_fraction,
        "spin_up_max_years" => settings.spin_up_max_years,
        "spin_up_tol" => settings.spin_up_tol,
        "spin_up_criterion" => String(settings.spin_up_criterion),
        "daily_forcing_file" => basename(String(settings.daily_forcing_file)),
    )
    interaction_matrix === nothing || (fingerprint["interaction_matrix"] = String(interaction_matrix))
    sigma === nothing || (fingerprint["thermal_sigma_multiplier"] = Float64(sigma))
    return fingerprint
end

function _stored_value_matches(actual, expected)
    if expected isa AbstractString
        return string(actual) == expected
    elseif expected isa Bool
        return Bool(actual) == expected
    elseif expected isa Integer
        return Int(actual) == Int(expected)
    end
    return isapprox(Float64(actual), Float64(expected); rtol=1e-9, atol=1e-12)
end

function _stored_matches(jld, fingerprint)
    for (key, expected) in fingerprint
        jldkey = key == "interaction_matrix" ? "interaction_matrix_type" : key
        haskey(jld, jldkey) || continue
        _stored_value_matches(jld[jldkey], expected) || return false
    end
    return true
end

"""
    sensitivity_run_status(run_dir, fingerprint)

`true` only when the run is complete *and* was produced with the same
science-affecting settings.  A run with a completion marker is compared against
that marker; a run without one (produced before markers existed) is accepted
only if its full export is present and the settings stored in its JLD2 match.
"""
function sensitivity_run_status(run_dir::AbstractString, fingerprint::AbstractDict)
    jld_path = joinpath(run_dir, "simulation_output.jld2")
    marker_path = joinpath(run_dir, SENSITIVITY_COMPLETE_MARKER)
    if isfile(jld_path) && isfile(marker_path)
        stored = jldopen(marker_path, "r") do f
            f["fingerprint"]
        end
        return stored == fingerprint
    end
    isfile(jld_path) || return false
    all(isfile(joinpath(run_dir, rel)) for rel in SENSITIVITY_REQUIRED_EXPORTS) || return false
    return jldopen(jld_path, "r") do jld
        _stored_matches(jld, fingerprint)
    end
end

"""
    mark_sensitivity_run_complete(run_dir, fingerprint)

Write the completion marker last, after the JLD2 and the full export succeeded.
"""
function mark_sensitivity_run_complete(run_dir::AbstractString, fingerprint::AbstractDict)
    jldsave(joinpath(run_dir, SENSITIVITY_COMPLETE_MARKER);
        fingerprint=fingerprint, completed_at=string(time()))
    return nothing
end

"""
    sensitivity_run_metadata(settings, model; script, upstream_cost, passability_scenario,
                             warming_end, heat_stress_k, spin_years, spin_converged,
                             obstacle_mode, extras)

Shared export provenance, so the two sweeps cannot drift.  `extras` carries the
script-specific keys (interaction matrix, thermal sigma).
"""
function sensitivity_run_metadata(settings, model; script::AbstractString,
        upstream_cost::Real, passability_scenario::AbstractString, warming_end::Real,
        heat_stress_k::Real, spin_years::Integer, spin_converged::Bool,
        obstacle_mode, extras::AbstractDict=Dict{String,Any}())
    metadata = Dict{String,Any}(
        "script" => String(script),
        "climate_scenario" => String(model.scenario), "gcm" => String(model.gcm),
        "start_year" => settings.start_year, "end_year" => settings.end_year,
        "upstream_cost" => Float64(upstream_cost),
        "passability_scenario" => String(passability_scenario),
        "warming_end_degc" => Float64(warming_end),
        "heat_stress_k" => Float64(heat_stress_k),
        "spin_up" => settings.spin_up,
        "spin_up_years" => Int(spin_years), "spin_up_converged" => Bool(spin_converged),
        "carrying_capacity_base_scaling" => settings.carrying_capacity_base_scaling,
        "carrying_capacity_scaling" => settings.carrying_capacity_scaling,
        "thermal_optima_fraction" => settings.thermal_optima_fraction,
        "obstacle_mode" => string(obstacle_mode),
    )
    merge!(metadata, extras)
    return metadata
end

"""
    observed_initial_state(u0_obs)

Fallback state used when the burn-in is disabled (`spin_up = false`).
"""
function observed_initial_state(u0_obs)
    return (state=vec(u0_obs), converged=false, years=0,
        last_basin_change=NaN, last_q95_change=NaN, last_relative_change=NaN)
end

# ---------------------------------------------------------------------------
# Summary index helpers
# ---------------------------------------------------------------------------

"""
    final_basin_metrics(run_dir, end_year)

Dict of the final-year basin-level metrics written to
`export/levels/level_basin.csv` (empty when the run has not been exported).
"""
function final_basin_metrics(run_dir::AbstractString, end_year::Int)
    path = joinpath(run_dir, "export", "levels", "level_basin.csv")
    isfile(path) || return Dict{String,Float64}()
    df = CSV.read(path, DataFrame)
    sub = df[df.year .== end_year, :]
    isempty(sub) && return Dict{String,Float64}()
    row = sub[1, :]
    out = Dict{String,Float64}()
    for col in names(df)
        val = row[Symbol(col)]
        val isa Real && (out[string(col)] = Float64(val))
    end
    return out
end

"""
    final_species_metrics(run_dir, species_code)

Final biomass, relative change and quasi-extinct fraction for one species from
`export/levels/quasi_extinction_summary.csv`.  Returns `(NaN, NaN, NaN)` when the
run has not been exported.
"""
function final_species_metrics(run_dir::AbstractString, species_code::AbstractString)
    path = joinpath(run_dir, "export", "levels", "quasi_extinction_summary.csv")
    isfile(path) || return (NaN, NaN, NaN)
    df = CSV.read(path, DataFrame)
    sub = df[df.species .== species_code, :]
    isempty(sub) && return (NaN, NaN, NaN)
    row = sub[1, :]
    return (Float64(row.final_biomass), Float64(row.relative_biomass_change),
        Float64(row.quasi_extinct_fraction))
end

"""
    empty_sensitivity_index()

Typed empty index DataFrame for a sensitivity sweep.
"""
function empty_sensitivity_index(; with_matrix::Bool=false)
    base = (
        climate_scenario=String[], gcm=String[], upstream_cost=Float64[],
        passability_scenario=String[], end_year=Int[], warming_end_degc=Float64[],
        basin_native_richness=Float64[], basin_native_extinction_risk=Float64[],
        basin_total_biomass=Float64[], basin_native_biomass=Float64[],
        basin_invasive_biomass=Float64[], basin_mean_temperature_c=Float64[],
        st_final_biomass=Float64[], st_relative_biomass_change=Float64[],
        st_quasi_extinct_fraction=Float64[], run_dir=String[])
    if with_matrix
        return DataFrame(; interaction_matrix=String[], thermal_sigma_multiplier=Float64[],
            base...)
    end
    return DataFrame(; base...)
end

"""
    index_row_values(run_dir, settings, model; upstream_cost, passability_scenario,
                     warming_end, matrix=nothing, sigma=nothing)

Assemble one summary-index row as a `Dict{String,Any}`.
"""
function index_row_values(run_dir::AbstractString, settings, model;
        upstream_cost::Real, passability_scenario::AbstractString, warming_end::Real,
        matrix::Union{Nothing,AbstractString}=nothing,
        sigma::Union{Nothing,Real}=nothing)
    basin = final_basin_metrics(run_dir, settings.end_year)
    st_final, st_rel, st_qe = final_species_metrics(run_dir, "ST")
    row = Dict{String,Any}(
        "climate_scenario" => model.scenario,
        "gcm" => model.gcm,
        "upstream_cost" => Float64(upstream_cost),
        "passability_scenario" => String(passability_scenario),
        "end_year" => settings.end_year,
        "warming_end_degc" => Float64(warming_end),
        "basin_native_richness" => get(basin, "mean_native_richness", NaN),
        "basin_native_extinction_risk" => get(basin, "mean_native_extinction_risk", NaN),
        "basin_total_biomass" => get(basin, "mean_total_biomass", NaN),
        "basin_native_biomass" => get(basin, "mean_native_biomass", NaN),
        "basin_invasive_biomass" => get(basin, "mean_invasive_biomass", NaN),
        "basin_mean_temperature_c" => get(basin, "mean_temperature_c", NaN),
        "st_final_biomass" => st_final,
        "st_relative_biomass_change" => st_rel,
        "st_quasi_extinct_fraction" => st_qe,
        "run_dir" => String(run_dir),
    )
    if matrix !== nothing
        row["interaction_matrix"] = String(matrix)
        row["thermal_sigma_multiplier"] = Float64(sigma)
    end
    return row
end

"""
    read_sensitivity_index(path; with_matrix=false)

Load an existing `runs_index.csv` into a DataFrame whose columns match
[`empty_sensitivity_index`](@ref) exactly and whose string columns are plain
`String` (CSV reading otherwise pools them as `InlineStrings.String15`, which
cannot hold labels longer than 15 characters such as `"invasive_favoring"`).
Missing columns are added as typed missing vectors so incremental appends keep
working across schema changes.
"""
function read_sensitivity_index(path::AbstractString; with_matrix::Bool=false)
    template = empty_sensitivity_index(; with_matrix=with_matrix)
    isfile(path) || return template
    raw = CSV.read(path, DataFrame)
    df = DataFrame()
    for col in names(template)
        T = eltype(template[!, col])
        if hasproperty(raw, col)
            v = raw[!, col]
            df[!, col] = T <: AbstractString ?
                String.(coalesce.(v, "")) : Vector{Union{Missing,T}}(coalesce.(v, missing))
        else
            df[!, col] = Vector{Union{Missing,T}}(missing, nrow(raw))
        end
    end
    return df
end

"""
    upsert_index_row!(index_df, row)

Replace any existing row with the same run key and append the new one, then
return the DataFrame (written by the caller after every run for resumability).
"""
function upsert_index_row!(index_df::DataFrame, row::AbstractDict)
    key_columns = ["climate_scenario", "gcm", "upstream_cost", "passability_scenario"]
    haskey(row, "interaction_matrix") && push!(key_columns, "interaction_matrix")
    haskey(row, "thermal_sigma_multiplier") && push!(key_columns, "thermal_sigma_multiplier")
    mask = trues(nrow(index_df))
    for col in key_columns
        cmp = index_df[!, Symbol(col)] .== row[col]
        mask .&= coalesce.(cmp, false)
    end
    !any(mask) || deleteat!(index_df, findall(mask))
    push!(index_df, row)
    return index_df
end
