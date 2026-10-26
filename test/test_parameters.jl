include(joinpath(@__DIR__, "..", "parameters.jl"))
using .SimulationParameters

@testset "SimulationParameters loader" begin
    params = SimulationParameters.load(joinpath(@__DIR__, "..", "parameters.toml"))
    SimulationParameters.require_sections(params,
        "general", "inputs", "obstacles", "species", "subcatchments",
        "scenarios", "run_model", "run_sensitivity",
        "run_sensitivity_report", "run_alt_interactions")

    @test params["general"]["days_per_year"] == 365

    # Obstacle/CEDEX defaults are the single source of truth for the scripts.
    @test params["obstacles"]["mode"] == "overlay"
    @test params["obstacles"]["matching_tolerance_m"] == 2000.0
    @test params["obstacles"]["upstream_passability"] == 0.1
    @test params["obstacles"]["downstream_passability"] == 0.5
    @test endswith(params["inputs"]["obstacles_file"], "Guadex.csv")
    @test endswith(params["inputs"]["cedex_var_file"], "2045_CEDEX_GUADALQUIVIR_VAR.csv")

    subcatchments = SimulationParameters.float_vector(params["subcatchments"]["all"])
    @test length(subcatchments) == 77
    @test all(isfinite, subcatchments)

    # WP0: ST (cold-water keystone) is native; AA/AAL/LR/MC are migratory and
    # deliberately excluded from both native and invasive richness.
    native = String.(params["species"]["native"])
    invasive = String.(params["species"]["invasive"])
    migratory = String.(get(params["species"], "migratory", String[]))
    @test length(native) == 10
    @test "ST" in native
    @test length(invasive) == 10
    @test sort(migratory) == ["AA", "AAL", "LR", "MC"]
    @test isempty(intersect(native, invasive))
    @test isempty(intersect(native, migratory))
    @test isempty(intersect(invasive, migratory))

    # Every scenario name referenced by an entry script must exist in the
    # scenario library (so typos fail loudly during configuration smoke tests).
    exploitation_library = params["scenarios"]["exploitation"]
    passability_library = params["scenarios"]["passability"]
    for (section, exploitation_names, passability_names) in (
        ("run_model", [params["run_model"]["exploitation_scenario"]], [params["run_model"]["passability_scenario"]]),
        ("run_sensitivity", params["run_sensitivity"]["exploitation_scenarios"], params["run_sensitivity"]["passability_scenarios"]),
        ("run_sensitivity_report", String[], params["run_sensitivity_report"]["passability_scenarios"]),
        ("run_alt_interactions", String[], params["run_alt_interactions"]["passability_scenarios"]),
    )
        for name in exploitation_names
            @test name in keys(exploitation_library)
        end
        for name in passability_names
            @test name in keys(passability_library)
        end
    end

    # Uniform-factor expansion matches the previous per-subcatchment dicts.
    exploitation = SimulationParameters.scenario_library(
        subcatchments, exploitation_library, params["run_sensitivity"]["exploitation_scenarios"])
    @test Set(keys(exploitation)) == Set(params["run_sensitivity"]["exploitation_scenarios"])
    @test all(==(1.0), values(exploitation["baseline"]))
    @test exploitation["high_exploitation"][1.1] == 0.5
    @test exploitation["high_exploitation"][41.0] == 0.5

    passability = SimulationParameters.scenario_library(
        subcatchments, passability_library, params["run_sensitivity_report"]["passability_scenarios"])
    @test passability["improved_passability"][1.1] == 1.5
    @test passability["blocked"][41.0] == 0.1

    # Per-subcatchment overrides are honoured on top of the uniform factor.
    overridden = SimulationParameters.scenario_dict(
        subcatchments, Dict("factor" => 1.0, "overrides" => Dict("1.1" => 0.7)))
    @test overridden[1.1] == 0.7
    @test overridden[41.0] == 1.0

    @test_throws Exception SimulationParameters.scenario_library(
        subcatchments, passability_library, ["does_not_exist"])
    @test_throws Exception SimulationParameters.expand_scenario(
        subcatchments, Dict("name" => "no_factor"))
end

@testset "committed k1x burn-in configuration" begin
    config_path = joinpath(@__DIR__, "..", "legacy", "parameters_climate_scenarios_k1x_burnin.toml")
    @test isfile(config_path)
    config = SimulationParameters.load(config_path)
    SimulationParameters.require_sections(config, "general", "inputs", "obstacles",
        "species", "subcatchments", "scenarios", "run_climate_scenarios", "temperature_stress")

    climate = config["run_climate_scenarios"]
    # Daily per-GCM forcing and no-warming control.
    @test climate["forcing_mode"] == "daily"
    @test climate["daily_forcing_per_gcm"] == true
    @test climate["control_run"] == true
    @test climate["baseline_period_start"] == 1986
    @test climate["baseline_period_end"] == 2005
    # 20-year horizon, K = 1x observed.
    @test climate["end_year"] - climate["start_year"] + 1 == 20
    @test climate["carrying_capacity_base_scaling"] == 1.0
    @test climate["carrying_capacity_scaling"] == 1.0
    @test climate["thermal_optima_fraction"] == 0.5
    # Burn-in: the September run converged at 373 years under the 600-year cap.
    @test climate["spin_up"] == true
    @test climate["spin_up_max_years"] == 600
    @test climate["spin_up_tol"] == 5.0e-5
    @test climate["spin_up_criterion"] == "basin"
    @test climate["spin_up_composition_tol"] == 1.0e-3
    @test climate["spin_up_min_years"] == 2
    @test climate["output_dir"] == "results/climate_scenarios_k1x_burnin"

    stress = config["temperature_stress"]
    @test stress["enabled"] == true
    @test stress["calibrate"] == false
    @test stress["k"] ≈ 3.7817692640400324e-5

    # The default file still ships the pre-September defaults other tests rely on.
    defaults = SimulationParameters.load(joinpath(@__DIR__, "..", "parameters.toml"))
    @test defaults["run_climate_scenarios"]["forcing_mode"] == "annual_mean"
    @test defaults["run_climate_scenarios"]["spin_up"] == false
    @test defaults["run_climate_scenarios"]["carrying_capacity_base_scaling"] == 10.0
    @test defaults["temperature_stress"]["enabled"] == false
    # E4: the shipped default stop rule is the stricter two-criteria rule, with
    # a composition tolerance and a minimum number of burn-in years.
    @test defaults["run_climate_scenarios"]["spin_up_criterion"] == "both"
    @test defaults["run_climate_scenarios"]["spin_up_composition_tol"] == 1.0e-2
    @test defaults["run_climate_scenarios"]["spin_up_min_years"] == 10
    @test defaults["obstacle_sensitivity"]["spin_up_criterion"] == "both"
    @test defaults["obstacle_sensitivity"]["spin_up_composition_tol"] == 1.0e-2
    @test defaults["obstacle_sensitivity"]["spin_up_min_years"] == 10
    # C3: the shipped defaults never select the interim route for the sweeps, so
    # their behaviour is unchanged (full burn-in).
    @test !haskey(defaults["obstacle_sensitivity"], "projection_route")
    @test defaults["obstacle_sensitivity"]["output_dir"] == "results/sensitivity_obstacles"
    @test !haskey(defaults["run_alt_interactions"], "output_dir")
    # C3: the default config does not select the interim route, so the full
    # burn-in (E4) remains the default behaviour everywhere.
    @test !haskey(defaults["run_climate_scenarios"], "projection_route")
end

@testset "corrected climate re-run configuration" begin
    config_path = joinpath(@__DIR__, "..", "parameters_climate_scenarios_corrected.toml")
    @test isfile(config_path)
    config = SimulationParameters.load(config_path)
    SimulationParameters.require_sections(config, "general", "inputs", "obstacles",
        "species", "subcatchments", "scenarios", "obstacle_sensitivity",
        "run_climate_scenarios", "temperature_stress")

    climate = config["run_climate_scenarios"]
    # Daily per-GCM forcing, corrected baseline window and no-warming control.
    @test climate["forcing_mode"] == "daily"
    @test climate["daily_forcing_per_gcm"] == true
    @test climate["control_run"] == true
    @test climate["baseline_period_start"] == 1986
    @test climate["baseline_period_end"] == 2005
    # E4 corrected stop rule: `both` (basin + composition), NOT the legacy
    # September `basin` rule.
    @test climate["spin_up"] == true
    @test climate["spin_up_criterion"] == "both"
    @test climate["spin_up_composition_tol"] == 1.0e-2
    @test climate["spin_up_min_years"] == 10
    # Corrected output root, distinct from BOTH the legacy reproduction
    # directory (`results/climate_scenarios_k1x_burnin`) and the pre-correction
    # September outputs (`results/climate_scenarios`), so neither is overwritten.
    @test climate["output_dir"] == "results/climate_scenarios_corrected"
    @test climate["carrying_capacity_base_scaling"] == 1.0
    # C3 interim projection route selected, with a short fixed spin-up and the
    # mandatory matched no-warming control.
    @test climate["projection_route"] == "interim_observed"
    @test 0 <= climate["interim_spin_up_years"] <= 5
    @test climate["control_run"] == true

    # Heat stress calibrated from the corrected water-temperature baseline
    # rather than the pinned September slope.
    stress = config["temperature_stress"]
    @test stress["enabled"] == true
    @test stress["calibrate"] == true

    # The climate-consistent sensitivity sweeps (obstacle_sensitivity) share the
    # corrected stop rule, so the two code paths cannot diverge.
    sensitivity = config["obstacle_sensitivity"]
    @test sensitivity["spin_up_criterion"] == "both"
    @test sensitivity["spin_up_composition_tol"] == 1.0e-2
    @test sensitivity["spin_up_min_years"] == 10
    @test sensitivity["heat_stress_calibrate"] == true
    # C3: the sweeps share the SAME interim projection route as the climate run,
    # with the corrected output roots that preserve the legacy September dirs.
    @test sensitivity["projection_route"] == "interim_observed"
    @test sensitivity["interim_spin_up_years"] == 3
    @test sensitivity["output_dir"] == "results/sensitivity_obstacles_corrected"
    @test config["run_alt_interactions"]["output_dir"] == "results/alt_interactions_corrected"
    # (c) The sweeps inherit the C4 on-path graph and the team biological options.
    @test SimulationParameters.connectivity_method(config) == :on_path
    @test SimulationParameters.biological_options(config).pool_capacity_mode == :cap

    # The legacy reproduction stays pinned to `basin` and differs on the stop
    # rule, so the two configurations cannot be confused.
    legacy = SimulationParameters.load(
        joinpath(@__DIR__, "..", "legacy", "parameters_climate_scenarios_k1x_burnin.toml"))
    @test legacy["run_climate_scenarios"]["spin_up_criterion"] == "basin"
    @test legacy["temperature_stress"]["calibrate"] == false
    @test climate["spin_up_criterion"] != legacy["run_climate_scenarios"]["spin_up_criterion"]
    # The legacy reproduction does NOT select the C3 route, so it keeps the
    # full-burn-in behaviour (the route is an explicit opt-in only).
    @test !haskey(legacy["run_climate_scenarios"], "projection_route")
end

@testset "biological-assumption option loader (E13-E16, E18)" begin
    defaults = SimulationParameters.load(joinpath(@__DIR__, "..", "parameters.toml"))
    bio = SimulationParameters.biological_options(defaults)
    @test bio.absence_growth_fraction == 0.1
    @test bio.pool_capacity_mode == :legacy
    @test bio.fishless_dificil_capacity == :legacy
    @test bio.nonreproducing_local_growth == :legacy
    @test bio.exclude_fishfarm_eel_records == false
    @test bio.salinity_envelope == false

    flat = SimulationParameters.biological_options_dict(defaults)
    @test flat["pool_capacity_mode"] == "legacy"
    @test flat["absence_growth_fraction"] == 0.1
    @test flat["salinity_envelope"] == false

    # An absent section resolves to exactly the same legacy defaults.
    @test SimulationParameters.biological_options(Dict{String,Any}()) == bio

    # Overrides parse and validate.
    parsed = SimulationParameters.biological_options(Dict("biological_options" => Dict(
        "absence_growth_fraction" => 1.0,
        "pool_capacity_mode" => "exclude",
        "fishless_dificil_capacity" => "near_zero",
        "nonreproducing_local_growth" => "zero",
        "exclude_fishfarm_eel_records" => true,
        "salinity_envelope" => true)))
    @test parsed.absence_growth_fraction == 1.0
    @test parsed.pool_capacity_mode == :exclude
    @test parsed.fishless_dificil_capacity == :near_zero
    @test parsed.nonreproducing_local_growth == :zero
    @test parsed.exclude_fishfarm_eel_records == true
    @test parsed.salinity_envelope == true

    @test_throws ErrorException SimulationParameters.biological_options(Dict(
        "biological_options" => Dict("pool_capacity_mode" => "bogus")))
    @test_throws ErrorException SimulationParameters.biological_options(Dict(
        "biological_options" => Dict("fishless_dificil_capacity" => "bogus")))
    @test_throws ErrorException SimulationParameters.biological_options(Dict(
        "biological_options" => Dict("nonreproducing_local_growth" => "bogus")))
    @test_throws ErrorException SimulationParameters.biological_options(Dict(
        "biological_options" => Dict("absence_growth_fraction" => -1.0)))

    # The corrected re-run config carries the TEAM-CHOSEN values (E13-E16, E18);
    # `parameters.toml` keeps the safe legacy defaults asserted above.
    corrected = SimulationParameters.load(
        joinpath(@__DIR__, "..", "parameters_climate_scenarios_corrected.toml"))
    chosen = SimulationParameters.biological_options(corrected)
    @test chosen.absence_growth_fraction == 1.0
    @test chosen.pool_capacity_mode == :cap
    @test chosen.fishless_dificil_capacity == :near_zero
    @test chosen.nonreproducing_local_growth == :zero
    @test chosen.exclude_fishfarm_eel_records == true
    @test chosen.salinity_envelope == true
    @test chosen != bio

    flat_chosen = SimulationParameters.biological_options_dict(corrected)
    @test flat_chosen["absence_growth_fraction"] == 1.0
    @test flat_chosen["pool_capacity_mode"] == "cap"
    @test flat_chosen["fishless_dificil_capacity"] == "near_zero"
    @test flat_chosen["nonreproducing_local_growth"] == "zero"
    @test flat_chosen["exclude_fishfarm_eel_records"] == true
    @test flat_chosen["salinity_envelope"] == true
end
