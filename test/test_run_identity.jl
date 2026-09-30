using Test
using SparseArrays
using Dates
using JLD2
using Guadex

# `sensitivity_core.jl` is a plain `include`d file (not part of the Guadex
# module); include it here so the burn-in cache key and run fingerprint can be
# exercised directly.  It only references names already available in this test
# environment (`using Guadex`, `using SparseArrays`, Statistics).
include(joinpath(@__DIR__, "..", "sensitivity_core.jl"))

function _digest_test_params()
    return Guadex.MetacommunityParams(
        2, 2, [1.0 0.3; -0.2 1.0], sparse([1, 2], [1, 2], [0.0, 0.0], 2, 2),
        [0.5, 1.0], [0.1 0.05; 0.02 0.1], [20.0, 18.0], [1.0, 1.0],
        [15.0, 14.0], [3.0, 3.5], [50.0, 40.0])
end

@testset "parameter digest" begin
    params = _digest_test_params()
    options = Dict{String,Any}("forcing_mode" => "daily", "baseline_period" => "1986-2005")
    base = Guadex.parameter_digest(params; options=options)
    @test base isa String
    @test length(base) == 16
    # Stable across repeated calls with identical inputs.
    @test base == Guadex.parameter_digest(params; options=options)
    # Option insertion order is irrelevant, and the code version participates.
    @test Guadex.parameter_digest(params; options=Dict("a" => 1, "b" => 2)) ==
        Guadex.parameter_digest(params; options=Dict("b" => 2, "a" => 1))
    @test Guadex.parameter_digest(params; code="aaaa") !=
        Guadex.parameter_digest(params; code="bbbb")
    # A changed option value changes the digest.
    @test Guadex.parameter_digest(params; options=Dict("k" => 1.0)) !=
        Guadex.parameter_digest(params; options=Dict("k" => 10.0))
    # Interaction matrix.
    changed_matrix = Guadex.set_interaction_matrix(params, [1.0 0.9; -0.2 1.0])
    @test Guadex.parameter_digest(changed_matrix) != Guadex.parameter_digest(params)
    # Scalar model parameters: heat-stress slope and carrying-capacity scaling.
    @test Guadex.parameter_digest(Guadex.set_heat_stress_rate(params, 0.02)) !=
        Guadex.parameter_digest(params)
    @test Guadex.parameter_digest(Guadex.scale_carrying_capacity(params, 3.0)) !=
        Guadex.parameter_digest(params)
    # Thermal optima and sigma.
    @test Guadex.parameter_digest(Guadex.set_thermal_optima(params, [10.0, 20.0])) !=
        Guadex.parameter_digest(params)
    @test Guadex.parameter_digest(Guadex.set_thermal_sigma_multiplier(params, 0.3)) !=
        Guadex.parameter_digest(params)
end

@testset "burn-in cache key includes the parameter digest" begin
    settings = (
        days_per_year = 365,
        carrying_capacity_base_scaling = 1.0,
        carrying_capacity_scaling = 1.0,
        thermal_optima_fraction = 0.5,
        heat_stress_enabled = true,
        heat_stress_k = 0.0,
        heat_stress_calibrate = true,
        heat_stress_max_loss = 0.05,
        spin_up_max_years = 600,
        spin_up_tol = 5.0e-5,
        spin_up_composition_tol = 1.0e-3,
        spin_up_min_years = 10,
        spin_up_criterion = :basin,
        baseline_period_start = 1986,
        baseline_period_end = 2005,
        daily_forcing_file = "no_such_forcing_file.csv",
        upstream_cost_default = 0.05,
    )
    args = (tag = "t", initial_state = [1.0, 2.0], baseline_temps = [10.0 11.0; 9.0 8.0])
    key_a = _spinup_cache_key(settings; param_digest="aaaa", args...)
    key_b = _spinup_cache_key(settings; param_digest="aaaa", args...)
    key_c = _spinup_cache_key(settings; param_digest="bbbb", args...)
    @test key_a == key_b
    @test key_a != key_c
end

@testset "run fingerprint includes the parameter digest" begin
    params = _digest_test_params()
    settings = (
        start_year = 2026, end_year = 2045,
        carrying_capacity_base_scaling = 1.0,
        carrying_capacity_scaling = 1.0,
        thermal_optima_fraction = 0.5,
        spin_up_max_years = 600,
        spin_up_tol = 5.0e-5,
        spin_up_composition_tol = 1.0e-3,
        spin_up_min_years = 10,
        spin_up_criterion = :basin,
        daily_forcing_file = "no_such_forcing_file.csv",
        baseline_period_start = 1986,
        baseline_period_end = 2005,
        spin_up = true,
        upstream_cost_default = 0.05,
    )
    options = Dict{String,Any}("obstacle_mode" => "overlay")
    model = (scenario = "ssp126", gcm = "ACCESS-CM2")
    fp = sensitivity_run_fingerprint(settings, model;
        upstream_cost=0.05, passability_scenario="baseline", params=params,
        digest_options=options)
    @test haskey(fp, "parameter_digest")
    # Changing one model parameter changes the fingerprint.
    changed = Guadex.set_heat_stress_rate(params, 0.02)
    fp_changed = sensitivity_run_fingerprint(settings, model;
        upstream_cost=0.05, passability_scenario="baseline", params=changed,
        digest_options=options)
    @test fp != fp_changed
    @test fp["parameter_digest"] != fp_changed["parameter_digest"]
    # Without params the legacy fingerprint shape is preserved.
    fp_legacy = sensitivity_run_fingerprint(settings, model;
        upstream_cost=0.05, passability_scenario="baseline")
    @test !haskey(fp_legacy, "parameter_digest")
end

@testset "sensitivity interim projection route wiring (C3)" begin
    minimal_raw(section) = Dict{String,Any}(
        "general" => Dict{String,Any}("days_per_year" => 365),
        "run_climate_scenarios" => Dict{String,Any}(),
        "obstacle_sensitivity" => section,
    )

    # Default: no key anywhere -> full burn-in (nothing else changes).
    d = load_sensitivity_settings(minimal_raw(Dict{String,Any}()))
    @test d.projection_route === :full_burnin

    # Explicit C3 interim route in [obstacle_sensitivity].
    s = load_sensitivity_settings(minimal_raw(Dict{String,Any}(
        "projection_route" => "interim_observed", "interim_spin_up_years" => 2)))
    @test s.projection_route === :interim_observed
    @test s.interim_spin_up_years == 2

    # A typo errors loudly instead of silently changing the route.
    @test_throws Exception load_sensitivity_settings(
        minimal_raw(Dict{String,Any}("projection_route" => "bogus")))

    # The route participates in the burn-in digest and the cache key, so a
    # full-burn-in equilibrium can never be reused for an interim run.
    complete = (
        spin_up = true,
        spin_up_max_years = 600, spin_up_tol = 5.0e-5,
        spin_up_composition_tol = 1.0e-2, spin_up_min_years = 10,
        spin_up_criterion = :both, carrying_capacity_base_scaling = 1.0,
        carrying_capacity_scaling = 1.0, thermal_optima_fraction = 0.5,
        heat_stress_enabled = false, heat_stress_k = 0.0,
        heat_stress_calibrate = false, heat_stress_max_loss = 0.05,
        upstream_cost_default = 0.05, daily_forcing_file = "no_such_file.csv",
        baseline_period_start = 1986, baseline_period_end = 2005,
        projection_route = :full_burnin, interim_spin_up_years = 0)
    inter = merge(complete, (projection_route = :interim_observed,
        interim_spin_up_years = 3))
    @test _burnin_digest_options(complete, "x") != _burnin_digest_options(inter, "x")
    @test _spinup_cache_key(complete; param_digest = "aaaa") !=
        _spinup_cache_key(inter; param_digest = "aaaa")

    # Integration: the interim route integrates EXACTLY the requested number of
    # baseline years from the observed start and bypasses the E4 stop rule.
    params = _digest_test_params()
    initial_state = [1.0, 2.0, 3.0, 4.0]
    baseline_dates = [Date(1986, 1, 1) + Day(k) for k in 0:729]
    baseline_temps = [12.0 + 0.01 * (k % 365) for _ in 1:2, k in 0:729]
    settings = load_sensitivity_settings(minimal_raw(Dict{String,Any}(
        "spin_up" => true, "spin_up_max_years" => 600, "spin_up_tol" => 5.0e-5,
        "spin_up_composition_tol" => 1.0e-2, "spin_up_min_years" => 10,
        "spin_up_criterion" => "both", "heat_stress_enabled" => false,
        "daily_forcing_file" => "no_such_file.csv",
        "projection_route" => "interim_observed", "interim_spin_up_years" => 2)))
    cache_dir = mktempdir()
    old_dir, old_reuse = get(ENV, "GUADEX_SPINUP_CACHE_DIR", nothing),
        get(ENV, "GUADEX_SPINUP_REUSE", nothing)
    ENV["GUADEX_SPINUP_CACHE_DIR"] = cache_dir
    ENV["GUADEX_SPINUP_REUSE"] = "0"
    try
        spin = get_or_compute_spinup(settings, params, initial_state,
            baseline_temps, baseline_dates; tag = "interim_test", extra = "x")
        @test spin.years == 2                       # exactly the fixed relaxation
        @test length(spin.state) == 4
        @test spin.state != initial_state           # integration actually ran
        @test !spin.converged                       # E4 stop rule bypassed
        @test !isempty(filter(f -> endswith(f, ".jld2"), readdir(cache_dir)))

        # Route provenance reaches the run metadata.
        meta = sensitivity_run_metadata(settings,
            (scenario = "ssp126", gcm = "TEST"); script = "test.jl",
            upstream_cost = 0.05, passability_scenario = "baseline",
            warming_end = 1.0, heat_stress_k = 0.0, spin_years = spin.years,
            spin_converged = spin.converged, obstacle_mode = :overlay)
        @test meta["projection_route"] == "interim_observed"
        @test meta["interim_spin_up_years"] == 2
    finally
        old_dir === nothing ? delete!(ENV, "GUADEX_SPINUP_CACHE_DIR") :
            (ENV["GUADEX_SPINUP_CACHE_DIR"] = old_dir)
        old_reuse === nothing ? delete!(ENV, "GUADEX_SPINUP_REUSE") :
            (ENV["GUADEX_SPINUP_REUSE"] = old_reuse)
    end
end

@testset "burn-in diagnostics flow into run metadata and index" begin
    settings = (
        start_year = 2026, end_year = 2045,
        carrying_capacity_base_scaling = 1.0,
        carrying_capacity_scaling = 1.0,
        thermal_optima_fraction = 0.5,
        spin_up = true,
        spin_up_max_years = 600,
        spin_up_tol = 5.0e-5,
        spin_up_composition_tol = 1.0e-3,
        spin_up_min_years = 10,
        spin_up_criterion = :both,
        daily_forcing_file = "no_such_forcing_file.csv",
        baseline_period_start = 1986,
        baseline_period_end = 2005,
        upstream_cost_default = 0.05,
    )
    model = (scenario = "ssp126", gcm = "ACCESS-CM2")
    spin = (years = 494, converged = true, last_basin_change = 3.0e-6,
        last_composition_change = 4.0e-4, last_composition_change_all_cells = 9.0e-3,
        composition_active_cells = 2500, criteria = [:basin, :composition],
        last_q95_change = 1.0e-3, last_relative_change = 1.0e-2)
    meta = sensitivity_run_metadata(settings, model;
        script = "test.jl", upstream_cost = 0.05, passability_scenario = "baseline",
        warming_end = 1.0, heat_stress_k = 0.0,
        spin_years = spin.years, spin_converged = spin.converged,
        spin_basin_change = spin.last_basin_change,
        spin_composition_change = spin.last_composition_change,
        spin_composition_change_all_cells = spin.last_composition_change_all_cells,
        spin_composition_active_cells = spin.composition_active_cells,
        spin_criteria = spin.criteria, obstacle_mode = :overlay)
    @test meta["spin_up_years"] == 494
    @test meta["spin_up_total_biomass_change"] == 3.0e-6
    @test meta["spin_up_composition_q95_change"] == 4.0e-4
    @test meta["spin_up_composition_q95_change_all_cells"] == 9.0e-3
    @test meta["spin_up_composition_active_cells"] == 2500
    @test meta["spin_up_criteria"] == ["basin", "composition"]
    @test meta["spin_up_composition_tol"] == 1.0e-3
    @test meta["spin_up_composition_active_floor"] == 0.1
    @test meta["spin_up_min_years"] == 10

    row = index_row_values(joinpath(@__DIR__, "no_such_run_dir"), settings, model;
        upstream_cost = 0.05, passability_scenario = "baseline", warming_end = 1.0,
        spin = spin)
    @test row["spin_up_years"] == 494
    @test row["spin_up_total_biomass_change"] == 3.0e-6
    @test row["spin_up_composition_q95_change"] == 4.0e-4
    @test row["spin_up_composition_q95_change_all_cells"] == 9.0e-3
    @test row["spin_up_composition_active_cells"] == 2500
    @test row["spin_up_composition_active_floor"] == 0.1
    @test row["spin_up_composition_tol"] == 1.0e-3
    @test row["spin_up_criteria"] == "basin;composition"
end
