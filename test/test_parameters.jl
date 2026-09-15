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
