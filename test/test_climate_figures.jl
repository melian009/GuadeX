using Test
using DataFrames
using CSV
using Statistics
using Guadex

# =============================================================================
# Synthetic four-level run tables (no simulation needed)
# =============================================================================

function _write_run(root, scenario, gcm, years, basin_values)
    level_dir = joinpath(root, scenario, gcm, "export", "levels")
    mkpath(level_dir)

    basin = DataFrame(
        basin=fill("ES050", length(years)),
        year=years,
        n_sites=fill(10, length(years)),
        mean_native_richness=basin_values,
        mean_native_extinction_risk=fill(0.1, length(years)),
        mean_total_biomass=basin_values .* 10.0,
    )
    CSV.write(joinpath(level_dir, "level_basin.csv"), basin)

    unit_years = repeat(years; inner=2)
    unit_values = Float64[]
    for v in basin_values
        push!(unit_values, v)
        push!(unit_values, v + 2)
    end
    subcatchment = DataFrame(
        subcatchment=repeat(["A", "B"]; outer=length(years)),
        year=unit_years,
        n_sites=fill(5, length(unit_years)),
        mean_native_richness=unit_values,
        mean_native_extinction_risk=fill(0.1, length(unit_years)),
        mean_total_biomass=unit_values .* 10.0,
    )
    CSV.write(joinpath(level_dir, "level_subcatchment.csv"), subcatchment)

    water_body = DataFrame(
        water_body=repeat(["W1", "W2"]; outer=length(years)),
        water_body_name=repeat(["Body 1", "Body 2"]; outer=length(years)),
        year=unit_years,
        n_sites=fill(5, length(unit_years)),
        mean_native_richness=unit_values,
        mean_native_extinction_risk=fill(0.1, length(unit_years)),
        mean_total_biomass=unit_values .* 10.0,
    )
    CSV.write(joinpath(level_dir, "level_water_body.csv"), water_body)

    sites = DataFrame(
        year=repeat(years; inner=2),
        CODIGO=["1.1.1", "1.1.2", "1.1.1", "1.1.2"],
        native_richness=[1.0, 3.0, 2.0, 4.0],
        invasive_richness=[1.0, 1.0, 1.0, 2.0],
        total_richness=[2.0, 4.0, 3.0, 6.0],
        native_biomass=[5.0, 7.0, 6.0, 9.0],
        invasive_biomass=[1.0, 2.0, 1.0, 3.0],
        total_biomass=[6.0, 9.0, 7.0, 12.0],
        native_richness_relative=[1.0, 1.0, 1.0, 1.0],
        native_extinction_risk=[0.0, 0.0, 0.0, 0.0],
        temperature_c=[15.0, 16.0, 16.0, 17.0],
        delta_temperature_c=[0.0, 0.0, 1.0, 1.0],
    )
    CSV.write(joinpath(level_dir, "level_sampling_point.csv"), sites)
    return level_dir
end

@testset "climate figures helpers" begin
    root = mktempdir()
    _write_run(root, "ssp126", "GCM1", [2026, 2045], [2.0, 4.0])
    _write_run(root, "ssp126", "GCM2", [2026, 2045], [2.0, 6.0])
    _write_run(root, "ssp370", "GCMX", [2026, 2027], [1.0, 2.0])
    CSV.write(joinpath(root, "runs_index.csv"), DataFrame(
        scenario=["ssp126", "ssp126", "ssp370"],
        gcm=["GCM1", "GCM2", "GCMX"],
        end_year=[2045, 2045, 2027],
    ))

    runs, expected = Guadex.discover_climate_runs(root)
    @test expected == 2045
    @test length(runs) == 3
    complete = [r for r in runs if r.complete]
    @test length(complete) == 2
    incomplete = only([r for r in runs if !r.complete])
    @test incomplete.scenario == "ssp370"
    @test incomplete.gcm == "GCMX"

    # Level tables load from the synthetic CSVs; missing files return nothing.
    level_dir = joinpath(root, "ssp126", "GCM1", "export", "levels")
    basin_table = Guadex.load_level_table(joinpath(level_dir, "level_basin.csv"), :basin)
    @test basin_table !== nothing
    @test issubset(Set([:year, :mean_native_richness]), Set(Symbol.(names(basin_table))))
    @test Guadex.load_level_table(joinpath(level_dir, "nope.csv"), :basin) === nothing

    gcm1 = only([r for r in runs if r.gcm == "GCM1"])
    @test Guadex.run_level_stats(gcm1, :basin, :native_richness).mean == [2.0, 4.0]
    sub = Guadex.run_level_stats(gcm1, :subcatchment, :native_richness)
    @test sub.p10[end] < sub.p50[end] < sub.p90[end]

    ens = Guadex.ensemble_level_stats(complete, :basin, :native_richness)
    @test ens.years == [2026, 2045]
    @test ens.median == [2.0, 5.0]
    @test ens.p10[end] ≈ 4.2 atol = 1e-9
    @test ens.p90[end] ≈ 5.8 atol = 1e-9
    @test ens.n_models == [2, 2]

    @test Guadex.climate_metric_column(:site, :native_richness) == :native_richness
    @test Guadex.climate_metric_column(:basin, :native_richness) == :mean_native_richness

    # Missing level column → all NaN, no error.
    unknown = Guadex.run_level_stats(gcm1, :basin, :not_a_metric)
    @test all(isnan, unknown.mean)
end

@testset "ensemble nested-level central line is not pinned" begin
    root = mktempdir()
    scenario, gcm = "ssp126", "GCM1"
    level_dir = joinpath(root, scenario, gcm, "export", "levels")
    mkpath(level_dir)
    years = [2026, 2045]

    # 6 of 10 sub-basins stay at 1.0 and 4 rise to 10.0, so the pooled median
    # stays pinned at the modal value while the level aggregate changes.
    sub_2026 = fill(1.0, 10)
    sub_2045 = vcat(fill(1.0, 6), fill(10.0, 4))
    CSV.write(joinpath(level_dir, "level_subcatchment.csv"), DataFrame(
        subcatchment=repeat(string.(1:10); outer=2),
        year=repeat(years; inner=10),
        n_sites=fill(1, 20),
        mean_native_richness=vcat(sub_2026, sub_2045),
        mean_native_extinction_risk=zeros(20),
        mean_total_biomass=vcat(sub_2026, sub_2045) .* 10.0,
    ))
    CSV.write(joinpath(level_dir, "level_basin.csv"), DataFrame(
        basin=fill("ES050", 2),
        year=years,
        n_sites=fill(10, 2),
        mean_native_richness=[1.6, 4.6],
        mean_native_extinction_risk=zeros(2),
        mean_total_biomass=[16.0, 46.0],
    ))
    CSV.write(joinpath(root, "runs_index.csv"),
        DataFrame(scenario=[scenario], gcm=[gcm], end_year=[2045]))

    runs, _ = Guadex.discover_climate_runs(root)
    run = only(runs)
    ens = Guadex.ensemble_level_stats([run], :subcatchment, :native_richness)
    # Each run is reduced to its mean across sub-basins, so the modal unit value
    # (1.0) never pins the line: 2026 mean = 1.0, 2045 mean = (6*1 + 4*10)/10.
    @test ens.mean[1] ≈ 1.0
    @test ens.mean[2] ≈ 4.6
    @test ens.median[2] ≈ 4.6
    @test ens.p10[2] ≈ 4.6
end

@testset "climate diagnostics data readers" begin
    root = mktempdir()
    _write_run(root, "ssp126", "GCM1", [2026, 2045], [2.0, 4.0])
    _write_run(root, "ssp126", "GCM2", [2026, 2045], [2.0, 6.0])
    CSV.write(joinpath(root, "runs_index.csv"), DataFrame(
        scenario=["ssp126", "ssp126"],
        gcm=["GCM1", "GCM2"],
        end_year=[2045, 2045],
    ))

    series = Guadex.read_climate_basin_series(root)
    @test length(series) == 2
    @test all(s -> s.scenario == "ssp126", series)
    gcm1 = only([s for s in series if s.gcm == "GCM1"])
    @test gcm1.years == [2026, 2045]
    @test gcm1.native_richness == [2.0, 4.0]
    @test gcm1.total_biomass == [20.0, 40.0]
    # Columns absent from the synthetic table read as NaN rather than failing.
    @test all(isnan, gcm1.delta_temperature_c)
    @test all(isnan, gcm1.invasive_richness)

    # An empty tree yields no series instead of an error.
    @test isempty(Guadex.read_climate_basin_series(mktempdir()))
end
