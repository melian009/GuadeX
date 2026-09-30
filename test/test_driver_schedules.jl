# =============================================================================
# Driver daily schedules reconstruct the corrected series (ITEM 1).
#
# `run_climate_experiments.jl` (`experiment_daily_schedule`) and
# `sensitivity_core.jl` (`load_climate_model_forcing`) must reference the daily
# anomaly to the SAME per-site level stored as
# `MetacommunityParams.temperatures` (the corrected `tw_baseline_mean`), as
# `run_climate_scenarios.jl` does.  Differencing against the forcing file's own
# `historical` mean instead leaves a residual offset.  This file exercises the
# shared invariant and the `sensitivity_core` driver directly.
# =============================================================================
using Test
using DataFrames
using CSV
using Dates
using Statistics

# The sensitivity driver helpers are plain functions defined by
# `sensitivity_core.jl`; include them so the driver code path itself is tested.
include(joinpath(@__DIR__, "..", "sensitivity_core.jl"))

@testset "driver daily schedules reconstruct the corrected series (ITEM 1)" begin
    sites = ["s1", "s2"]
    # Corrected model level (the corrected `tw_baseline_mean`); deliberately
    # different from the forcing file's own `historical` window mean below.
    levels = [12.0, 15.0]

    hist_dates = [Date(1986, 1, 1) + Day(k) for k in 0:7299]   # 1986-2005
    sim_dates = [Date(2026, 1, 1) + Day(k) for k in 0:7299]    # 2026-2045
    hist = DataFrame(
        date=hist_dates, scenario=fill("historical", length(hist_dates)),
        s1=fill(10.0, length(hist_dates)), s2=fill(13.0, length(hist_dates)))
    scen = DataFrame(
        date=sim_dates, scenario=fill("ssp585", length(sim_dates)),
        s1=[12.0 + 0.001 * k for k in 1:length(sim_dates)],
        s2=[15.0 + 0.002 * k for k in 1:length(sim_dates)])

    tmp = mktempdir()
    base = joinpath(tmp, "forcing.csv")
    per_gcm = joinpath(tmp, "forcing_ssp585_TESTGCM.csv")
    CSV.write(per_gcm, vcat(hist, scen))

    src = Guadex.load_daily_forcing_any(per_gcm)
    fdates, ftemps = Guadex.wide_forcing_matrix(src, sites; scenario="ssp585")
    hdates, htemps = Guadex.wide_forcing_matrix(src, sites; scenario="historical")

    forcing = load_climate_model_forcing(base, sites, "ssp585", "TESTGCM";
        baseline_start=1986, baseline_end=2005,
        baseline_means=levels, days_per_year=365.0,
        year_labels=collect(2026:2045))

    # The driver records the level it actually used (the corrected model level).
    @test forcing.baseline_means == levels

    # level + anomaly reconstructs the absolute daily series on sample days.
    for i in eachindex(sites), day in (1, 2, 500, length(fdates))
        @test levels[i] +
            Guadex.temperature_delta(forcing.schedule, i, Float64(day - 1)) ≈
            ftemps[i, day]
    end

    # `run_climate_experiments.jl` builds its schedule with exactly this call
    # (same function, same `baseline_means`), whether it reads a per-GCM
    # candidate or the ensemble-median fallback, so assert it directly too.
    sched_exp, means_exp = Guadex.daily_temperature_schedule(ftemps;
        dates=fdates, baseline_start=1986, baseline_end=2005,
        baseline_means=levels)
    @test means_exp == levels
    for i in eachindex(sites), day in (1, 1000, length(fdates))
        @test levels[i] +
            Guadex.temperature_delta(sched_exp, i, Float64(day - 1)) ≈
            ftemps[i, day]
    end

    # A schedule differenced against the file's own `historical` mean would
    # leave exactly the corrected level's offset (2 degC here), so the two
    # schedules must differ by that offset at the first node.
    legacy, legacy_means = Guadex.daily_temperature_schedule(ftemps; dates=fdates,
        baseline_start=1986, baseline_end=2005,
        baseline_temps=htemps, baseline_dates=hdates)
    @test legacy_means == [10.0, 13.0]
    @test isapprox(Guadex.temperature_delta(legacy, 1, 0.0) -
                   Guadex.temperature_delta(forcing.schedule, 1, 0.0),
        2.0; atol=1e-9)
end
