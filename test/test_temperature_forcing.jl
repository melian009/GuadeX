using Test
using DataFrames
using CSV
using Dates
using Statistics
using SparseArrays
using LinearAlgebra
using DifferentialEquations
using Guadex

# =============================================================================
# 1. Empirical thermal limits (WP3)
# =============================================================================
@testset "parse_temperature_range_and_sigma limits" begin
    p = Guadex.parse_temperature_range_and_sigma("4 to 20")
    @test p.optimum == 12.0
    @test p.sigma ≈ 16.0 / 6.0
    @test p.lower == 4.0
    @test p.upper == 20.0

    p2 = Guadex.parse_temperature_range_and_sigma("")
    @test p2.lower == -Inf
    @test p2.upper == Inf
end

@testset "load_species_characteristics exposes limits" begin
    df = Guadex.load_species_characteristics(SPECIES_CHARS_FILE)
    @test "thermal_lower" in names(df)
    @test "thermal_upper" in names(df)
    st = df[lowercase.(string.(df.SP)) .== "st", :][1, :]
    @test st.thermal_lower == 4.0
    @test st.thermal_upper == 20.0
    @test st.thermal_optimum == 12.0
end

@testset "optimum_sweep_optima" begin
    chars = DataFrame(SP=["St", "Ab"], thermal_lower=[4.0, 8.0], thermal_upper=[20.0, 30.0])
    cold = Guadex.optimum_sweep_optima(chars, ["ST", "AB"]; margin=2.0, fraction=0.0)
    mid = Guadex.optimum_sweep_optima(chars, ["ST", "AB"]; margin=2.0, fraction=0.5)
    warm = Guadex.optimum_sweep_optima(chars, ["ST", "AB"]; margin=2.0, fraction=1.0)
    @test cold ≈ [6.0, 10.0]
    @test mid ≈ [12.0, 19.0]
    @test warm ≈ [18.0, 28.0]
    @test_throws Exception Guadex.optimum_sweep_optima(chars, ["ST"]; fraction=2.0)
end

# =============================================================================
# 2. Heat-stress calibration and exposure (WP3)
# =============================================================================
@testset "calibrate_heat_stress_rate" begin
    k = Guadex.calibrate_heat_stress_rate([1.0, 4.0]; max_annual_loss=0.05)
    @test k ≈ -log(0.95) / 4.0
    @test Guadex.calibrate_heat_stress_rate([0.0, 0.0]) == 0.0
    @test_throws Exception Guadex.calibrate_heat_stress_rate([1.0]; max_annual_loss=1.5)
end

@testset "exceedance_energy and exposure_days" begin
    temps = [28.0, 29.0, 30.0, 20.0]
    @test Guadex.exceedance_energy(temps, 28.0; days_per_year=4.0) ≈ 5.0
    @test Guadex.exceedance_energy(temps, 100.0) == 0.0
    @test Guadex.exposure_days(temps, 28.0) == 2
    @test Guadex.exposure_days(temps, Inf) == 0
end

@testset "exposure_table" begin
    sites = ["a", "b"]
    species = ["ST", "AB"]
    dates = [Date(2026, 1, 1), Date(2026, 1, 2), Date(2027, 1, 1)]
    temps = [30.0 20.0 30.0; 10.0 10.0 10.0]
    tbl = Guadex.exposure_table(sites, species, temps, dates;
        upper_limits=[20.0, 30.0])
    @test nrow(tbl) == 2 * 2 * 2
    row = tbl[(tbl.CODIGO .== "a") .& (tbl.species .== "ST") .& (tbl.year .== 2026), :][1, :]
    @test row.exposure_days == 1
    @test row.exceedance_energy ≈ 100.0
end

# =============================================================================
# 2c. C7 — realised applied warming anomaly
# =============================================================================
@testset "realised warming anomaly (C7)" begin
    # Site × year annual-mean anomalies: sites (1.0, 2.0) in 2026, (3.0, 4.0)
    # in 2027.  The window means are the site-and-window means.
    warming = [1.0 3.0; 2.0 4.0]
    years = [2026, 2027]
    @test Guadex.mean_warming_anomaly(warming, years;
        first_year=2026, last_year=2027) ≈ mean([1.0, 2.0, 3.0, 4.0])
    @test Guadex.mean_warming_anomaly(warming, years;
        first_year=2027, last_year=2027) ≈ 3.5

    # The schedule overload recovers the same window mean from a daily series.
    dates = vcat(fill(Date(2026, 1, 1), 2), fill(Date(2027, 1, 1), 2))
    deltas = [1.0 1.0 3.0 3.0; 2.0 2.0 4.0 4.0]
    schedule = Guadex.TemperatureSchedule(deltas, 1.0)
    summary = Guadex.realised_warming_anomaly(schedule, dates;
        early_window=(2026, 2027), late_window=(2027, 2027))
    @test summary.mean_2026_2045 ≈ 2.5
    @test summary.mean_2036_2045 ≈ 3.5
    @test summary.n_sites == 2
    @test summary.n_years == 2
end

# =============================================================================
# 2d. C7 — exposure restricted to baseline-established sites
# =============================================================================
@testset "established_exposure_summary (C7)" begin
    sites = ["s1", "s2", "s3"]
    species = ["ST"]
    dates = [Date(2026, 1, 1), Date(2026, 1, 2), Date(2026, 1, 3)]
    # s2 is the hottest site (3 days above 20 °C) but is NOT baseline
    # established for ST; s1 has 1 day and s3 has 0.
    temps = [21.0 19.0 19.0;
             40.0 40.0 40.0;
             19.0 19.0 19.0]
    exposure = Guadex.exposure_table(sites, species, temps, dates; upper_limits=[20.0])
    all_sites_max = maximum(exposure.exposure_days)

    species_metrics = DataFrame(
        year=fill(2026, 3), CODIGO=sites, species=fill("ST", 3),
        baseline_established=[true, false, true])
    summary = Guadex.established_exposure_summary(exposure, species_metrics)
    row = summary[1, :]
    @test row.site_set == "baseline_established"
    @test row.n_sites_total == 3
    @test row.n_sites_selected == 2
    @test row.selected_sites == "s1;s3"
    @test row.max_exposure_days == 1.0
    @test row.mean_exposure_days ≈ 0.5
    # The established-site value differs from the all-sites value by design.
    @test row.max_exposure_days != all_sites_max
    @test all_sites_max == 3.0

    # No established information → explicit all-sites fallback.
    fallback = Guadex.established_exposure_summary(exposure, nothing)
    @test fallback[1, :site_set] == "all_sites"
    @test fallback[1, :n_sites_selected] == 3
    @test fallback[1, :max_exposure_days] == 3.0
end

# =============================================================================
# 2b. Daily forcing file loaders (WP1)
# =============================================================================
@testset "daily forcing loaders" begin
    df = DataFrame(
        site_id=["a", "a", "b", "b"],
        scenario=fill("ssp126", 4),
        date=[Date(2026, 1, 1), Date(2026, 1, 2), Date(2026, 1, 1), Date(2026, 1, 2)],
        tw_ensemble_median=[10.0, 11.0, 12.0, 13.0],
    )
    path = joinpath(mktempdir(), "forcing.csv")
    CSV.write(path, df)
    loaded = Guadex.load_daily_temperature_forcing(path)
    @test loaded.date isa Vector{Date}
    dates, temps = Guadex.daily_forcing_matrix(loaded, ["a", "b"]; scenario="ssp126")
    @test dates == [Date(2026, 1, 1), Date(2026, 1, 2)]
    @test temps == [10.0 11.0; 12.0 13.0]

    @test_throws Exception Guadex.load_daily_temperature_forcing(
        joinpath(mktempdir(), "missing.csv"))

    # Wide layout (WP1 primary product for the ODE).
    wide = DataFrame(
        date=[Date(2026, 1, 1), Date(2026, 1, 2)],
        scenario=fill("ssp126", 2),
        a=[10.0, 11.0], b=[12.0, 13.0],
    )
    @test Guadex.is_wide_daily_forcing(wide)
    @test !Guadex.is_wide_daily_forcing(df)
    wdates, wtemps = Guadex.wide_forcing_matrix(wide, ["a", "b"]; scenario="ssp126")
    @test wdates == [Date(2026, 1, 1), Date(2026, 1, 2)]
    @test wtemps == [10.0 11.0; 12.0 13.0]
    @test Guadex.wide_baseline_means(wide, ["a", "b"]; baseline_scenario="ssp126") ≈ [10.5, 12.5]
end

# =============================================================================
# 3. Daily / annual-mean schedule builders (WP2)
# =============================================================================
@testset "daily_temperature_schedule" begin
    dates = [Date(1986, 1, 1), Date(1986, 1, 2), Date(1986, 6, 1), Date(1986, 6, 2)]
    temps = [10.0 12.0 14.0 16.0]
    schedule, baseline = Guadex.daily_temperature_schedule(temps;
        dates=dates, baseline_start=1986, baseline_end=2005)
    @test baseline == [13.0]
    @test schedule.days_per_year == 1.0
    @test vec(schedule.deltas) ≈ [-3.0, -1.0, 1.0, 3.0]

    # Separate baseline series (WP1 'historical' rows).
    base_dates = [Date(2000, 1, 1), Date(2000, 1, 2)]
    base_temps = [8.0 12.0]
    s2, b2 = Guadex.daily_temperature_schedule(temps; dates=dates,
        baseline_start=1986, baseline_end=2005,
        baseline_temps=base_temps, baseline_dates=base_dates)
    @test b2 == [10.0]
    @test vec(s2.deltas) ≈ [0.0, 2.0, 4.0, 6.0]
end

@testset "baseline_climatology_schedule" begin
    dates = [Date(2000, 1, 1), Date(2000, 1, 2), Date(2000, 1, 3), Date(2000, 1, 4)]
    temps = [10.0 12.0 14.0 16.0]
    schedule = Guadex.baseline_climatology_schedule(temps, dates, 2; days_per_year=4.0)
    @test size(schedule.deltas) == (1, 8)
    @test vec(schedule.deltas)[1:4] ≈ [-3.0, -1.0, 1.0, 3.0]
    @test vec(schedule.deltas)[5:8] ≈ [-3.0, -1.0, 1.0, 3.0]
end

@testset "annual_mean_deltas_by_year" begin
    dates = [Date(2026, 1, 1), Date(2026, 7, 1), Date(2027, 1, 1), Date(2027, 7, 1)]
    schedule = Guadex.TemperatureSchedule([0.0 2.0 4.0 6.0], 1.0)
    years, mat = Guadex.annual_mean_deltas_by_year(schedule, dates)
    @test years == [2026, 2027]
    @test mat == [1.0 5.0]
end

# =============================================================================
# 4. Spin-up and carrying-capacity helpers (WP4)
# =============================================================================
@testset "spin_up converges to carrying capacity" begin
    n_sites, n_species = 1, 1
    params = Guadex.MetacommunityParams(
        1, 1, [0.0;;], sparse([1], [1], [0.0], 1, 1),
        [1.0], [0.05;;], [20.0], [1.0], [20.0], [5.0], [100.0])
    result = Guadex.spin_up(params; initial_state=[1.0], max_years=60, tol=1e-8)
    @test result.converged
    @test result.site_biomass[1] > 90.0
    @test result.site_biomass[1] <= 100.0 + 1e-6
    @test length(result.history) == result.years
end

@testset "scale_carrying_capacity / set_thermal_optima / set_heat_stress_rate" begin
    params = Guadex.MetacommunityParams(
        2, 2, [0.0 0.0; 0.0 0.0], sparse([1, 2], [1, 2], [0.0, 0.0], 2, 2),
        [1.0, 1.0], [0.1 0.1; 0.1 0.1], [20.0, 20.0], [1.0, 1.0],
        [15.0, 15.0], [3.0, 3.0], [50.0, 50.0])

    scaled = Guadex.scale_carrying_capacity(params, 3.0)
    @test scaled.carrying_capacity == [150.0, 150.0]
    @test scaled.heat_stress_rate == 0.0

    shifted = Guadex.set_thermal_optima(params, [10.0, 20.0])
    @test shifted.thermal_optima == [10.0, 20.0]

    stressed = Guadex.set_heat_stress_rate(shifted, 0.02)
    @test stressed.heat_stress_rate == 0.02
    @test stressed.thermal_optima == [10.0, 20.0]
end

# =============================================================================
# 4b. C5/E1: one corrected water-temperature series for the whole model
# =============================================================================
const WT_BASELINE_FILE = joinpath(@__DIR__, "..", "guadex_tw", "outputs", "tables",
    "water_temp_baseline_guadex_sites.csv")
const WT_WIDE_FILE = joinpath(@__DIR__, "..", "guadex_tw", "outputs", "tables",
    "water_temp_daily_guadex_sites_wide.csv")

@testset "corrected site water-temperature baseline (C5/E1)" begin
    @test isfile(WT_BASELINE_FILE)
    lookup = Guadex.load_site_water_temperature_baseline(WT_BASELINE_FILE)
    @test length(lookup) == 776
    @test haskey(lookup, "1.30.20")   # extra baseline row, excluded from the model

    # Exact-column selection: a table without `tw_baseline_mean` must error, so a
    # wrong column (e.g. TEMP_MEDIA_SC) can never be selected silently.
    mktempdir() do tmp
        wrong = joinpath(tmp, "wrong.csv")
        CSV.write(wrong, DataFrame(site_id=["a"], TEMP_MEDIA_SC=[11.0]))
        @test_throws Exception Guadex.load_site_water_temperature_baseline(wrong)
    end

    # Trout sites: sites with a non-zero baseline brown-trout density.  The
    # corrected mean is now ≈12.83 °C (pre-E22 pre-correction value 11.7 °C).
    dens = CSV.read(DENSITY_FILE, DataFrame)
    trout = [string(r.CODIGO) for r in eachrow(dens)
             if r.CODIGO in keys(lookup) && !ismissing(r.ST_DEN) && r.ST_DEN > 0.0]
    @test length(trout) == 36
    trout_mean = mean([lookup[s] for s in trout])
    @test isapprox(trout_mean, 12.830498685857533; atol=0.01)
end

@testset "model temperature equals corrected daily series (C5/E1)" begin
    lookup = Guadex.load_site_water_temperature_baseline(WT_BASELINE_FILE)
    site_df = Guadex.load_site_data(CONNECTIVITY_FILE, ENVIRONMENTAL_FILE)
    sites = String.(site_df.CODIGO)
    sample = sites[[1, 100, 400]]
    levels = [lookup[s] for s in sample]

    cols = vcat([:date, :scenario], Symbol.(sample))
    df = CSV.read(WT_WIDE_FILE, DataFrame; select=cols)
    df.date = Date.(df.date)

    fdates, ftemps = Guadex.wide_forcing_matrix(df, sample; scenario="ssp585")
    bdates, btemps = Guadex.wide_forcing_matrix(df, sample; scenario="historical")

    # The corrected baseline is exactly the 1986-2005 mean of the historical
    # daily series (no double correction, no offset).
    bmask = [1986 <= year(d) <= 2005 for d in bdates]
    for i in eachindex(sample)
        @test isapprox(mean(btemps[i, bmask]), levels[i]; atol=1e-6)
    end

    # Reference the anomaly to the SAME per-site mean used as the model level.
    sched, bmeans = Guadex.daily_temperature_schedule(ftemps;
        dates=fdates, baseline_start=1986, baseline_end=2005,
        baseline_means=levels)
    @test bmeans == levels

    # level + anomaly reconstructs the corrected daily series at sample days.
    for i in eachindex(sample), day in (1, 500, 3000, length(fdates))
        t = Float64(day - 1)   # the schedule stores one node per day here
        model_T = levels[i] + Guadex.temperature_delta(sched, i, t)
        @test isapprox(model_T, ftemps[i, day]; atol=1e-9)
    end
end

# =============================================================================
# 5. Quasi-extinction (WP5)
# =============================================================================
@testset "quasi_extinction_flags" begin
    # Baseline 10, q = 0.1 → cutoff 1.0 (above the 0.1 presence threshold).
    series = [10.0, 10.0, 0.5, 0.5, 0.5, 2.0, 0.5, 0.5, 0.5]
    flags = Guadex.quasi_extinction_flags(series; threshold=0.1,
        baseline_density=10.0, q=0.1, persistence=3)
    @test flags == [false, false, false, false, true, false, false, false, true]
    @test Guadex.time_to_quasi_extinction(flags; year_labels=2000:2008) == 2004

    # Presence threshold dominates when the baseline is small.
    flags2 = Guadex.quasi_extinction_flags([0.05, 0.05, 0.05];
        threshold=0.1, baseline_density=0.2, q=0.1, persistence=3)
    @test flags2 == [false, false, true]
    @test ismissing(Guadex.time_to_quasi_extinction([false, false];
        year_labels=[1, 2]))
end
