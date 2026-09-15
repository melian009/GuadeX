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
