using Test
using DataFrames
using CSV
using Statistics
using SparseArrays
using Guadex

# =============================================================================
# 1. species_indices
# =============================================================================
@testset "species_indices" begin
    species = ["AB", "GH", "SP", "MS"]
    @test Guadex.species_indices(species, ["AB", "SP"]) == [1, 3]
    @test Guadex.species_indices(species, ["gh", "MS"]) == [2, 4]
    @test Guadex.species_indices(species, ["XX"]) == Int[]
end

# =============================================================================
# 2. Site -> level vectors and crosswalk
# =============================================================================
@testset "site_level_vectors" begin
    crosswalk = DataFrame(
        CODIGO=["1.1.2", "1.1.3", "2.1.1"],
        CODIGO_S=["1.1", "1.1", "2.1"],
        ID_masa=["ES050MSPF001", "ES050MSPF001", "ES050MSPF002"],
        water_body_name=["Body A", "Body A", "Body B"],
        subzone=["1", "1", "2"], zone=["0001", "0001", "0002"], system=["S1", "S1", "S2"],
    )
    sites = ["1.1.2", "1.1.3", "2.1.1"]
    site_df = DataFrame(CODIGO=["1.1.2", "1.1.3", "2.1.1"], CODIGO_S=["1.1", "1.1", "2.1"])
    levels = Guadex.site_level_vectors(sites, site_df, crosswalk)

    @test levels.subcatchment == ["1.1", "1.1", "2.1"]
    @test levels.water_body == ["ES050MSPF001", "ES050MSPF001", "ES050MSPF002"]
    @test levels.water_body_name[3] == "Body B"
    @test all(levels.basin .== "ES050")

    # Unknown sites fall back to explicit sentinels.
    fallback = Guadex.site_level_vectors(["9.9.9"], nothing, crosswalk)
    @test fallback.water_body == [Guadex.GUADEX_UNASSIGNED_WATER_BODY]
    @test fallback.subcatchment == [Guadex.GUADEX_UNKNOWN_SUBCATCHMENT]
end

# =============================================================================
# 3. compute_site_metrics
# =============================================================================
function synthetic_timeseries(n_sites, n_species)
    times = [0.0, 182.5, 365.0, 547.5, 730.0]
    states = Vector{Vector{Float64}}()
    for (k, t) in enumerate(times)
        mat = zeros(n_sites, n_species)
        if k == 1 || k == 3
            # Site 1 native present at t = 0 and at the end of year 1 (t = 365),
            # matching the new end-of-year reporting offsets.
            mat[1, 1] = 5.0
            mat[2, 2] = 3.0   # site 2 invasive present
        else
            mat[2, 2] = 3.0
            k == 5 && (mat[1, 2] = 1.0)
        end
        push!(states, vec(mat))
    end
    return times, states
end

@testset "compute_site_metrics" begin
    sites = ["1.1.2", "1.1.3"]
    species = ["AB", "GH"]
    levels = (
        subcatchment=["1.1", "1.1"],
        water_body=["ES050MSPF001", "ES050MSPF001"],
        water_body_name=["Body A", "Body A"],
        zone=["0001", "0001"], subzone=["1", "1"], basin=["ES050", "ES050"],
    )
    times, states = synthetic_timeseries(2, 2)
    warming = [0.0 1.0; 0.0 2.0]

    df = Guadex.compute_site_metrics(times, states, sites, levels, species,
        ["AB"], ["GH"];
        days_per_year=365.0, year_offsets=[1, 2], year_labels=[2026, 2027],
        temperature_baseline=[15.0, 20.0], warming=warming, threshold=0.1)

    @test nrow(df) == 4
    @test sort(unique(df.year)) == [2026, 2027]

    r = df[(df.year .== 2026) .& (df.CODIGO .== "1.1.2"), :][1, :]
    @test r.native_richness == 1
    @test r.total_biomass == 5.0
    @test r.realised_richness_loss == 0.0
    @test r.temperature_c == 15.0

    r = df[(df.year .== 2027) .& (df.CODIGO .== "1.1.2"), :][1, :]
    @test r.native_richness == 0
    @test r.native_richness_relative == 0.0
    @test r.realised_richness_loss == 1.0
    @test r.temperature_c == 16.0
    @test r.delta_temperature_c == 1.0

    r2 = df[(df.year .== 2027) .& (df.CODIGO .== "1.1.3"), :][1, :]
    @test r2.invasive_richness == 1
    @test r2.realised_richness_loss == 0.0
    @test r2.temperature_c == 22.0
end

# =============================================================================
# 3b. Baseline temperature is reported when no warming matrix is supplied
# =============================================================================
@testset "temperature baseline without warming" begin
    sites = ["1.1.2", "1.1.3"]
    species = ["AB", "GH"]
    levels = (
        subcatchment=["1.1", "1.1"],
        water_body=["ES050MSPF001", "ES050MSPF001"],
        water_body_name=["Body A", "Body A"],
        zone=["0001", "0001"], subzone=["1", "1"], basin=["ES050", "ES050"],
    )
    times, states = synthetic_timeseries(2, 2)
    df = Guadex.compute_site_metrics(times, states, sites, levels, species,
        ["AB"], ["GH"];
        year_offsets=[1], year_labels=[2026],
        temperature_baseline=[15.0, 20.0], warming=nothing)
    @test all(isfinite.(df.temperature_c))
    @test df[df.CODIGO .== "1.1.2", :temperature_c] == [15.0]
    @test df[df.CODIGO .== "1.1.3", :temperature_c] == [20.0]
    @test all(df.delta_temperature_c .== 0.0)
end

# =============================================================================
# 3c. Directional connectivity is not transposed
# =============================================================================
@testset "site_connectivity_metrics directionality" begin
    sites = ["a", "b"]
    distances = sparse([1, 2], [2, 1], [1000.0, 1000.0], 2, 2)
    # movement 1 -> 2 uses dams[2, 1] = 0.2; movement 2 -> 1 uses dams[1, 2] = 0.5
    dams = [1.0 0.5; 0.2 1.0]
    conn = Guadex.site_connectivity_metrics(sites, dams, distances)
    @test conn.mean_in_passability ≈ [0.5, 0.2]
    @test conn.mean_out_passability ≈ [0.2, 0.5]
    @test conn.barrier_links_in == [1, 1]
    @test conn.barrier_links_out == [1, 1]
end

# =============================================================================
# 4. aggregate_metrics
# =============================================================================
@testset "aggregate_metrics" begin
    sites = ["1.1.2", "1.1.3", "2.1.1"]
    species = ["AB", "GH"]
    levels = (
        subcatchment=["1.1", "1.1", "2.1"],
        water_body=["ES050MSPF001", "ES050MSPF001", "ES050MSPF002"],
        water_body_name=["Body A", "Body A", "Body B"],
        zone=["0001", "0001", "0002"], subzone=["1", "1", "2"],
        basin=["ES050", "ES050", "ES050"],
    )
    times = [0.0, 365.0, 730.0]
    states = [vec([5.0 0.0; 3.0 0.0; 0.0 2.0]), vec([5.0 0.0; 3.0 0.0; 0.0 2.0]),
        vec([5.0 0.0; 3.0 0.0; 0.0 2.0])]
    df = Guadex.compute_site_metrics(times, states, sites, levels, species, ["AB"], ["GH"];
        year_offsets=[1, 2], year_labels=[2026, 2027])

    basin = Guadex.aggregate_metrics(df; level=:basin)
    @test nrow(basin) == 2
    @test all(basin.n_sites .== 3)
    row = basin[basin.year .== 2026, :][1, :]
    @test row.mean_native_richness ≈ (1 + 1 + 0) / 3
    @test row.mean_total_biomass ≈ (5 + 3 + 2) / 3

    sub = Guadex.aggregate_metrics(df; level=:subcatchment)
    @test nrow(sub) == 4  # two subcatchments × two years
    s11 = sub[(sub.subcatchment .== "1.1") .& (sub.year .== 2026), :][1, :]
    @test s11.n_sites == 2
    @test s11.mean_total_biomass ≈ (5 + 3) / 2

    wb = Guadex.aggregate_metrics(df; level=:water_body)
    @test "water_body_name" in names(wb)
    @test nrow(wb) == 4
end

# =============================================================================
# 5. Temperature schedule interpolation
# =============================================================================
@testset "TemperatureSchedule" begin
    schedule = Guadex.TemperatureSchedule([0.0 1.0 2.0; 0.0 2.0 4.0], 365.0)
    @test Guadex.temperature_delta(schedule, 1, 0.0) ≈ 0.0
    @test Guadex.temperature_delta(schedule, 1, 365.0) ≈ 1.0
    @test Guadex.temperature_delta(schedule, 1, 547.5) ≈ 1.5
    @test Guadex.temperature_delta(schedule, 2, 365.0) ≈ 2.0
    # Beyond the last node → held constant.
    @test Guadex.temperature_delta(schedule, 2, 5000.0) ≈ 4.0
end

# =============================================================================
# 6. basin_warming_curve from projection-shaped rows
# =============================================================================
@testset "basin_warming_curve" begin
    projections = DataFrame(
        site_id=["A", "B", "A", "B", "A", "B"],
        scenario=fill("ssp126", 6),
        gcm=fill("GCM1", 6),
        period_start=[2021, 2021, 2036, 2036, 2041, 2041],
        period_end=[2040, 2040, 2055, 2055, 2070, 2070],
        delta_tw_mean=[0.4, 0.6, 0.8, 1.0, 1.2, 1.4],
    )
    years, curve = Guadex.basin_warming_curve(projections, "ssp126", "GCM1", 2026, 2045)
    @test years[1] == 2026
    @test curve[1] == 0.0
    @test all(diff(curve) .>= -1e-12)
    @test curve[end] > 0.5
    @test issorted(years)
end

# =============================================================================
# 7. export_run_outputs end-to-end (files + viewer content)
# =============================================================================
@testset "export_run_outputs" begin
    sites = ["1.1.2", "1.1.3", "2.1.1"]
    species = ["AB", "GH"]
    site_df = DataFrame(CODIGO=sites, CODIGO_S=["1.1", "1.1", "2.1"])
    times, states = synthetic_timeseries(3, 2)

    dir = mktempdir()
    result = Guadex.export_run_outputs(dir;
        sol_t=times, sol_u=states, sites=sites, species=species, site_df=site_df,
        native_species=["AB"], invasive_species=["GH"],
        report_year_offsets=[1, 2], report_year_labels=[2026, 2027],
        crosswalk_path=joinpath(@__DIR__, "..", "data", "site_waterbody_crosswalk.csv"),
        run_metadata=Dict("scenario" => "test"))

    for rel in ["levels/level_sampling_point.csv", "levels/level_subcatchment.csv",
                "levels/level_water_body.csv", "levels/level_basin.csv",
                "viewer/guadex_results_realised_richness_loss_timeseries.csv",
                "viewer/guadex_results_metrics.json",
                "viewer/guadex_results_realised_richness_loss.json",
                # minor #14: species-level series and a renderable aggregate
                "viewer/species_AB_density_timeseries.json",
                "viewer/species_GH_relative_density_timeseries.json",
                "viewer/level_basin_mean_native_richness_timeseries.json",
                "viewer/level_basin_mean_native_richness_by_site_timeseries.json",
                "run_metadata.json"]
        @test isfile(joinpath(dir, rel))
        @test filesize(joinpath(dir, rel)) > 0
    end

    # A per-metric CSV is a single-variable time series (CODIGO, step, value),
    # never the old 15-metric × 20-year mixed slider.
    csv = CSV.read(joinpath(dir, "viewer",
        "guadex_results_realised_richness_loss_timeseries.csv"), DataFrame)
    @test names(csv) == ["CODIGO", "step", "value"]
    @test Set(csv.CODIGO) == Set(sites)

    # Species-level and by-site aggregate files use the viewer's timeseries
    # shape, keyed by CODIGO.
    species_ts = read(joinpath(dir, "viewer", "species_AB_density_timeseries.json"), String)
    @test occursin("\"steps\"", species_ts)
    @test occursin("\"1.1.2\"", species_ts)
    agg_ts = read(joinpath(dir, "viewer",
        "level_basin_mean_native_richness_by_site_timeseries.json"), String)
    @test occursin("\"steps\"", agg_ts)
    @test occursin("\"1.1.2\"", agg_ts)

    meta = read(joinpath(dir, "run_metadata.json"), String)
    @test occursin("sampling_point", meta)
    @test occursin("\"scenario\": \"test\"", meta)
    # A2: the old `crosswalk_water_bodies` counted sites, not water bodies.
    # Two of the three synthetic sites resolve to distinct water bodies; the
    # third is unassigned, so water-body levels = 2 assigned + 1 unassigned.
    @test !occursin("crosswalk_water_bodies", meta)
    @test occursin("\"sites_with_assigned_water_body\": 2", meta)
    @test occursin("\"n_assigned_water_bodies\": 2", meta)
    @test occursin("\"n_water_bodies\": 3", meta)

    @test nrow(result.site_metrics) == 6
end

# =============================================================================
# 8. Species-level abundance, occupancy and quasi-extinction (WP5)
# =============================================================================
@testset "compute_species_metrics and quasi-extinction" begin
    dates = 2026:2030
    times = collect(0.0:365.0:1825.0)
    states = [vec([d;;]) for d in [10.0, 10.0, 10.0, 0.5, 0.5, 0.5]]

    sp = Guadex.compute_species_metrics(times, states, ["s1"], ["ST"];
        year_offsets=[1, 2, 3, 4, 5], year_labels=dates,
        baseline_species_density=[10.0;;],
        quasi_extinction_q=0.1, quasi_extinction_persistence=3)

    @test nrow(sp) == 5
    @test sp[sp.year .== 2028, :quasi_extinct] == [false]
    @test sp[sp.year .== 2030, :quasi_extinct] == [true]
    @test sp[1, :baseline_density] == 10.0
    @test sp[1, :baseline_established] == true
    @test all(sp.baseline_established)
    @test sp[1, :relative_density] == 1.0
    @test sp[1, :time_to_quasi_extinction] == 2030

    summary = Guadex.quasi_extinction_summary(sp)
    row = summary[1, :]
    @test row.species == "ST"
    # Occupancy (presence threshold 0.1) stays 1: that is exactly why WP5 adds
    # the stricter quasi-extinction metric on top of it.
    @test row.occupancy_final == 1.0
    @test row.n_sites == 1
    @test row.n_sites_established == 1
    @test row.quasi_extinct_fraction == 1.0
    @test row.median_time_to_quasi_extinction == 2030.0
    @test row.relative_biomass_change ≈ 0.5 / 10.0 - 1.0
end

@testset "compute_site_metrics rebased on spun-up baseline" begin
    species = ["ST", "GH"]
    levels = (
        subcatchment=["1.1"], water_body=["ES050MSPF001"],
        water_body_name=["Body A"], zone=["0001"], subzone=["1"], basin=["ES050"],
    )
    times = collect(0.0:365.0:1095.0)
    states = [vec([10.0 5.0]), vec([10.0 5.0]), vec([0.5 5.0]), vec([0.5 5.0])]
    baseline = [10.0 5.0]

    df = Guadex.compute_site_metrics(times, states, ["s1"], levels, species, ["ST"], ["GH"];
        year_offsets=[1, 2, 3], year_labels=[2026, 2027, 2028],
        baseline_species_density=baseline, quasi_extinction_q=0.1,
        quasi_extinction_persistence=2)

    r = df[(df.year .== 2027) .& (df.CODIGO .== "s1"), :][1, :]
    @test r.native_biomass_relative ≈ 0.05
    @test r.native_occupancy ≈ 1.0
    @test r.invasive_occupancy ≈ 1.0
    @test r.native_quasi_extinct == 0
    @test r.native_baseline_richness == 1

    r2 = df[(df.year .== 2028) .& (df.CODIGO .== "s1"), :][1, :]
    @test r2.native_quasi_extinct == 1
    @test r2.native_quasi_extinct_fraction ≈ 1.0
end

# =============================================================================
# 8b. E6 — quasi-extinction only where the species was established
# =============================================================================
@testset "quasi_extinction_flags never flags non-established sites (E6)" begin
    # Baseline density at or below the presence threshold → never flagged.
    @test Guadex.quasi_extinction_flags(zeros(5);
        threshold=0.1, baseline_density=0.0, persistence=3) == falses(5)
    @test Guadex.quasi_extinction_flags(zeros(5);
        threshold=0.1, baseline_density=0.05, persistence=3) == falses(5)
    # An explicit mask wins over the baseline value.
    @test Guadex.quasi_extinction_flags(zeros(5);
        threshold=0.1, baseline_density=10.0, persistence=3,
        established=false) == falses(5)
    # Established sites keep the original behaviour (three consecutive low years
    # after two low-density free years → flag from the third low year).
    @test Guadex.quasi_extinction_flags([10.0, 10.0, 0.0, 0.0, 0.0];
        threshold=0.1, baseline_density=10.0, persistence=3) ==
        [false, false, false, false, true]
end

@testset "quasi-extinction fraction over established sites (E6)" begin
    # 3 sites, one species: established at sites 1 and 3 (density 10), absent
    # at site 2 (baseline 0).  Site 1 collapses, sites 2 and 3 do not.
    times = collect(0.0:365.0:1825.0)
    baseline = reshape([10.0, 0.0, 10.0], 3, 1)
    series = (
        [10.0, 10.0, 0.0, 0.0, 0.0, 0.0],  # site 1: collapses below cutoff 1.0
        zeros(6),                          # site 2: never occupied
        fill(10.0, 6),                     # site 3: stable
    )
    states = [vec([series[1][k], series[2][k], series[3][k]]) for k in 1:6]

    sp = Guadex.compute_species_metrics(times, states, ["s1", "s2", "s3"], ["ST"];
        year_offsets=[1, 2, 3, 4, 5], year_labels=2026:2030,
        baseline_species_density=baseline, quasi_extinction_persistence=3)

    @test all(sp[sp.CODIGO .== "s1", :baseline_established])
    @test all(.!sp[sp.CODIGO .== "s2", :baseline_established])
    final = sp[sp.year .== 2030, :]
    @test final[final.CODIGO .== "s1", :quasi_extinct] == [true]
    @test final[final.CODIGO .== "s2", :quasi_extinct] == [false]
    @test final[final.CODIGO .== "s3", :quasi_extinct] == [false]

    row = Guadex.quasi_extinction_summary(sp)[1, :]
    @test row.n_sites == 3
    @test row.n_sites_established == 2
    # 1 of the 2 established sites went quasi-extinct; the never-occupied site
    # must not inflate the numerator or the denominator.
    @test row.quasi_extinct_fraction ≈ 1 / 2
end

@testset "species absent everywhere in the baseline → fraction 0 (E6)" begin
    times = collect(0.0:365.0:1825.0)
    states = [vec(zeros(2, 1)) for _ in times]
    sp = Guadex.compute_species_metrics(times, states, ["s1", "s2"], ["ST"];
        year_offsets=[1, 2, 3, 4, 5], year_labels=2026:2030,
        baseline_species_density=zeros(2, 1), quasi_extinction_persistence=3)

    @test all(.!sp.baseline_established)
    @test all(.!sp.quasi_extinct)
    row = Guadex.quasi_extinction_summary(sp)[1, :]
    @test row.n_sites == 2
    @test row.n_sites_established == 0
    # No established site → documented convention is 0.0, not NaN and not 1.0.
    @test row.quasi_extinct_fraction == 0.0
    @test !isnan(row.quasi_extinct_fraction)
end

@testset "trout-style fraction is not 1 - occupancy (E6)" begin
    # Established at 3 of 4 sites and nowhere collapsing: the old code flagged
    # the never-occupied 4th site from year 3 and returned 1 - occupancy = 0.25.
    times = collect(0.0:365.0:1825.0)
    baseline = reshape([10.0, 10.0, 10.0, 0.0], 4, 1)
    states = [vec(fill(10.0, 4)) for _ in times]
    for st in states
        st[4] = 0.0
    end

    sp = Guadex.compute_species_metrics(times, states,
        ["s1", "s2", "s3", "s4"], ["ST"];
        year_offsets=[1, 2, 3, 4, 5], year_labels=2026:2030,
        baseline_species_density=baseline, quasi_extinction_persistence=3)

    row = Guadex.quasi_extinction_summary(sp)[1, :]
    @test row.occupancy_final ≈ 0.75
    @test 1.0 - row.occupancy_final ≈ 0.25
    @test row.n_sites_established == 3
    @test row.quasi_extinct_fraction == 0.0
    @test row.quasi_extinct_fraction != 1.0 - row.occupancy_final
end

# =============================================================================
# 9. Annual snapshot alignment (E5): end-of-year rows, full horizon, forcing
# =============================================================================
@testset "annual snapshot alignment" begin
    days_per_year = 365.0
    n_years = 20
    t_end = n_years * days_per_year
    interval = days_per_year / 12.0
    times = unique(vcat(collect(0.0:interval:t_end), t_end))
    # Density equals elapsed time, so the initial and end-of-horizon states are
    # distinguishable and each snapshot can be attributed to its exact time.
    states = [vec([Float64(t)]) for t in times]

    sites = ["s1"]
    species = ["ST"]
    levels = (
        subcatchment=["1.1"], water_body=["ES050MSPF001"],
        water_body_name=["Body A"], zone=["0001"], subzone=["1"], basin=["ES050"],
    )
    offsets = collect(1:n_years)
    labels = collect(2026:2045)

    df = Guadex.compute_site_metrics(times, states, sites, levels, species,
        ["ST"], String[];
        days_per_year=days_per_year, year_offsets=offsets, year_labels=labels,
        baseline_species_density=[10.0;;])

    # (a) The first reported row is the end of year 2026 (t = 365), not the
    # pre-forcing initial state (t = 0).
    @test df[df.year .== 2026, :total_biomass][1] == 365.0
    @test df[df.year .== 2026, :total_biomass][1] != states[1][1]

    # (b) The last row is the state at the end of the simulation horizon.
    @test df[df.year .== 2045, :total_biomass][1] == t_end
    @test df[end, :year] == 2045

    # The snapshot index resolves to exactly the intended saveat grid point.
    for offset in offsets
        @test times[Guadex._snapshot_index(times, offset * days_per_year)] ==
            offset * days_per_year
    end

    # ITEM 2: the shared helper the standalone exporter now uses derives the same
    # 1-based end-of-year grid (offset k = state at k * days_per_year), so a
    # re-export cannot fall back to the old 0-based grid.
    @test Guadex._report_offsets(times, days_per_year) == collect(1:n_years)
    @test first(Guadex._report_offsets(times, days_per_year)) == 1

    # A grid point an ulp below the target must not shift the snapshot onto the
    # following month.
    perturbed = [0.0, 365.0 - 1.0e-9, 730.0]
    @test Guadex._snapshot_index(perturbed, 365.0) == 2

    # The row labelled Y reports year Y's warming, and 2045 maps to the final
    # forcing node (the schedule is held there for the last year).
    warming = reshape(Float64.(1:n_years), 1, n_years)
    dfw = Guadex.compute_site_metrics(times, states, sites, levels, species,
        ["ST"], String[];
        days_per_year=days_per_year, year_offsets=offsets, year_labels=labels,
        temperature_baseline=[15.0], warming=warming,
        baseline_species_density=[10.0;;])
    @test dfw[dfw.year .== 2026, :delta_temperature_c][1] == 1.0
    @test dfw[dfw.year .== 2045, :delta_temperature_c][1] == 20.0
    @test dfw[dfw.year .== 2045, :temperature_c][1] == 35.0
end

@testset "JSON number serialisation" begin
    @test Guadex._json_number(3.0) == "3"
    @test Guadex._json_number(3.5) == "3.5"
    @test Guadex._json_number(Inf) == "null"
    @test Guadex._json_number(NaN) == "null"
    # A relative metric divided by a near-zero burn-in baseline can be huge;
    # serialisation must emit a JSON float rather than throw InexactError.
    @test Guadex._json_number(5.3e164) == string(5.3e164)
    @test Guadex._json_number(2.0^60) == string(2.0^60)
end
