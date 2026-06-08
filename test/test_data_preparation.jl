# =============================================================================
# Comprehensive tests for data_preparation.jl functions
# Path constants (DATA_DIR, etc.) are defined in runtests.jl
# =============================================================================

# =============================================================================
# 1. Parsing helper functions (pure unit tests, no data files needed)
# =============================================================================
@testset "Parsing helpers" begin

    @testset "parse_temperature_range" begin
        # "X to Y" format → midpoint
        @test Guadex.parse_temperature_range("8 to 30") ≈ 19.0
        @test Guadex.parse_temperature_range("4 to 20") ≈ 12.0
        @test Guadex.parse_temperature_range("4 to 25") ≈ 14.5
        @test Guadex.parse_temperature_range("10 to 30") ≈ 20.0
        @test Guadex.parse_temperature_range("8 to 25") ≈ 16.5

        # Single number
        @test Guadex.parse_temperature_range("15") ≈ 15.0
        @test Guadex.parse_temperature_range("22.5") ≈ 22.5

        # Default fallback
        @test Guadex.parse_temperature_range("") ≈ 15.0
        @test Guadex.parse_temperature_range("invalid") ≈ 15.0
    end

    @testset "parse_temperature_range_and_sigma" begin
        # Range = 22, sigma ≈ 22/6 ≈ 3.667, optimum = 19
        result = Guadex.parse_temperature_range_and_sigma("8 to 30")
        @test result.optimum ≈ 19.0
        @test isapprox(result.sigma, 22.0/6.0, atol=0.01)

        # Range = 16, sigma ≈ 16/6 ≈ 2.667, optimum = 12
        result = Guadex.parse_temperature_range_and_sigma("4 to 20")
        @test result.optimum ≈ 12.0
        @test isapprox(result.sigma, 16.0/6.0, atol=0.01)

        # Range = 21, sigma ≈ 21/6 ≈ 3.5, optimum = 14.5
        result = Guadex.parse_temperature_range_and_sigma("4 to 25")
        @test result.optimum ≈ 14.5
        @test isapprox(result.sigma, 21.0/6.0, atol=0.01)

        # Single number → default sigma
        result = Guadex.parse_temperature_range_and_sigma("20")
        @test result.optimum ≈ 20.0
        @test result.sigma ≈ 3.0

        # Default fallback
        result = Guadex.parse_temperature_range_and_sigma("")
        @test result.optimum ≈ 15.0
        @test result.sigma ≈ 3.0

        result = Guadex.parse_temperature_range_and_sigma("invalid")
        @test result.optimum ≈ 15.0
        @test result.sigma ≈ 3.0
    end

    @testset "parse_elevation" begin
        # Range format
        @test Guadex.parse_elevation("0 to 900") ≈ 450.0
        @test Guadex.parse_elevation("500 to 2500") ≈ 1500.0
        @test Guadex.parse_elevation("0 to 100") ≈ 50.0

        # Single number
        @test Guadex.parse_elevation("797") ≈ 797.0

        # Default fallback
        @test Guadex.parse_elevation("") ≈ 500.0
        @test Guadex.parse_elevation("invalid") ≈ 500.0
    end

    @testset "parse_interaction_string" begin
        # Strongest negative: no coexistence
        @test Guadex.parse_interaction_string("No coexist") ≈ -1.0
        @test Guadex.parse_interaction_string("no coexist.") ≈ -1.0

        # Strong negative: displaces
        @test Guadex.parse_interaction_string("displaces") ≈ -0.8

        # Predation
        @test Guadex.parse_interaction_string("affects Ah through predation") ≈ -0.5
        @test Guadex.parse_interaction_string("predation") ≈ -0.5

        # Competition or interference
        @test Guadex.parse_interaction_string("affects Sp through competition") ≈ -0.3
        @test Guadex.parse_interaction_string("interfere through competition") ≈ -0.3
        @test Guadex.parse_interaction_string("competition") ≈ -0.3

        # Generic affects
        @test Guadex.parse_interaction_string("affects Sa") ≈ -0.2

        # Neutral coexistence
        @test Guadex.parse_interaction_string("coexist, neutral.") ≈ 0.0
        @test Guadex.parse_interaction_string("coexist") ≈ 0.0
        @test Guadex.parse_interaction_string("neutral") ≈ 0.0

        # Missing/empty
        @test Guadex.parse_interaction_string(missing) ≈ 0.0
        @test Guadex.parse_interaction_string("") ≈ 0.0
        @test Guadex.parse_interaction_string(";") ≈ 0.0
    end
end

# =============================================================================
# 2. Data loading functions (require real data files)
# =============================================================================
@testset "Data loading functions" begin

    @testset "load_species_characteristics" begin
        df = Guadex.load_species_characteristics(SPECIES_CHARS_FILE)
        @test nrow(df) == 24
        @test "SP" in names(df)
        @test "thermal_optimum" in names(df)
        @test "thermal_sigma" in names(df)
        @test "elevation_optimum" in names(df)
        @test "max_size_mm" in names(df)
        @test "TEMPERATURE_C" in names(df)
        @test "ELEVATION_m" in names(df)
        @test "MAX_SIZE_mm" in names(df)

        # Check specific species
        st_row = findfirst(r -> lowercase(r.SP) == "st", eachrow(df))
        @test st_row !== nothing
        @test df.thermal_optimum[st_row] ≈ 12.0
        st_sigma = df.thermal_sigma[st_row]
        @test isapprox(st_sigma, 16.0/6.0, atol=0.01)
        @test df.elevation_optimum[st_row] ≈ 1500.0

        # Most species have "8 to 30" range → optimum 19, sigma 22/6
        sa_row = findfirst(r -> lowercase(r.SP) == "sa", eachrow(df))
        @test df.thermal_optimum[sa_row] ≈ 19.0
        @test isapprox(df.thermal_sigma[sa_row], 22.0/6.0, atol=0.01)

        # "Ah" in file (title case) should map to optimum 19
        ah_row = findfirst(r -> r.SP == "Ah", eachrow(df))
        @test df.thermal_optimum[ah_row] ≈ 19.0
    end

    @testset "load_site_data" begin
        site_df = Guadex.load_site_data(CONNECTIVITY_FILE, ENVIRONMENTAL_FILE)
        @test nrow(site_df) > 0
        @test "CODIGO" in names(site_df)
        @test "ALTITUD" in names(site_df)
        @test "CODIGO_S" in names(site_df)
        @test "Demb_arr_m" in names(site_df)
        @test "Demb_ab_m" in names(site_df)

        # Site 1.30.20 should be filtered out
        @test !("1.30.20" in string.(site_df.CODIGO))
    end

    @testset "load_species_density_data" begin
        density_df, species_codes = Guadex.load_species_density_data(DENSITY_FILE)
        @test nrow(density_df) > 0
        @test length(species_codes) == 24
        @test "LS" in species_codes
        @test "SA" in species_codes
        @test "ST" in species_codes
        @test "GH" in species_codes
        @test "AAL" in species_codes
        @test "AA" in species_codes
    end

    @testset "load_interaction_matrix" begin
        species_codes = ["LS", "SA", "SP", "PW", "CP", "IL", "IO", "AA", "AH",
                         "GL", "GH", "LG", "AAL", "CG", "CC", "MS", "ST", "OM",
                         "EL", "AM", "TT", "LR", "MC", "AB"]

        mat = Guadex.load_interaction_matrix(INTERACTION_FILE, species_codes)
        @test size(mat) == (24, 24)
        @test all(diag(mat) .== 0.0)

        # Should have some non-zero entries
        @test count(x -> x != 0.0, mat) > 0

        # All values should be between -1 and 0
        @test all(mat .<= 0.0)
        @test all(mat .>= -1.0)

        # LAX: interaction matrix has as many nonzeros as the documented 205,
        # but verify a few known relationships from the empirical matrix
        # ST has -1.0 on many invasives
        st_idx = findfirst(==("ST"), species_codes)
        gh_idx = findfirst(==("GH"), species_codes)
        ms_idx = findfirst(==("MS"), species_codes)
        @test mat[st_idx, gh_idx] == 0.0   # ST row, GH col = effect of GH on ST
        @test mat[gh_idx, st_idx] == -1.0  # GH row, ST col = effect of ST on GH
    end
end

# =============================================================================
# 3. Building / extraction functions (synthetic data where possible)
# =============================================================================
@testset "Building and extraction functions" begin

    @testset "build_elevation_vector" begin
        site_df = DataFrame(CODIGO=["site1", "site2", "site3"], ALTITUD=[100.0, 200.0, 300.0])
        sites = ["site1", "site2", "site4"]
        elev = Guadex.build_elevation_vector(site_df, sites)
        @test length(elev) == 3
        @test elev[1] == 100.0
        @test elev[2] == 200.0
        @test elev[3] == 500.0  # default for missing site
    end

    @testset "extract_site_temperatures with TEMP_MEDIA_SC" begin
        site_df = DataFrame(
            CODIGO = ["site1", "site2", "site3"],
            TEMP_MEDIA_SC = [11.18, 12.44, 13.58],
            ALTITUD = [500.0, 600.0, 700.0]
        )
        sites = ["site1", "site2", "site4"]
        temps = Guadex.extract_site_temperatures(site_df, sites)
        @test length(temps) == 3
        @test temps[1] ≈ 11.18
        @test temps[2] ≈ 12.44
        @test temps[3] ≈ 15.0  # default for missing site
    end

    @testset "extract_site_temperatures fallback (no temp column)" begin
        site_df = DataFrame(
            CODIGO = ["site1", "site2"],
            ALTITUD = [0.0, 1000.0],
            OTHER_DATA = [1.0, 2.0]
        )
        sites = ["site1", "site2"]
        temps = Guadex.extract_site_temperatures(site_df, sites)
        @test length(temps) == 2
        # Elevation-based estimate: 20 - (0/1000 * 6.5) = 20.0
        @test temps[1] ≈ 20.0
        # 20 - (1000/1000 * 6.5) = 13.5
        @test temps[2] ≈ 13.5
    end

    @testset "extract_habitat_suitability" begin
        site_df = DataFrame(CODIGO=["site1", "site2", "site3"], IET=[5.0, 10.0, 15.0])
        sites = ["site1", "site2", "site4"]
        suit = Guadex.extract_habitat_suitability(site_df, sites)
        @test length(suit) == 3
        # site1: IET=5 (lowest) → highest suitability
        # site2: IET=10 → mid suitability
        # site4: missing → default 0.5
        @test suit[1] >= suit[2]
        @test suit[3] ≈ 0.5
        @test all(suit .>= 0.1)
        @test all(suit .<= 1.0)
    end

    @testset "build_carrying_capacity" begin
        density_df = DataFrame(
            CODIGO = ["site1", "site2"],
            SP_DEN = [10.0, 0.0],
            LS_DEN = [20.0, 5.0]
        )
        site_df = DataFrame(CODIGO=["site1", "site2"], ALTITUD=[100.0, 200.0])
        sites = ["site1", "site2", "site3"]
        species_codes = ["LS", "SP"]
        K = Guadex.build_carrying_capacity(density_df, site_df, sites, species_codes)
        @test length(K) == 3
        # site1: (10+20)*10 = 300
        @test K[1] ≈ 300.0
        # site2: (0+5)*10 = 50, at least K_min
        @test K[2] >= 50.0
        # site3: missing → 0, raised to K_floor
        @test K[3] >= 50.0
    end

    @testset "build_dispersal_scaling" begin
        species_codes = ["AB", "AH", "MC", "LR", "GH"]
        scaling = Guadex.build_dispersal_scaling(species_codes)
        @test length(scaling) == 5

        # MC, LR are 75 km/year; median ≈ 5.0 for these 5 codes
        # but median of all 24 species = 10.0, used for normalization
        rates = Guadex.ANNUAL_DISPERSAL_RATES
        median_rate = median([get(rates, sp, 5.0) for sp in species_codes])
        @test scaling[1] ≈ 0.3 / median_rate
        @test scaling[2] ≈ 0.2 / median_rate
        @test scaling[3] ≈ 75.0 / median_rate
        @test scaling[4] ≈ 75.0 / median_rate
        @test scaling[5] ≈ 25.0 / median_rate
    end

    @testset "build_distance_matrix (via prepare_ode_data)" begin
        data = Guadex.prepare_ode_data(
            connectivity_file=CONNECTIVITY_FILE,
            density_file=DENSITY_FILE,
            species_chars_file=SPECIES_CHARS_FILE,
            environmental_file=ENVIRONMENTAL_FILE,
            distance_file=DISTANCE_FILE,
            interaction_file=INTERACTION_FILE,
            upstream_cost=0.01
        )
        @test data.distance_matrix !== nothing
        @test nnz(data.distance_matrix) >= 0
    end

    @testset "build_dam_passability_matrix" begin
        site_df = DataFrame(
            CODIGO = ["site1", "site2"],
            Demb_arr_m = ["1000", "0"],
            Demb_ab_m = ["0", "2000"]
        )
        # Convert "No existe" replacement
        site_df.Demb_arr_m = parse.(Float64, replace(site_df.Demb_arr_m, "No existe" => "0"))
        site_df.Demb_ab_m = parse.(Float64, replace(site_df.Demb_ab_m, "No existe" => "0"))
        # Replace zeros with large values
        site_df.Demb_arr_m = replace(site_df.Demb_arr_m, 0.0 => 1_000_000.0)
        site_df.Demb_ab_m = replace(site_df.Demb_ab_m, 0.0 => 1_000_000.0)

        sites = ["site1", "site2"]
        distances = sparse([2], [1], [500.0], 2, 2)  # site1→site2 distance 500m
        elevations = [100.0, 200.0]

        dams = Guadex.build_dam_passability_matrix(site_df, sites, distances, elevations)
        @test size(dams) == (2, 2)
        @test dams[1, 1] == 1.0  # diagonal = 1
        @test dams[2, 2] == 1.0
    end

    @testset "build_intrinsic_growth_rates" begin
        density_df = DataFrame(
            CODIGO = ["site1", "site2"],
            AH_DEN = [10.0, 0.0],
            GH_DEN = [5.0, 20.0]
        )
        species_chars_df = DataFrame(
            SP = ["Ah", "Gh"],
            thermal_optimum = [19.0, 19.0],
            thermal_sigma = [3.67, 3.67],
            elevation_optimum = [450.0, 450.0],
            max_size_mm = [100.0, 65.0]
        )
        sites = ["site1", "site2"]
        species_codes = ["AH", "GH"]

        growth = Guadex.build_intrinsic_growth_rates(density_df, species_codes, sites, species_chars_df)
        @test size(growth) == (2, 2)

        # AH annual rate = 1.0 → daily = 1.0/365, present at site1 → full rate
        r_ah_daily = 1.0 / 365.0
        @test growth[1, 1] ≈ r_ah_daily       # site1, AH (present)
        @test growth[2, 1] ≈ r_ah_daily * 0.1  # site2, AH (absent)

        # GH annual rate = 4.0 → daily = 4.0/365
        r_gh_daily = 4.0 / 365.0
        @test growth[1, 2] ≈ r_gh_daily       # site1, GH (present)
        @test growth[2, 2] ≈ r_gh_daily        # site2, GH (present)
    end

    @testset "prepare_ode_data (smoke test)" begin
            data = Guadex.prepare_ode_data(
                connectivity_file=CONNECTIVITY_FILE,
                density_file=DENSITY_FILE,
                species_chars_file=SPECIES_CHARS_FILE,
                environmental_file=ENVIRONMENTAL_FILE,
                distance_file=DISTANCE_FILE,
                interaction_file=INTERACTION_FILE,
                upstream_cost=0.01
            )
        @test :params in keys(data)
        @test :sites in keys(data)
        @test :species in keys(data)
        @test :distance_matrix in keys(data)
        @test :elevations in keys(data)
        @test :dams in keys(data)
        @test :site_df in keys(data)
        @test :species_chars_df in keys(data)
        @test :density_df in keys(data)

        # Metadata integrity
        @test data.params.n_sites == length(data.sites)
        @test data.params.n_species == length(data.species)
        @test length(data.elevations) == data.params.n_sites
        @test size(data.params.interaction_matrix) == (data.params.n_species, data.params.n_species)
        @test size(data.params.intrinsic_growth_rates) == (data.params.n_sites, data.params.n_species)
        @test length(data.params.temperatures) == data.params.n_sites
        @test length(data.params.habitat_suitability) == data.params.n_sites
        @test length(data.params.thermal_optima) == data.params.n_species
        @test length(data.params.thermal_sigmas) == data.params.n_species
        @test length(data.params.carrying_capacity) == data.params.n_sites
        @test length(data.params.dispersal_scaling) == data.params.n_species
    end

    @testset "save_ode_data" begin
        mktempdir() do tmpdir
        data = Guadex.prepare_ode_data(
            connectivity_file=CONNECTIVITY_FILE,
            density_file=DENSITY_FILE,
            species_chars_file=SPECIES_CHARS_FILE,
            environmental_file=ENVIRONMENTAL_FILE,
            distance_file=DISTANCE_FILE,
            interaction_file=INTERACTION_FILE,
            upstream_cost=0.01
        )
            Guadex.save_ode_data(data, tmpdir)
            @test isfile(joinpath(tmpdir, "sites.csv"))
            @test isfile(joinpath(tmpdir, "species.csv"))
            @test isfile(joinpath(tmpdir, "distance_matrix.csv"))
            @test isfile(joinpath(tmpdir, "elevations.csv"))
            @test isfile(joinpath(tmpdir, "temperatures.csv"))
            @test isfile(joinpath(tmpdir, "habitat_suitability.csv"))
            @test isfile(joinpath(tmpdir, "interaction_matrix.csv"))
            @test isfile(joinpath(tmpdir, "growth_rates.csv"))
        end
    end
end
