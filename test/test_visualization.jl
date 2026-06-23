# =============================================================================
# Smoke tests for visualization functions
# =============================================================================

@testset "Visualization smoke tests" begin
    # Create minimal data for visualization tests
    data = Guadex.prepare_ode_data(
        connectivity_file=CONNECTIVITY_FILE,
        density_file=DENSITY_FILE,
        species_chars_file=SPECIES_CHARS_FILE,
        environmental_file=ENVIRONMENTAL_FILE,
        distance_file=DISTANCE_FILE,
        interaction_file=INTERACTION_FILE,
        upstream_cost=0.01
    )
    ndays = 5

    # Run a tiny simulation for plotting
    density_cols = [Symbol("$(sp)_DEN") for sp in data.species]
    density_df_filtered = filter(row -> row.CODIGO in data.sites, data.density_df)
    u0 = Matrix(density_df_filtered[:, density_cols])
    u0 = max.(u0, 0.0)
    u0_flat = vec(u0)

    prob = ODEProblem(
        Guadex.metacommunity_ode!, u0_flat,
        (0.0, Float64(ndays)), data.params
    )
    sol = solve(prob, Tsit5(), saveat=0:1.0:Float64(ndays),
                reltol=1e-4, abstol=1e-4)

    @testset "plot_avg_total_biomass" begin
        fig = Guadex.plot_avg_total_biomass(sol, data.sites, data.species)
        @test fig !== nothing
    end

    @testset "plot_avg_species_richness" begin
        fig = Guadex.plot_avg_species_richness(sol, data.sites, data.species)
        @test fig !== nothing
    end

    @testset "plot_total_biomass" begin
        fig = Guadex.plot_total_biomass(sol, data.sites, data.species)
        @test fig !== nothing
    end

    @testset "plot_species_richness" begin
        fig = Guadex.plot_species_richness(sol, data.sites, data.species)
        @test fig !== nothing
    end

    @testset "plot_combined_analysis" begin
        fig = Guadex.plot_combined_analysis(
            sol, data.site_df, data.sites, data.species,
            data.distance_matrix
        )
        @test fig !== nothing
    end

    @testset "plot_sites_map" begin
        fig = Guadex.plot_sites_map(data.site_df; color_by=:ALTITUD)
        @test fig !== nothing
    end

    @testset "save_figure" begin
        mktempdir() do tmpdir
            fig = Guadex.plot_avg_total_biomass(sol, data.sites, data.species)
            Guadex.save_figure(fig, joinpath(tmpdir, "test.png"))
            @test isfile(joinpath(tmpdir, "test.png"))
        end
    end
end
