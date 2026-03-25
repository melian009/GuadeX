"""
    Visualization Module for ODE Model Results

This module provides functions to visualize:
1. ODE simulation results (population dynamics over time)
2. Site locations on a map with coordinates
3. Connectivity network between sites

All plots use Makie.jl for high-quality rendering.
"""

"""
    plot_ode_solution(sol, sites, species; kwargs)

Plot the population dynamics from ODE solution.

# Arguments
- `sol`: DifferentialEquations solution object
- `sites`: Vector of site codes
- `species`: Vector of species codes
- `max_species_to_plot`: Maximum number of species to plot (default: 10)
- `max_sites_to_plot`: Maximum number of sites to plot (default: 6)
- `figure_size`: Tuple for figure dimensions (default: (1200, 800))

# Returns
- A Figure object with subplots showing population dynamics
"""
function plot_ode_solution(sol, sites, species;
    max_species_to_plot::Int = 10,
    max_sites_to_plot::Int = 6,
    figure_size::Tuple = (1200, 800))

    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    # Reshape solution to (n_sites, n_species, n_timepoints)
    u_reshaped = [reshape(sol.u[i], n_sites, n_species) for i in 1:n_timepoints]

    # Determine layout
    n_sp = min(max_species_to_plot, n_species)
    n_st = min(max_sites_to_plot, n_sites)
    n_rows = ceil(Int, n_sp / 3)
    n_cols = min(3, n_sp)

    fig = Figure(figure_size = figure_size)

    time = sol.t

    # Get species to plot (first n_sp)
    species_to_plot = species[1:n_sp]

    for (sp_idx, sp_code) in enumerate(species_to_plot)
        ax = Axis(fig[floor(Int, (sp_idx-1)/n_cols) + 1, (sp_idx-1) % n_cols + 1],
            title = "$sp_code - Population Dynamics",
            xlabel = "Time",
            ylabel = "Population Density")

        # Plot each site for this species
        for site_idx in 1:n_st
            pop_values = [u_reshaped[t][site_idx, sp_idx] for t in 1:n_timepoints]
            lines!(ax, time, pop_values, label = sites[site_idx], linewidth = 1.5)
        end

        if n_st > 0
            axislegend(ax, position = :rt)
        end
    end

    # Add summary title
    Label(fig[0, :], "ODE Simulation Results: Population Dynamics Over Time", 
        fontsize = 16)

    return fig
end

"""
    plot_total_biomass(sol, sites, species; kwargs)

Plot total biomass across all species over time.

# Arguments
- `sol`: DifferentialEquations solution object
- `sites`: Vector of site codes
- `species`: Vector of species codes

# Returns
- A Figure object showing total biomass trajectories
"""
function plot_total_biomass(sol, sites, species;
    figure_size::Tuple = (900, 600))

    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    # Calculate total biomass per site at each timepoint
    total_biomass = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        total_biomass[t, :] .= sum(u_matrix, dims = 2)[:]
    end

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Total Biomass Over Time by Site",
        xlabel = "Time",
        ylabel = "Total Biomass")

    # Color gradient for sites
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :gray]

    for site_idx in 1:n_sites
        lines!(ax, time, total_biomass[:, site_idx],
            label = sites[site_idx],
            color = colors[(site_idx-1) % length(colors) + 1],
            linewidth = 2)
    end

    axislegend(ax, position = :rt)

    return fig
end

"""
    plot_species_richness(sol, sites, species; kwargs)

Plot species richness (number of species with population > threshold) over time.

# Arguments
- `sol`: DifferentialEquations solution object
- `sites`: Vector of site codes
- `species`: Vector of species codes
- `threshold`: Population threshold for species presence (default: 0.1)

# Returns
- A Figure object showing species richness trajectories
"""
function plot_species_richness(sol, sites, species;
    threshold::Float64 = 0.1,
    figure_size::Tuple = (900, 600))

    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    # Calculate richness per site at each timepoint
    richness = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        for site_idx in 1:n_sites
            richness[t, site_idx] = sum(u_matrix[site_idx, :] .> threshold)
        end
    end

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Species Richness Over Time by Site (threshold = $threshold)",
        xlabel = "Time",
        ylabel = "Number of Species")

    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :brown, :pink, :gray]

    for site_idx in 1:n_sites
        lines!(ax, time, richness[:, site_idx],
            label = sites[site_idx],
            color = colors[(site_idx-1) % length(colors) + 1],
            linewidth = 2)
    end

    axislegend(ax, position = :rt)

    return fig
end

"""
    plot_sites_map(site_df; kwargs)

Plot sites on a map using their UTM coordinates.

# Arguments
- `site_df`: DataFrame with site data including UTMX, UTMY, CODIGO columns
- `color_by`: Column to color points by (default: :ALTITUD)
- `figure_size`: Tuple for figure dimensions (default: (1000, 800))
- `markersize`: Size of site markers (default: 15)

# Returns
- A Figure object showing site locations
"""
function plot_sites_map(site_df;
    color_by::Symbol = :ALTITUD,
    figure_size::Tuple = (1000, 800),
    markersize::Int = 15)

    # Extract coordinates
    x = Float64.(site_df.UTMX)
    y = Float64.(site_df.UTMY)

    # Get coloring values
    color_values = Float64.(site_df[!, color_by])

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Site Locations ($(color_by))",
        xlabel = "UTM X (meters)",
        ylabel = "UTM Y (meters)")

    # Create scatter plot with colorbar
    sc = scatter!(ax, x, y,
        color = color_values,
        markersize = markersize,
        colormap = :viridis,
        strokecolor = :black,
        strokewidth = 1)

    # Add colorbar
    Colorbar(fig[1, 2], sc, label = string(color_by))

    # Add site labels (only for subset if many sites)
    n_sites = length(x)
    label_fraction = max(1, floor(Int, n_sites / 20))
    for i in 1:n_sites
        if i % label_fraction == 0 || i <= 10
            text!(ax, x[i], y[i], text = string(site_df.CODIGO[i]),
                fontsize = 7, align = (:center, :bottom), offset = (0, 5))
        end
    end

    return fig
end

"""
    plot_site_connectivity_map(site_df, distance_matrix, sites; kwargs)

Plot sites on a map showing connectivity network.

# Arguments
- `site_df`: DataFrame with site data including UTMX, UTMY, CODIGO columns
- `distance_matrix`: Sparse distance matrix between sites
- `sites`: Vector of site codes (matching site_df order)
- `max_connections_to_show`: Maximum number of connections to display (default: 100)
- `figure_size`: Tuple for figure dimensions (default: (1200, 900))
- `connection_alpha`: Transparency of connection lines (default: 0.3)

# Returns
- A Figure object showing site locations with connectivity network
"""
function plot_site_connectivity_map(site_df, distance_matrix, sites;
    max_connections_to_show::Int = 100,
    figure_size::Tuple = (1200, 900),
    connection_alpha::Float64 = 0.3)

    # Create site to index mapping
    site_to_idx = Dict(s => i for (i, s) in enumerate(sites))

    # Filter site_df to only include sites in our list
    filtered_df = filter(row -> row.CODIGO in sites, site_df)

    # Extract coordinates for filtered sites
    x = Float64.(filtered_df.UTMX)
    y = Float64.(filtered_df.UTMY)
    site_codes = String.(filtered_df.CODIGO)
    elevations = Float64.(filtered_df.ALTITUD)

    n_sites = length(site_codes)

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Site Network Connectivity",
        xlabel = "UTM X (meters)",
        ylabel = "UTM Y (meters)")

    # Plot elevation-colored scatter
    sc = scatter!(ax, x, y,
        color = elevations,
        markersize = 20,
        colormap = :terrain,
        strokecolor = :black,
        strokewidth = 2)

    Colorbar(fig[1, 2], sc, label = "Elevation (m)")

    # Get non-zero connections from distance matrix
    I, J, V = findnz(distance_matrix)

    # Filter and limit connections
    connection_count = 0
    for idx in 1:length(I)
        if connection_count >= max_connections_to_show
            break
        end

        i, j = I[idx], J[idx]
        if i > n_sites || j > n_sites
            continue
        end

        # Get coordinates
        x1, y1 = x[j], y[j]  # origin
        x2, y2 = x[i], y[i]  # destination

        # Connection strength based on distance (inverse)
        strength = min(1.0, 1000.0 / (V[idx] + 1))

        lines!(ax, [x1, x2], [y1, y2],
            color = (:blue, connection_alpha * strength),
            linewidth = 0.5 + strength * 2)

        connection_count += 1
    end

    # Add site labels for subset
    label_fraction = max(1, floor(Int, n_sites / 25))
    for i in 1:n_sites
        if i % label_fraction == 0 || i <= 10
            text!(ax, x[i], y[i], text = site_codes[i],
                fontsize = 7, align = (:center, :bottom), offset = (0, 5))
        end
    end

    # Add legend text
    text!(ax, 0.02, 0.98, text = "Lines show dispersal connections\n(Thickness = connection strength)",
        position = (minimum(x) + (maximum(x) - minimum(x)) * 0.02, maximum(y) * 0.98),
        fontsize = 10, align = (:left, :top), space = :relative)

    return fig
end

"""
    plot_subcatchment_network(site_df, sites; kwargs)

Plot sites colored by subcatchment showing network structure.

# Arguments
- `site_df`: DataFrame with site data including UTMX, UTMY, CODIGO, CODIGO_S columns
- `sites`: Vector of site codes (matching site_df order)
- `figure_size`: Tuple for figure dimensions (default: (1200, 900))

# Returns
- A Figure object showing subcatchment network structure
"""
function plot_subcatchment_network(site_df, sites;
    figure_size::Tuple = (1200, 900))

    # Filter to sites in our list
    filtered_df = filter(row -> row.CODIGO in sites, site_df)

    # Get unique subcatchments
    subcatchments = unique(filtered_df.CODIGO_S)
    n_subcatchments = length(subcatchments)

    # Create subcatchment colors
    colors = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :brown,
              :pink, :gray, :yellow, :teal, :navy, :olive, :maroon, :coral]

    # Map subcatchments to colors
    sc_to_color = Dict(sc => colors[(i-1) % length(colors) + 1]
                       for (i, sc) in enumerate(subcatchments))

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Site Network by Subcatchment",
        xlabel = "UTM X (meters)",
        ylabel = "UTM Y (meters)")

    # Plot each subcatchment
    for sc in subcatchments
        sc_df = filter(row -> row.CODIGO_S == sc, filtered_df)
        x = Float64.(sc_df.UTMX)
        y = Float64.(sc_df.UTMY)

        scatter!(ax, x, y,
            color = sc_to_color[sc],
            markersize = 15,
            label = "Subcatchment $sc",
            strokecolor = :black,
            strokewidth = 1)

        # Connect sites within subcatchment (sorted by distance to river)
        if nrow(sc_df) > 1
            # Sort by distance to Guadalquivir
            sorted_df = sort(sc_df, :("Dist.Guadalq.(m)"))
            sx = Float64.(sorted_df.UTMX)
            sy = Float64.(sorted_df.UTMY)
            lines!(ax, sx, sy, color = sc_to_color[sc], linewidth = 1.5, alpha = 0.5)
        end
    end

    axislegend(ax, position = :rt, fontsize = 8)

    return fig
end

"""
    plot_combined_analysis(sol, site_df, sites, species, distance_matrix; kwargs)

Create a comprehensive combined visualization.

# Arguments
- `sol`: DifferentialEquations solution object
- `site_df`: DataFrame with site data
- `sites`: Vector of site codes
- `species`: Vector of species codes
- `distance_matrix`: Sparse distance matrix
- `figure_size`: Tuple for figure dimensions (default: (1400, 1000))

# Returns
- A Figure object with multiple subplots
"""
function plot_combined_analysis(sol, site_df, sites, species, distance_matrix;
    figure_size::Tuple = (1400, 1000))

    n_sites = length(sites)
    n_species = length(species)

    fig = Figure(figure_size = figure_size)

    # Layout: 2x2 grid
    # Top-left: Total biomass over time
    # Top-right: Species richness over time
    # Bottom: Site map with connectivity

    # Top-left: Total Biomass
    ax1 = Axis(fig[1, 1],
        title = "Total Biomass Over Time",
        xlabel = "Time",
        ylabel = "Total Biomass")

    time = sol.t
    total_biomass = [sum(reshape(sol.u[t], n_sites, n_species)) for t in 1:length(sol.t)]
    lines!(ax1, time, total_biomass, linewidth = 2, color = :blue)

    # Top-right: Species Richness
    ax2 = Axis(fig[1, 2],
        title = "Total Species Richness Over Time",
        xlabel = "Time",
        ylabel = "Richness")

    richness = [sum(reshape(sol.u[t], n_sites, n_species) .> 0.1) for t in 1:length(sol.t)]
    lines!(ax2, time, richness, linewidth = 2, color = :green)

    # Bottom: Site Map with Connectivity
    ax3 = Axis(fig[2, 1:2],
        title = "Site Network",
        xlabel = "UTM X (meters)",
        ylabel = "UTM Y (meters)")

    # Filter site_df
    filtered_df = filter(row -> row.CODIGO in sites, site_df)
    x = Float64.(filtered_df.UTMX)
    y = Float64.(filtered_df.UTMY)
    elevations = Float64.(filtered_df.ALTITUD)

    # Plot sites
    sc = scatter!(ax3, x, y,
        color = elevations,
        markersize = 15,
        colormap = :terrain,
        strokecolor = :black,
        strokewidth = 1)

    Colorbar(fig[2, 3], sc, label = "Elevation (m)")

    # Plot connections (subset)
    I, J, V = findnz(distance_matrix)
    max_shown = min(100, length(I))
    for idx in 1:max_shown
        i, j = I[idx], J[idx]
        if i <= length(x) && j <= length(x)
            lines!(ax3, [x[j], x[i]], [y[j], y[i]],
                color = (:blue, 0.2),
                linewidth = 0.5)
        end
    end

    return fig
end

"""
    save_figure(fig, filename::String; kwargs)

Save a Makie figure to file.

# Arguments
- `fig`: Makie Figure object
- `filename`: Output filename (with extension)
- `resolution`: Tuple for image resolution (default: (1200, 900))
"""
function save_figure(fig, filename::String; resolution::Tuple = (1200, 900))
    Makie.save(filename, fig, resolution = resolution)
    println("Figure saved to: $filename")
end
