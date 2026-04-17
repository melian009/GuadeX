"""
    Visualization Module for ODE Model Results

This module provides functions to visualize:
1. ODE simulation results (population dynamics over time)
2. Site locations on a map with coordinates
3. Connectivity network between sites

All plots use Makie.jl with GeoMakie.jl for geographic mapping.
"""


"""
    utm_to_latlon(easting, northing, zone::Int=30)

Convert UTM coordinates to WGS84 latitude/longitude.
Approximate conversion for UTM Zone 30N (covers western Europe, Spain).

# Arguments
- `easting`: UTM easting (meters)
- `northing`: UTM northing (meters)
- `zone`: UTM zone number (default: 30)

# Returns
- Tuple of (latitude, longitude) in decimal degrees
"""
function utm_to_latlon(easting, northing, zone::Int=30)
    # WGS84 parameters
    a = 6378137.0  # semi-major axis
    f = 1/298.257223563  # flattening
    k0 = 0.9996  # scale factor

    # Derived parameters
    e = sqrt(2*f - f^2)
    e2 = e^2 / (1 - e^2)

    # False easting and northing
    FE = 500000.0
    FN = (zone >= 33) ? 0.0 : 0.0  # Southern hemisphere offset

    # Remove false easting/northing
    x = easting - FE
    y = northing - FN

    # Footprint latitude
    M = y / k0
    mu = M / (a * (1 - e^2/4 - 3*e^4/64 - 5*e^6/256))

    e1 = (1 - sqrt(1 - e^2)) / (1 + sqrt(1 - e^2))

    phi1 = mu +
           (3*e1/2 - 27*e1^3/32) * sin(2*mu) +
           (21*e1^2/16 - 55*e1^4/32) * sin(4*mu) +
           (151*e1^3/96) * sin(6*mu) +
           (1097*e1^4/512) * sin(8*mu)

    # Latitude
    sin_phi1 = sin(phi1)
    cos_phi1 = cos(phi1)
    tan_phi1 = tan(phi1)

    N1 = a / sqrt(1 - e^2 * sin_phi1^2)
    T1 = tan_phi1^2
    C1 = e2 * cos_phi1^2
    R1 = a * (1 - e^2) / ((1 - e^2 * sin_phi1^2)^1.5)
    D = x / (N1 * k0)

    lat = phi1 - (N1 * tan_phi1 / R1) * (
        D^2/2 -
        (5 + 3*T1 + 10*C1 - 4*C1^2 - 9*e2) * D^4/24 +
        (61 + 90*T1 + 298*C1 + 45*T1^2 - 252*e2 - 3*C1^2) * D^6/720
    )

    # Longitude
    lon = (D -
        (1 + 2*T1 + C1) * D^3/6 +
        (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*e2 + 24*T1^2) * D^5/120
    ) / cos_phi1

    # Adjust for zone
    lon = lon + (zone - 1) * 6 - 180 + 3  # Central meridian

    # Convert to degrees
    lat_deg = lat * 180 / π
    lon_deg = lon * 180 / π

    return (lat_deg, lon_deg)
end

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
            axislegend(ax, position = (1, 0.5), fontsize = 8)
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

    axislegend(ax, position = (1, 0.5), fontsize = 8)

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

    axislegend(ax, position = (1, 0.5), fontsize = 8)

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
function plot_sites_map(site_df; color_by::Symbol = :ALTITUD, figure_size::Tuple = (1000, 800), markersize::Int = 15)
    utm_x = Float64.(site_df.UTMX)
    utm_y = Float64.(site_df.UTMY)

    lats = Float64[]
    lons = Float64[]
    for (ex, ey) in zip(utm_x, utm_y)
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(lats, lat)
        push!(lons, lon)
    end

    color_values = Float64.(site_df[!, color_by])
    n_sites = length(lats)

    fig = Figure(figure_size = figure_size)

    ax = GeoAxis(fig[1, 1];
        title = "Site Locations in Guadalquivir Basin ($(color_by))",
        xlabel = "Longitude",
        ylabel = "Latitude",
        dest = "+proj=latlong",
        lonlims = extrema(lons) .+ (-0.5, 0.5),
        latlims = extrema(lats) .+ (-0.5, 0.5))

    land = GeoMakie.land()
    poly!(ax, land; color = (:lightgray, 0.3), strokecolor = :gray, strokewidth = 0.5)
    lines!(ax, GeoMakie.coastlines(); color = :darkgray, linewidth = 1)

    sc = scatter!(ax, lons, lats;
        color = color_values,
        markersize = markersize,
        colormap = :viridis,
        strokecolor = :black,
        strokewidth = 1)

    Colorbar(fig[1, 2], sc, label = string(color_by))

    label_fraction = max(1, floor(Int, n_sites / 20))
    for i in 1:n_sites
        if i % label_fraction == 0 || i <= 10
            text!(ax, lons[i], lats[i], text = string(site_df.CODIGO[i]),
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
    connection_alpha::Float64 = 0.6)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)

    # Create mapping from site code to position in filtered_df (1-indexed)
    site_to_pos = Dict{String, Int}()
    for (pos, row) in enumerate(eachrow(filtered_df))
        site_to_pos[row.CODIGO] = pos
    end

    # Extract UTM coordinates
    utm_x = Float64.(filtered_df.UTMX)
    utm_y = Float64.(filtered_df.UTMY)

    lats = Float64[]
    lons = Float64[]
    for (ex, ey) in zip(utm_x, utm_y)
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(lats, lat)
        push!(lons, lon)
    end

    elevations = Float64.(filtered_df.ALTITUD)
    site_codes = String.(filtered_df.CODIGO)
    n_sites = length(site_codes)

    fig = Figure(figure_size = figure_size)

    ax = GeoAxis(fig[1, 1];
        title = "Site Network Connectivity - Guadalquivir Basin",
        xlabel = "Longitude",
        ylabel = "Latitude",
        dest = "+proj=latlong",
        lonlims = extrema(lons) .+ (-0.5, 0.5),
        latlims = extrema(lats) .+ (-0.5, 0.5))

    land = GeoMakie.land()
    poly!(ax, land; color = (:lightgray, 0.3), strokecolor = :gray, strokewidth = 0.5)
    lines!(ax, GeoMakie.coastlines(); color = :darkgray, linewidth = 1)

    sc = scatter!(ax, lons, lats;
        color = elevations,
        markersize = 25,
        colormap = :terrain,
        strokecolor = :black,
        strokewidth = 2)

    Colorbar(fig[1, 2], sc, label = "Elevation (m)")

    # Get non-zero connections from distance matrix
    I, J, V = findnz(distance_matrix)

    # Filter connections to only those between sites in our filtered list
    valid_connections = []
    for idx in 1:length(I)
        i, j = I[idx], J[idx]

        # Get site codes for these indices (in the original sites order)
        if i > length(sites) || j > length(sites)
            continue
        end
        site_i = sites[i]
        site_j = sites[j]

        # Check if both sites are in our filtered set
        if haskey(site_to_pos, site_i) && haskey(site_to_pos, site_j)
            pos_i = site_to_pos[site_i]
            pos_j = site_to_pos[site_j]
            push!(valid_connections, (pos_i, pos_j, V[idx]))
        end
    end

    # Sort by distance (shorter connections first, more important)
    sort!(valid_connections, by = x -> x[3])

    # Limit and plot connections
    connection_count = 0
    for (pos_i, pos_j, dist) in valid_connections
        if connection_count >= max_connections_to_show
            break
        end
        x1, y1 = lons[pos_j], lats[pos_j]
        x2, y2 = lons[pos_i], lats[pos_i]
        strength = min(1.0, 5000.0 / (dist + 1))
        lines!(ax, [x1, x2], [y1, y2],
            color = (:blue, connection_alpha * strength),
            linewidth = 1.0 + strength * 3)
        connection_count += 1
    end

    # Add site labels for subset
    label_fraction = max(1, floor(Int, n_sites / 25))
    for i in 1:n_sites
        if i % label_fraction == 0 || i <= 10
            text!(ax, lons[i], lats[i], text = site_codes[i],
                fontsize = 7, align = (:center, :bottom), offset = (0, 5))
        end
    end

    return fig
end

"""
    plot_subcatchment_network(site_df, sites; kwargs)
"""
function plot_subcatchment_network(site_df, sites; figure_size::Tuple = (1200, 900))
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

    all_utm_x = Float64.(filtered_df.UTMX)
    all_utm_y = Float64.(filtered_df.UTMY)
    all_lats = Float64[]
    all_lons = Float64[]
    for (ex, ey) in zip(all_utm_x, all_utm_y)
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(all_lats, lat)
        push!(all_lons, lon)
    end

    fig = Figure(figure_size = figure_size)

    ax = GeoAxis(fig[1, 1];
        title = "Site Network by Subcatchment - Guadalquivir Basin",
        xlabel = "Longitude",
        ylabel = "Latitude",
        dest = "+proj=latlong",
        lonlims = extrema(all_lons) .+ (-0.5, 0.5),
        latlims = extrema(all_lats) .+ (-0.5, 0.5))

    land = GeoMakie.land()
    poly!(ax, land; color = (:lightgray, 0.3), strokecolor = :gray, strokewidth = 0.5)
    lines!(ax, GeoMakie.coastlines(); color = :darkgray, linewidth = 1)

    for sc in subcatchments
        sc_df = filter(row -> row.CODIGO_S == sc, filtered_df)

        lats = Float64[]
        lons = Float64[]
        for (ex, ey) in zip(Float64.(sc_df.UTMX), Float64.(sc_df.UTMY))
            lat, lon = utm_to_latlon(ex, ey, 30)
            push!(lats, lat)
            push!(lons, lon)
        end

        scatter!(ax, lons, lats;
            color = sc_to_color[sc],
            markersize = 15,
            label = "Subcatchment $sc",
            strokecolor = :black,
            strokewidth = 1)

        if nrow(sc_df) > 1
            sorted_df = sort(sc_df, :("Dist.Guadalq.(m)"))
            slats = Float64[]
            slons = Float64[]
            for (ex, ey) in zip(Float64.(sorted_df.UTMX), Float64.(sorted_df.UTMY))
                lat, lon = utm_to_latlon(ex, ey, 30)
                push!(slats, lat)
                push!(slons, lon)
            end
            lines!(ax, slons, slats, color = sc_to_color[sc], linewidth = 1.5, alpha = 0.5)
        end
    end

    axislegend(ax, position = (1, 0.5), fontsize = 8)
    return fig
end

"""
    plot_combined_analysis(sol, site_df, sites, species, distance_matrix; kwargs)
"""
function plot_combined_analysis(sol, site_df, sites, species, distance_matrix; figure_size::Tuple = (1400, 1000))
    n_sites = length(sites)
    n_species = length(species)

    fig = Figure(figure_size = figure_size)

    ax1 = Axis(fig[1, 1], title = "Total Biomass Over Time", xlabel = "Time", ylabel = "Total Biomass")
    time = sol.t
    total_biomass = [sum(reshape(sol.u[t], n_sites, n_species)) for t in 1:length(sol.t)]
    lines!(ax1, time, total_biomass, linewidth = 2, color = :blue)

    ax2 = Axis(fig[1, 2], title = "Total Species Richness Over Time", xlabel = "Time", ylabel = "Richness")
    richness = [sum(reshape(sol.u[t], n_sites, n_species) .> 0.1) for t in 1:length(sol.t)]
    lines!(ax2, time, richness, linewidth = 2, color = :green)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)

    lats = Float64[]
    lons = Float64[]
    for (ex, ey) in zip(Float64.(filtered_df.UTMX), Float64.(filtered_df.UTMY))
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(lats, lat)
        push!(lons, lon)
    end

    elevations = Float64.(filtered_df.ALTITUD)

    ax3 = GeoAxis(fig[2, 1:2];
        title = "Site Network - Guadalquivir Basin",
        xlabel = "Longitude",
        ylabel = "Latitude",
        dest = "+proj=latlong",
        lonlims = extrema(lons) .+ (-0.5, 0.5),
        latlims = extrema(lats) .+ (-0.5, 0.5))

    land = GeoMakie.land()
    poly!(ax3, land; color = (:lightgray, 0.3), strokecolor = :gray, strokewidth = 0.5)
    lines!(ax3, GeoMakie.coastlines(); color = :darkgray, linewidth = 1)

    sc = scatter!(ax3, lons, lats;
        color = elevations,
        markersize = 15,
        colormap = :terrain,
        strokecolor = :black,
        strokewidth = 1)

    Colorbar(fig[2, 3], sc, label = "Elevation (m)")

    # Create site position mapping for connectivity
    site_to_pos = Dict{String, Int}()
    for (pos, row) in enumerate(eachrow(filtered_df))
        site_to_pos[row.CODIGO] = pos
    end

    # Plot connections (subset)
    I, J, V = findnz(distance_matrix)
    valid_connections = []
    for idx in 1:length(I)
        i, j = I[idx], J[idx]
        if i > length(sites) || j > length(sites)
            continue
        end
        site_i = sites[i]
        site_j = sites[j]
        if haskey(site_to_pos, site_i) && haskey(site_to_pos, site_j)
            pos_i = site_to_pos[site_i]
            pos_j = site_to_pos[site_j]
            push!(valid_connections, (pos_i, pos_j, V[idx]))
        end
    end

    # Sort and limit connections
    sort!(valid_connections, by = c -> c[3])
    max_shown = min(100, length(valid_connections))

    for idx in 1:max_shown
        pos_i, pos_j, dist = valid_connections[idx]
        strength = min(1.0, 5000.0 / (dist + 1))
        lines!(ax3, [lons[pos_j], lons[pos_i]], [lats[pos_j], lats[pos_i]],
            color = (:blue, 0.6 * strength),
            linewidth = 1.0 + strength * 2)
    end

    return fig
end

"""
    plot_avg_total_biomass(sol, sites, species; kwargs)

Plot average total biomass across all sites over time with ±1 std band.

# Arguments
- `sol`: DifferentialEquations solution object
- `sites`: Vector of site codes
- `species`: Vector of species codes

# Returns
- A Figure object showing average total biomass trajectory with uncertainty band
"""
function plot_avg_total_biomass(sol, sites, species;
    figure_size::Tuple = (900, 600))

    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    total_biomass = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        total_biomass[t, :] .= sum(u_matrix, dims = 2)[:]
    end

    avg_biomass = vec(mean(total_biomass, dims = 2))
    std_biomass = vec(std(total_biomass, dims = 2))
    upper = avg_biomass .+ std_biomass
    lower = avg_biomass .- std_biomass
    lower = max.(lower, 0.0)

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Average Total Biomass Over Time (± 1 SD)",
        xlabel = "Time",
        ylabel = "Average Total Biomass")

    band!(ax, time, lower, upper, color = (:blue, 0.2))
    lines!(ax, time, avg_biomass, color = :blue, linewidth = 2)

    return fig
end

"""
    plot_avg_species_richness(sol, sites, species; kwargs)

Plot average species richness across all sites over time with ±1 std band.

# Arguments
- `sol`: DifferentialEquations solution object
- `sites`: Vector of site codes
- `species`: Vector of species codes
- `threshold`: Population threshold for species presence (default: 0.1)

# Returns
- A Figure object showing average species richness trajectory with uncertainty band
"""
function plot_avg_species_richness(sol, sites, species;
    threshold::Float64 = 0.1,
    figure_size::Tuple = (900, 600))

    n_sites = length(sites)
    n_species = length(species)
    n_timepoints = length(sol.t)

    richness = zeros(n_timepoints, n_sites)

    for t in 1:n_timepoints
        u_matrix = reshape(sol.u[t], n_sites, n_species)
        for site_idx in 1:n_sites
            richness[t, site_idx] = sum(u_matrix[site_idx, :] .> threshold)
        end
    end

    avg_richness = vec(mean(richness, dims = 2))
    std_richness = vec(std(richness, dims = 2))
    upper = avg_richness .+ std_richness
    lower = max.(avg_richness .- std_richness, 0.0)

    time = sol.t

    fig = Figure(figure_size = figure_size)
    ax = Axis(fig[1, 1],
        title = "Average Species Richness Over Time (± 1 SD, threshold = $threshold)",
        xlabel = "Time",
        ylabel = "Average Number of Species")

    band!(ax, time, lower, upper, color = (:green, 0.2))
    lines!(ax, time, avg_richness, color = :green, linewidth = 2)

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
