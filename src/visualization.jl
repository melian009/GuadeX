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
        limits = (extrema(lons) .+ (-0.5, 0.5), extrema(lats) .+ (-0.5, 0.5)))

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
        limits = (extrema(lons) .+ (-0.5, 0.5), extrema(lats) .+ (-0.5, 0.5)))

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
        limits = (extrema(all_lons) .+ (-0.5, 0.5), extrema(all_lats) .+ (-0.5, 0.5)))

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
        limits = (extrema(lons) .+ (-0.5, 0.5), extrema(lats) .+ (-0.5, 0.5)))

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
- `size`: Tuple for image size (default: (1200, 900))
"""
function save_figure(fig, filename::String; size::Tuple = (1200, 900))
    Makie.save(filename, fig, size = size)
    println("Figure saved to: $filename")
end

# =============================================================================
# --- Sensitivity Visualization: Internal Helpers ---
# =============================================================================

function _compute_richness_matrix(sol, n_sites, n_species, species_indices; threshold=0.1)
    n_t = length(sol.t)
    richness = zeros(Float64, n_sites, n_t)
    for t in 1:n_t
        u_mat = reshape(sol.u[t], n_sites, n_species)
        for i in 1:n_sites
            richness[i, t] = sum(u_mat[i, species_indices] .> threshold)
        end
    end
    return richness
end

function _compute_richness_changes(sol, n_sites, n_species, gidx, row_to_sol, year0_idx, years; days_per_year=365)
    year0_mat = reshape(sol.u[year0_idx], n_sites, n_species)
    results = Vector{Float64}[]
    for yr in years
        target_idx = year0_idx + round(Int, yr * days_per_year)
        yr_mat = reshape(sol.u[target_idx], n_sites, n_species)
        changes = Float64[]
        for sol_i in row_to_sol
            if sol_i === nothing
                push!(changes, 0.0)
            else
                r0 = sum(year0_mat[sol_i, gidx] .> 0.1)
                ry = sum(yr_mat[sol_i, gidx] .> 0.1)
                push!(changes, ry - r0)
            end
        end
        push!(results, changes)
    end
    return results
end

function _diverging_cmap()
    return cgrad([RGBf(0.129, 0.4, 0.675), RGBf(0.97, 0.97, 0.97), RGBf(0.698, 0.094, 0.169)])
end

# =============================================================================
# --- Sensitivity Visualization: Figure 1A — Per-Site Richness Change Map ---
# =============================================================================

"""
    plot_richness_change_per_site(sol, species, sites, site_df, native_idx, invasive_idx, years, output_path; days_per_year=365)

Generate a per-site map of species richness change relative to year 0.
Layout: 2 rows (Invasive / Native) × N columns (one per comparison year).
Each site is colored by Δ richness using a Blue-White-Red diverging scale.
"""
function plot_richness_change_per_site(sol, species, sites, site_df, native_idx, invasive_idx, years, output_path; days_per_year=365)
    n_sites = length(sites)
    n_species = length(species)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)
    row_to_sol = [findfirst(==(String(row.CODIGO)), sites) for row in eachrow(filtered_df)]

    lats = Float64[]
    lons = Float64[]
    for (ex, ey) in zip(Float64.(filtered_df.UTMX), Float64.(filtered_df.UTMY))
        lat, lon = utm_to_latlon(ex, ey, 30)
        push!(lats, lat)
        push!(lons, lon)
    end

    year0_idx = 1
    groups = [("Invasive Species", invasive_idx), ("Native Species", native_idx)]

    all_changes = Float64[]
    group_data = Dict{Int, Vector{Vector{Float64}}}()
    for (gi, (_, gidx)) in enumerate(groups)
        chgs = _compute_richness_changes(sol, n_sites, n_species, gidx, row_to_sol, year0_idx, years; days_per_year=days_per_year)
        group_data[gi] = chgs
        for c in chgs
            append!(all_changes, c)
        end
    end

    max_abs = maximum(abs.(all_changes))
    max_abs = max_abs == 0.0 ? 1.0 : max_abs
    clims = (-max_abs, max_abs)
    diverging_cmap = _diverging_cmap()

    n_rows = 2
    n_cols = length(years)
    cb_width = 70
    fig = Figure(size = (380 * n_cols + cb_width, 380 * n_rows))

    for (gi, (gname, _)) in enumerate(groups)
        chgs = group_data[gi]
        for (ci, yr) in enumerate(years)
            ax = GeoAxis(fig[gi, ci];
                title = "$gname — Year $yr",
                dest = "+proj=latlong",
                limits = (extrema(lons) .+ (-0.5, 0.5), extrema(lats) .+ (-0.5, 0.5)))

            land = GeoMakie.land()
            poly!(ax, land; color = (:lightgray, 0.2), strokecolor = :gray, strokewidth = 0.3)
            lines!(ax, GeoMakie.coastlines(); color = :darkgray, linewidth = 0.5)

            scatter!(ax, lons, lats;
                color = chgs[ci],
                colormap = diverging_cmap,
                colorrange = clims,
                markersize = 10,
                strokecolor = :black,
                strokewidth = 0.5)
        end
    end

    Colorbar(fig[1:2, n_cols + 1], limits = clims, colormap = diverging_cmap,
        label = "Δ Richness", vertical = true, width = 25)

    Label(fig[0, :], "Per-Site Richness Change from Year 0", fontsize = 14, font = :bold)

    save_figure(fig, output_path)
end

# =============================================================================
# --- Sensitivity Visualization: Figure 1B — Per-Subcatchment Richness Change Map ---
# =============================================================================

"""
    plot_richness_change_per_subcatchment(sol, species, sites, site_df, native_idx, invasive_idx, years, output_path; days_per_year=365)

Generate a per-subcatchment map of mean species richness change relative to year 0.
Sites are aggregated to their subcatchment (mean Δ richness) and plotted at the
subcatchment centroid. Layout: 2 rows (Invasive / Native) × N columns (years).
"""
function plot_richness_change_per_subcatchment(sol, species, sites, site_df, native_idx, invasive_idx, years, output_path; days_per_year=365)
    n_sites = length(sites)
    n_species = length(species)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)
    row_to_sol = [findfirst(==(String(row.CODIGO)), sites) for row in eachrow(filtered_df)]

    sc_groups = Dict{String, Vector{Int}}()
    for (fi, row) in enumerate(eachrow(filtered_df))
        sc = string(row.CODIGO_S)
        if !haskey(sc_groups, sc)
            sc_groups[sc] = Int[]
        end
        push!(sc_groups[sc], fi)
    end

    sc_list = sort(collect(keys(sc_groups)))
    sc_lats = Float64[]
    sc_lons = Float64[]
    for sc in sc_list
        fis = sc_groups[sc]
        ux = mean(Float64.(filtered_df.UTMX[fis]))
        uy = mean(Float64.(filtered_df.UTMY[fis]))
        lat, lon = utm_to_latlon(ux, uy, 30)
        push!(sc_lats, lat)
        push!(sc_lons, lon)
    end

    year0_idx = 1
    groups = [("Invasive Species", invasive_idx), ("Native Species", native_idx)]

    all_changes = Float64[]
    group_data = Dict{Int, Vector{Float64}}()
    for (gi, (_, gidx)) in enumerate(groups)
        sc_changes = Float64[]
        for yr in years
            target_idx = year0_idx + round(Int, yr * days_per_year)
            u0_mat = reshape(sol.u[year0_idx], n_sites, n_species)
            uy_mat = reshape(sol.u[target_idx], n_sites, n_species)
            for sc in sc_list
                chgs = Float64[]
                for fi in sc_groups[sc]
                    sol_i = row_to_sol[fi]
                    if sol_i !== nothing
                        r0 = sum(u0_mat[sol_i, gidx] .> 0.1)
                        ry = sum(uy_mat[sol_i, gidx] .> 0.1)
                        push!(chgs, ry - r0)
                    end
                end
                push!(sc_changes, isempty(chgs) ? 0.0 : mean(chgs))
            end
        end
        group_data[gi] = sc_changes
        append!(all_changes, sc_changes)
    end

    max_abs = maximum(abs.(all_changes))
    max_abs = max_abs == 0.0 ? 1.0 : max_abs
    clims = (-max_abs, max_abs)
    diverging_cmap = _diverging_cmap()

    n_rows = 2
    n_cols = length(years)
    cb_width = 70
    n_sc = length(sc_list)
    fig = Figure(size = (380 * n_cols + cb_width, 380 * n_rows))

    for (gi, (gname, _)) in enumerate(groups)
        sc_changes_all = group_data[gi]
        for (ci, yr) in enumerate(years)
            start_idx = (ci - 1) * n_sc + 1
            sc_changes = sc_changes_all[start_idx:start_idx + n_sc - 1]

            ax = GeoAxis(fig[gi, ci];
                title = "$gname — Year $yr",
                dest = "+proj=latlong",
                limits = (extrema(sc_lons) .+ (-0.5, 0.5), extrema(sc_lats) .+ (-0.5, 0.5)))

            land = GeoMakie.land()
            poly!(ax, land; color = (:lightgray, 0.2), strokecolor = :gray, strokewidth = 0.3)
            lines!(ax, GeoMakie.coastlines(); color = :darkgray, linewidth = 0.5)

            scatter!(ax, sc_lons, sc_lats;
                color = sc_changes,
                colormap = diverging_cmap,
                colorrange = clims,
                markersize = 14,
                strokecolor = :black,
                strokewidth = 1)

            for (si, sc) in enumerate(sc_list)
                text!(ax, sc_lons[si], sc_lats[si], text = sc,
                    fontsize = 6, align = (:center, :bottom), offset = (0, 8))
            end
        end
    end

    Colorbar(fig[1:2, n_cols + 1], limits = clims, colormap = diverging_cmap,
        label = "Mean Δ Richness", vertical = true, width = 25)

    Label(fig[0, :], "Per-Subcatchment Richness Change from Year 0", fontsize = 14, font = :bold)

    save_figure(fig, output_path)
end

# =============================================================================
# --- Sensitivity Visualization: Figure 2 — Richness Timeseries 2×2 Grid ---
# =============================================================================

"""
    plot_richness_timeseries_grid(sol, species, sites, site_df, native_idx, invasive_idx, n_sites, n_species, output_path; days_per_year=365)

Generate a 2×2 grid of species richness time series:
  Top-left:  Invasive species, per-site (thin semi-transparent lines + mean ± SD band)
  Top-right: Native species, per-site
  Bottom-left:  Invasive species, per-subcatchment (colored lines + legend + mean overlay)
  Bottom-right: Native species, per-subcatchment
Time axis is in years with a horizontal dashed line at t₀ mean richness.
"""
function plot_richness_timeseries_grid(sol, species, sites, site_df, native_idx, invasive_idx, n_sites, n_species, output_path; days_per_year=365)
    time = sol.t ./ days_per_year
    n_t = length(sol.t)

    filtered_df = filter(row -> row.CODIGO in sites, site_df)
    row_to_sol = [findfirst(==(String(row.CODIGO)), sites) for row in eachrow(filtered_df)]
    n_filtered = nrow(filtered_df)

    sc_map = Dict{String, Vector{Int}}()
    for (fi, row) in enumerate(eachrow(filtered_df))
        sc = string(row.CODIGO_S)
        if !haskey(sc_map, sc)
            sc_map[sc] = Int[]
        end
        push!(sc_map[sc], fi)
    end
    sc_list = sort(collect(keys(sc_map)))
    n_sc = length(sc_list)

    quad_cmap = Makie.resample_cmap(:tab20, max(n_sc, 4))

    native_mat = _compute_richness_matrix(sol, n_sites, n_species, native_idx)
    inv_mat = _compute_richness_matrix(sol, n_sites, n_species, invasive_idx)

    t0_native_mean = mean(native_mat[:, 1])
    t0_inv_mean = mean(inv_mat[:, 1])

    fig = Figure(size = (1400, 1000))

    # --- 2A: Per-Site Invasive ---
    axA = Axis(fig[1, 1];
        title = "Invasive Species Richness (Per-Site)",
        xlabel = "Time (years)", ylabel = "Species Richness (N)")
    for fi in 1:n_filtered
        sol_i = row_to_sol[fi]
        if sol_i !== nothing
            lines!(axA, time, inv_mat[sol_i, :]; color = (:red, 0.2), linewidth = 0.8)
        end
    end
    site_mean = vec(mean(inv_mat, dims = 1))
    site_std = vec(std(inv_mat, dims = 1))
    upper = site_mean .+ site_std
    lower = max.(site_mean .- site_std, 0.0)
    band!(axA, time, lower, upper; color = (:red, 0.15))
    lines!(axA, time, site_mean; color = :darkred, linewidth = 2.5)
    hlines!(axA, [t0_inv_mean]; color = :black, linestyle = :dash, linewidth = 1)

    # --- 2B: Per-Site Native ---
    axB = Axis(fig[1, 2];
        title = "Native Species Richness (Per-Site)",
        xlabel = "Time (years)", ylabel = "Species Richness (N)")
    for fi in 1:n_filtered
        sol_i = row_to_sol[fi]
        if sol_i !== nothing
            lines!(axB, time, native_mat[sol_i, :]; color = (:blue, 0.2), linewidth = 0.8)
        end
    end
    site_mean_n = vec(mean(native_mat, dims = 1))
    site_std_n = vec(std(native_mat, dims = 1))
    upper_n = site_mean_n .+ site_std_n
    lower_n = max.(site_mean_n .- site_std_n, 0.0)
    band!(axB, time, lower_n, upper_n; color = (:blue, 0.15))
    lines!(axB, time, site_mean_n; color = :darkblue, linewidth = 2.5)
    hlines!(axB, [t0_native_mean]; color = :black, linestyle = :dash, linewidth = 1)

    # --- 2C: Per-Subcatchment Invasive ---
    axC = Axis(fig[2, 1];
        title = "Invasive Species Richness (Per-Subcatchment)",
        xlabel = "Time (years)", ylabel = "Species Richness (N)")
    for (si, sc) in enumerate(sc_list)
        fis = sc_map[sc]
        vals = zeros(n_t)
        for t in 1:n_t
            v = Float64[]
            for fi in fis
                sol_i = row_to_sol[fi]
                if sol_i !== nothing
                    push!(v, inv_mat[sol_i, t])
                end
            end
            vals[t] = isempty(v) ? 0.0 : mean(v)
        end
        lines!(axC, time, vals; color = quad_cmap[si], linewidth = 1.5, label = sc)
    end
    lines!(axC, time, site_mean; color = :black, linewidth = 2.5, linestyle = :dash)
    hlines!(axC, [t0_inv_mean]; color = :black, linestyle = :dot, linewidth = 1)
    axislegend(axC; position = :rt, fontsize = 6, nbanks = 3)

    # --- 2D: Per-Subcatchment Native ---
    axD = Axis(fig[2, 2];
        title = "Native Species Richness (Per-Subcatchment)",
        xlabel = "Time (years)", ylabel = "Species Richness (N)")
    for (si, sc) in enumerate(sc_list)
        fis = sc_map[sc]
        vals = zeros(n_t)
        for t in 1:n_t
            v = Float64[]
            for fi in fis
                sol_i = row_to_sol[fi]
                if sol_i !== nothing
                    push!(v, native_mat[sol_i, t])
                end
            end
            vals[t] = isempty(v) ? 0.0 : mean(v)
        end
        lines!(axD, time, vals; color = quad_cmap[si], linewidth = 1.5, label = sc)
    end
    lines!(axD, time, site_mean_n; color = :black, linewidth = 2.5, linestyle = :dash)
    hlines!(axD, [t0_native_mean]; color = :black, linestyle = :dot, linewidth = 1)
    axislegend(axD; position = :rt, fontsize = 6, nbanks = 3)

    Label(fig[0, :], "Species Richness Time Series", fontsize = 16, font = :bold)

    save_figure(fig, output_path)
end
