using Guadex
using DataFrames

println("Checking site ordering and dam indexing...")

data = prepare_ode_data()

# Check site ordering
sites = data.sites
density_df = data.density_df
density_df_filtered = filter(row -> row.CODIGO in sites, density_df)

println("Number of sites in connectivity: ", length(sites))
println("Number of sites in filtered density: ", nrow(density_df_filtered))

if length(sites) == nrow(density_df_filtered)
    matches = all(sites .== density_df_filtered.CODIGO)
    println("Sites match exactly in order: ", matches)
    if !matches
        println("First 5 sites in connectivity: ", sites[1:5])
        println("First 5 sites in filtered density: ", density_df_filtered.CODIGO[1:5])
    end
else
    println("Mismatch in number of sites!")
end

# Check dam indexing
dams = data.dams
n_sites = length(sites)
symmetric = all(dams[i, j] == dams[j, i] for i in 1:n_sites, j in 1:n_sites)
println("Dam matrix is symmetric: ", symmetric)

if !symmetric
    # Find a non-symmetric pair
    for j in 1:n_sites, i in 1:n_sites
        if dams[i, j] != dams[j, i]
            println("Non-symmetric pair found: sites[", i, "]=", sites[i], ", sites[", j, "]=", sites[j])
            println("dams[", i, ", ", j, "] (passability ", sites[j], " -> ", sites[i], ") = ", dams[i, j])
            println("dams[", j, ", ", i, "] (passability ", sites[i], " -> ", sites[j], ") = ", dams[j, i])
            break
        end
    end
end
