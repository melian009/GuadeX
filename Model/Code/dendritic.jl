#----------------------------------------------------------------------------
#Coexistence of native-exotic metacommunities in disturbed dendritic networks
#Melian@EAWAG Mar 2018 -- Guadalquivir basin dendritic data
#Melian@UCordoba April 2023 -- Guadalquivir latitudinal-plots
#Ali advise on structure April 2023
#----------------------------------------------------------------------------
using DataFrames
using CSV
using VegaLite

# PRESENCE MATRIX PER SITE (NOT STANDARIZED BY SAMPLING EFFORT)
cooccur = CSV.read("FishSizeMatrix.csv", DataFrame)
colnames = names(cooccur)

# LON LAT PER SITE
utm = CSV.read("ConnectivityUTM.csv", DataFrame)
colnames = names(utm)

# BUILD REGIONALIZATION (MANY SECTORS)
# https://dataframes.juliadata.org/stable/man/joins/
df2_2 = utm[:, [:CODIGO, :UTMX, :UTMY]]
merged_df = leftjoin(cooccur, df2_2, on = :CODIGO)
colnames = names(merged_df)

#LOOP SECTORS
MINX, MAXX = extrema(merged_df[:, :UTMX])
MINY, MAXY = extrema(merged_df[:, :UTMY])

function equally_spaced_coordinates(lon_min, lon_max, lat_min, lat_max, divisions)
    lon_step = (lon_max - lon_min) / divisions
    lat_step = (lat_max - lat_min) / divisions

    lons = [lon_min + i * lon_step for i in 0:divisions]
    lats = [lat_min + i * lat_step for i in 0:divisions]

    return lons, lats
end

divisions = 10
lons, lats = equally_spaced_coordinates(MINX, MAXX, MINY, MAXY, divisions)


total_unique_native_species = length(unique(merged_df[merged_df.ESPECIE_1 .== "Nativo", :ESPECIE]))

outputdf = DataFrame()
for r in 1:divisions
    rows = (lons[r] .<= merged_df.UTMX .< lons[r+1]) .&& (merged_df.ESPECIE_1 .== "Nativo")
    newdf = merged_df[rows, :]
    newdfgroup = groupby(newdf, :CODIGO)
    newdfcombined = combine(newdfgroup, :ESPECIE => (x -> length(unique(x)) / total_unique_native_species))
    finaldf = unique(leftjoin(newdfcombined, newdf[:, [:CODIGO, :UTMX, :UTMY]], on = :CODIGO))
    finaldf[:, :range] .= r
    outputdf = vcat(outputdf, finaldf)
end


# PLOT

p0 = outputdf |> @vlplot(
    :area,
    transform = [{density = Symbol("ESPECIE_function"), groupby = [:range], extent = [0, 0.65]}],
    x = {field = :value, type = :quantitative,
        axis = {title = "Fraction of native species"},
    },
    y = {field = :density, type = :quantitative, stack = :zero,
        axis = {title = "Density"},
    },
    color = {field = :range, type = :nominal, title = "Longitude range"},
)

save("stacked_density_plot_of_fraction_of_native_specity_per_longitude_section.pdf", p0)

p = outputdf |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)


p2 = outputdf[outputdf.range .== 1, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [180_000, 220_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p3 = outputdf[outputdf.range .== 2, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [215_000, 260_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p4 = outputdf[outputdf.range .== 3, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [250_000, 300_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p5 = outputdf[outputdf.range .== 4, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [290_000, 340_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)


p6 = outputdf[outputdf.range .== 5, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [330_000, 390_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p7 = outputdf[outputdf.range .== 6, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [380_000, 430_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p8 = outputdf[outputdf.range .== 7, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [420_000, 470_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p9 = outputdf[outputdf.range .== 8, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [460_000, 510_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p9 = outputdf[outputdf.range .== 9, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [500_000, 540_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)

p10 = outputdf[outputdf.range .== 10, :] |> @vlplot(
    :point,
    x = {field = :UTMX, type = :quantitative, scale = {type = :sqrt, domain = [530_000, 580_000]}},
    y = {field = Symbol("ESPECIE_function")},
    color = {field = :range, type = :nominal},
    shape = {field = :range, type = :nominal}
)


