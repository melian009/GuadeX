# This is not part of the Guadex package, but a separate module for data exploration.

using DataFrames
using CSV


# Get species counts and their native/exotic status from the FishSizeMatrix.csv file
fishsizematrixfile = "data\\FishSizeMatrix.csv"
fsm_df = CSV.read(fishsizematrixfile, DataFrame)
count_df = combine(groupby(fsm_df, [:ESPECIE, :ESPECIE_1]), nrow => :count)

# Get species initials and status dictionary
exotic_species = count_df[count_df.ESPECIE_1.=="Exotica", :]
native_species = count_df[count_df.ESPECIE_1.=="Nativo", :]
function species_initials(species_name::AbstractString)
    words = split(species_name)
    if length(words) >= 2
        return join(uppercase.(first.(words[1:2])))
    else
        return uppercase(first.(words))
    end
end
exotic_species.initials = species_initials.(exotic_species.ESPECIE)
native_species.initials = species_initials.(native_species.ESPECIE)
# Full species names mapping
species_names = Dict(exotic_species.initials[i] => exotic_species.ESPECIE[i] for i in 1:nrow(exotic_species))
species_names = merge(species_names, Dict(native_species.initials[i] => native_species.ESPECIE[i] for i in 1:nrow(native_species)))
# Species status mapping
species_status = Dict(exotic_species.initials[i] => "Exotica" for i in 1:nrow(exotic_species))
species_status = merge(species_status, si => "Nativo" for i in 1:nrow(native_species)))