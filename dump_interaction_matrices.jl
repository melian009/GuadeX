using Pkg; Pkg.activate(".");
using DataFrames, CSV, Guadex, Random

# =============================================================================
# Generate markdown dumps of the three interaction matrices
# =============================================================================

const NATIVE_SPECIES = Set(["AB", "AH", "SP", "PW", "LS", "SA", "IL", "CP", "IO", "ST", "LR", "MC"])
const INVASIVE_SPECIES = Set(["GH", "MS", "LG", "CC", "CG", "AM", "OM", "EL", "GL", "TT", "AAL"])

function classify(sp::String)
    sp in NATIVE_SPECIES && return "Native"
    sp in INVASIVE_SPECIES && return "Invasive"
    return "Other"
end

function make_random_interaction_matrix(n_species::Int, val_range::Float64=1.0)
    Random.seed!(42)
    mat = zeros(n_species, n_species)
    for i in 1:n_species, j in 1:n_species
        i == j && continue
        mat[i, j] = -rand() * val_range
    end
    return mat
end

function make_invasive_favoring_matrix(n_species::Int, native_idx::Vector{Int}, invasive_idx::Vector{Int})
    mat = zeros(n_species, n_species)
    for inv in invasive_idx
        for nat in native_idx
            mat[nat, inv] = -0.8
        end
    end
    for inv in invasive_idx
        for inv2 in invasive_idx
            inv == inv2 && continue
            mat[inv2, inv] = -0.3
        end
    end
    for nat in native_idx
        for nat2 in native_idx
            nat == nat2 && continue
            mat[nat2, nat] = -0.1
        end
    end
    return mat
end

function count_nonzero(mat)
    count(!iszero, mat)
end

function write_matrix_md(filename::String, label::String, description::String,
    mat::Matrix{Float64}, species::Vector{String})
    open(filename, "w") do io
        n = length(species)
        nz = count_nonzero(mat)
        total_cells = n * n - n  # excluding diagonal
        nnz_diag = count(i -> !iszero(mat[i,i]), 1:n)

        write(io, "# Interaction Matrix: $label\n\n")
        write(io, "$description\n\n")
        write(io, "**Dimensions:** $(n) × $(n) species  \n")
        write(io, "**Non-zero off-diagonal entries:** $nz / $total_cells ($(round(nz/total_cells*100, digits=1))%)  \n")
        write(io, "**Convention:** `mat[recipient, actor]` = effect of actor on recipient  \n")
        write(io, "**Diagonal:** $(nnz_diag > 0 ? "non-zero ($nnz_diag entries)" : "zero (self-regulation via logistic carrying capacity)")  \n\n")

        write(io, "## Species Legend\n\n")
        write(io, "| Code | Status |\n")
        write(io, "|------|--------|\n")
        for sp in species
            write(io, "| $sp | $(classify(sp)) |\n")
        end
        write(io, "\n## Interaction Matrix\n\n")
        write(io, "Columns = actors, Rows = recipients. Negative = suppression/predation, zero = neutral.\n\n")
        write(io, "```\n")

        # Header row with species codes + classification
        write(io, "     ")
        for j in 1:n
            c = classify(species[j])
            abbr = c[1:1]  # N, I, O
            write(io, rpad("$(species[j])($abbr)", 8))
        end
        write(io, "\n")

        for i in 1:n
            c = classify(species[i])
            abbr = c[1:1]
            write(io, rpad("$(species[i])($abbr)", 4))
            for j in 1:n
                val = mat[i, j]
                if i == j
                    write(io, rpad(".", 8))
                else
                    write(io, rpad(string(round(val, digits=1)), 8))
                end
            end
            write(io, "\n")
        end
        write(io, "```\n\n")
    end
    println("Wrote $filename")
end

# =============================================================================
# Main
# =============================================================================
data_base = prepare_ode_data(upstream_cost=0.05)
species = data_base.species
n_species = length(species)
println("Species order ($(n_species)): $species")

native_idx = [i for (i, sp) in enumerate(species) if classify(sp) == "Native"]
invasive_idx = [i for (i, sp) in enumerate(species) if classify(sp) == "Invasive"]
println("Native indices: $native_idx")
println("Invasive indices: $invasive_idx")

original_mat = data_base.params.interaction_matrix
random_mat = make_random_interaction_matrix(n_species, abs(minimum(original_mat)))
invasive_fav_mat = make_invasive_favoring_matrix(n_species, native_idx, invasive_idx)

mkpath("docs")
write_matrix_md("docs/interaction_matrix_original.md",
    "Original (Empirical Guadalquivir)",
    "Empirical interaction matrix from the Guadalquivir River basin fish metacommunity. " *
    "Loaded from `data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv`. " *
    "Descriptive interaction strings were converted to numeric scores (−1 = no coexistence, −0.8 = exclusion/predation, −0.5 = competition, −0.3 = mild effect, 0 = neutral).",
    original_mat, species)

write_matrix_md("docs/interaction_matrix_random.md",
    "Random",
    "Fully randomized interaction matrix. Every off-diagonal entry is a random value " *
    "drawn uniformly from [−1.0, 0.0]. Seed fixed to 42 for reproducibility. " *
    "This serves as a null hypothesis: with no ecological structure, what richness pattern emerges?",
    random_mat, species)

write_matrix_md("docs/interaction_matrix_invasive_favoring.md",
    "Invasive-Favoring",
    "Invasives strongly suppress natives (−0.8), natives have zero effect on invasives. " *
    "Invasives compete moderately with each other (−0.3). Natives compete weakly with each other (−0.1). " *
    "Unclassified species (e.g., AA) are neutral (0.0). Diagonal is zero. " *
    "This is the inverse of the empirical pattern observed in the original matrix.",
    invasive_fav_mat, species)
