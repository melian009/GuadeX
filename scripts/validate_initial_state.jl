# =============================================================================
# validate_initial_state.jl
#
# Post-spin-up validation of the interim (C3) projection route against the
# 2006-2009 survey (issue 4.6).
#
# The interim route saves the spun-up state as `baseline_density` in each run's
# `export/levels/species_timeseries.csv`.  This script loads the no-warming
# control run (so the comparison is against the control's own spun-up state, not
# a scenario), compares presence/absence with the observed survey per species,
# and reports the control change 2026 -> 2045 per species.
#
# Usage:
#   julia --project=. scripts/validate_initial_state.jl
#
# Env overrides:
#   GUADEX_VALIDATION_RUN   run directory (default the corrected control run)
#   GUADEX_VALIDATION_OUT   output CSV (default under results/)
#   GUADEX_VALIDATION_THRESHOLD  model presence threshold (default 0.1)
# =============================================================================

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Guadex
using DataFrames
using CSV
using Statistics

const ROOT = dirname(@__DIR__)
const RUN_DIR = get(ENV, "GUADEX_VALIDATION_RUN",
    joinpath(ROOT, "results", "climate_scenarios_corrected", "control", "baseline"))
const OUT_DIR = get(ENV, "GUADEX_VALIDATION_OUT",
    joinpath(ROOT, "results", "climate_scenarios_corrected", "validation"))
const THRESHOLD = parse(Float64, get(ENV, "GUADEX_VALIDATION_THRESHOLD", "0.1"))
const DENSITY_FILE = joinpath(ROOT, "data", "BIOTIC", "FishDensity_and_Juveniles_Matrix.csv")
const SPECIES_TS = joinpath(RUN_DIR, "export", "levels", "species_timeseries.csv")

isfile(SPECIES_TS) || error("missing control time series: $SPECIES_TS")

ts = CSV.read(SPECIES_TS, DataFrame)
years = sort(unique(ts.year))
first_year, last_year = years[1], years[end]

# Site order from the run, so observed and modelled matrices align.
sites = unique(String.(ts.CODIGO))
species = sort(unique(String.(ts.species)))

density_df, _ = Guadex.load_species_density_data(DENSITY_FILE)
observed = Guadex.observed_density_matrix(density_df, sites, species)

baseline_rows = ts[ts.year .== first_year, :]
base = zeros(length(sites), length(species))
site_index = Dict(site => i for (i, site) in enumerate(sites))
species_index = Dict(sp => i for (i, sp) in enumerate(species))
for r in eachrow(baseline_rows)
    i = site_index[String(r.CODIGO)]
    s = species_index[String(r.species)]
    base[i, s] = Float64(r.baseline_density)
end

# Average ranks (ties) for a Spearman correlation.
function _ranks(x::AbstractVector)
    n = length(x)
    order = sortperm(x)
    r = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && x[order[j + 1]] == x[order[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            r[order[k]] = avg
        end
        i = j + 1
    end
    return r
end

function _spearman(a::AbstractVector, b::AbstractVector)
    length(a) < 3 && return NaN
    return cor(_ranks(a), _ranks(b))
end

# Control change 2026 -> 2045 per species.
first_rows = ts[ts.year .== first_year, :]
last_rows = ts[ts.year .== last_year, :]
first_total = combine(groupby(first_rows, :species), :density => sum => :total)
last_total = combine(groupby(last_rows, :species), :density => sum => :total)
start_total = Dict(String(r.species) => Float64(r.total) for r in eachrow(first_total))
final_total = Dict(String(r.species) => Float64(r.total) for r in eachrow(last_total))
change = Dict{String,NamedTuple{(:start, :final, :relative),Tuple{Float64,Float64,Float64}}}()
for sp in species
    start = get(start_total, sp, 0.0)
    final = get(final_total, sp, 0.0)
    change[sp] = (start=start, final=final,
                  relative=start > 0 ? (final - start) / start : NaN)
end

rows = NamedTuple[]
for (s, sp) in enumerate(species)
    obs_present = observed[:, s] .> 0.0
    mod_present = base[:, s] .> THRESHOLD
    tp = count(obs_present .& mod_present)
    fn = count(obs_present .& .!mod_present)
    fp = count(.!obs_present .& mod_present)
    tn = count(.!obs_present .& .!mod_present)
    sensitivity = (tp + fn) > 0 ? tp / (tp + fn) : NaN
    specificity = (tn + fp) > 0 ? tn / (tn + fp) : NaN
    tss = (isnan(sensitivity) || isnan(specificity)) ? NaN : sensitivity + specificity - 1
    mask = obs_present
    rho = count(mask) >= 3 ? _spearman(observed[mask, s], base[mask, s]) : NaN
    c = change[sp]
    push!(rows, (
        species=sp,
        n_sites=length(sites),
        observed_occupancy=count(obs_present) / length(sites),
        modelled_occupancy=count(mod_present) / length(sites),
        tp=tp, fn=fn, fp=fp, tn=tn,
        sensitivity=sensitivity, specificity=specificity, tss=tss,
        abundance_rank_rho=rho,
        control_start=Float64(c.start), control_final=Float64(c.final),
        control_relative_change=c.relative,
    ))
end

result = DataFrame(rows)
mkpath(OUT_DIR)
out_csv = joinpath(OUT_DIR, "initial_state_validation.csv")
CSV.write(out_csv, result)

println("Validation of the post-spin-up state: $(RUN_DIR)")
println("  schema: observed presence = survey density > 0; modelled presence = spun-up density > $THRESHOLD")
println("  sites = $(length(sites)), species = $(length(species)), years = $first_year-$last_year")
println("  pooled sensitivity = ",
    round(sum(result.tp) / max(sum(result.tp) + sum(result.fn), 1), digits = 3),
    ", specificity = ",
    round(sum(result.tn) / max(sum(result.tn) + sum(result.fp), 1), digits = 3))
println()
println("species  obs_occ  mod_occ   TSS    rho    dN(26-45)")
for r in eachrow(result)
    println(rpad(r.species, 8), lpad(round(r.observed_occupancy, digits = 3), 8),
        lpad(round(r.modelled_occupancy, digits = 3), 9),
        lpad(round(r.tss, digits = 3), 7),
        lpad(round(r.abundance_rank_rho, digits = 2), 7),
        lpad(isnan(r.control_relative_change) ? "n/a" :
             string(round(100 * r.control_relative_change, digits = 1), "%"), 11))
end
println("\nWrote $out_csv")
