# =============================================================================
# Shared output-root resolution for the GuadeX report plotting scripts.
#
# Each plotting script reads the output root of one or more completed sweeps.
# The roots are overridable from the environment so the same scripts can plot
# either the legacy (pre-correction) results or the corrected re-run without
# editing the script body:
#
#   GUADEX_CLIMATE_OUTPUT_DIR      default results/climate_scenarios_k1x_burnin
#   GUADEX_SENSITIVITY_OUTPUT_DIR  default results/sensitivity_obstacles
#   GUADEX_ALT_OUTPUT_DIR          default results/sensitivity_obstacles/alt_interactions
#
# A relative value (including every default) is resolved against `base` when one
# is supplied (the repository root for callers in `scripts/`), otherwise against
# the working directory, exactly reproducing the previous hard-coded constants.
# =============================================================================

"""
    report_output_root(envvar, default; base="")

Return the output root named by the environment variable `envvar`, or `default`
when the variable is unset or empty.  Relative paths are joined to `base` when
`base` is non-empty and otherwise left relative to the working directory.
"""
function report_output_root(envvar::AbstractString, default::AbstractString;
        base::AbstractString="")
    value = strip(get(ENV, envvar, ""))
    isempty(value) && (value = default)
    isabspath(value) && return normpath(value)
    return isempty(base) ? normpath(value) : normpath(joinpath(base, value))
end
