# =============================================================================
# Full GCM x SSP ensemble at effective carrying capacity = 1x observed density,
# per-site daily (seasonal) forcing, and a long baseline burn-in before any
# scenario starts.
#
# Effective K multiplier = carrying_capacity_base_scaling * carrying_capacity_scaling
#                        = 1.0 * 1.0 = 1x observed.
#
# Run from the repository root:
#   pwsh -NoProfile -File scripts/run_k1x_ensemble.ps1
# =============================================================================

$ErrorActionPreference = "Stop"
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)

$env:GUADEX_CLIMATE_OUTPUT_DIR          = "results/climate_scenarios_k1x_burnin"
$env:GUADEX_CLIMATE_FORCING_MODE        = "daily"
$env:GUADEX_CLIMATE_DAILY_FILE          = "guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide.csv"
$env:GUADEX_CLIMATE_DAILY_PER_GCM       = "1"
$env:GUADEX_CLIMATE_SPIN_UP             = "1"
$env:GUADEX_CLIMATE_SPIN_UP_YEARS       = if ($env:GUADEX_CLIMATE_SPIN_UP_YEARS) { $env:GUADEX_CLIMATE_SPIN_UP_YEARS } else { "1000" }
$env:GUADEX_CLIMATE_SPIN_UP_TOL         = if ($env:GUADEX_CLIMATE_SPIN_UP_TOL) { $env:GUADEX_CLIMATE_SPIN_UP_TOL } else { "5.0e-5" }
$env:GUADEX_CLIMATE_SPIN_UP_CRITERION   = if ($env:GUADEX_CLIMATE_SPIN_UP_CRITERION) { $env:GUADEX_CLIMATE_SPIN_UP_CRITERION } else { "basin" }
$env:GUADEX_CLIMATE_K_BASE              = "1.0"
$env:GUADEX_CLIMATE_K_SCALING           = "1.0"
# Re-run any existing outputs so a config change cannot silently resume stale runs.
$env:GUADEX_CLIMATE_FORCE               = if ($env:GUADEX_CLIMATE_FORCE) { $env:GUADEX_CLIMATE_FORCE } else { "1" }
$env:GUADEX_CLIMATE_HEAT_STRESS         = "1"
$env:GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE = "1"
$env:GUADEX_CLIMATE_CONTROL             = "1"
$env:GUADEX_CLIMATE_PLOT                = "1"

Write-Host "Ensemble: effective K = $($env:GUADEX_CLIMATE_K_BASE) x $($env:GUADEX_CLIMATE_K_SCALING) observed; forcing=$($env:GUADEX_CLIMATE_FORCING_MODE); burn-in=$($env:GUADEX_CLIMATE_SPIN_UP_YEARS) yr"
julia --project=. run_climate_scenarios.jl
