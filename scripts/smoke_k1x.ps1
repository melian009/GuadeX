# Fast plumbing check for the K=1x / daily-forcing / burn-in pipeline.
# Short burn-in and a single run; skip figures.  Not a scientific run.
$ErrorActionPreference = "Stop"
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)

$env:GUADEX_CLIMATE_OUTPUT_DIR          = "results/_smoke_k1x"
$env:GUADEX_CLIMATE_FORCING_MODE        = "daily"
$env:GUADEX_CLIMATE_DAILY_FILE          = "guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide.csv"
$env:GUADEX_CLIMATE_DAILY_PER_GCM       = "1"
$env:GUADEX_CLIMATE_SPIN_UP             = "1"
$env:GUADEX_CLIMATE_SPIN_UP_YEARS       = "3"
$env:GUADEX_CLIMATE_SPIN_UP_TOL         = "1.0e-3"
$env:GUADEX_CLIMATE_SPIN_UP_CRITERION   = "basin"
$env:GUADEX_CLIMATE_K_BASE              = "1.0"
$env:GUADEX_CLIMATE_K_SCALING           = "1.0"
$env:GUADEX_CLIMATE_HEAT_STRESS         = "1"
$env:GUADEX_CLIMATE_HEAT_STRESS_CALIBRATE = "1"
$env:GUADEX_CLIMATE_CONTROL             = "1"
$env:GUADEX_CLIMATE_MAX_RUNS            = "1"
$env:GUADEX_CLIMATE_PLOT                = "0"

julia --project=. run_climate_scenarios.jl
