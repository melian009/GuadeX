# Using the updated 2045 inputs in simulations

This note describes the current, deliberately explicit integration of the CEDEX 2045 tables and the 2026 obstacle inventory. The implementation is in `src/data_preparation.jl`; the simulation entry points use it from `scripts/run_model.jl`, `scripts/run_sensitivity.jl`, `run_sensitivity_report.jl`, and `run_alt_interactions.jl`.

## Run parameters live in `parameters.toml`

All run settings for the four entry scripts — the updated input paths, obstacle overlay settings, species groups, subcatchments, scenario factors, sweep grids, and simulation lengths — are read from the single parameter file **`parameters.toml`** in the repository root. The entry scripts no longer hard-code these values; they load them through **`parameters.jl`** (module `SimulationParameters`) at startup and print which parameter file was used.

* **Change a parameter:** edit `parameters.toml` and rerun.
* **Track an experiment:** copy the file to e.g. `parameters_experiment_a.toml`, edit the copy, and run with `GUADEX_PARAMETERS_FILE=<path>`.
* **Validate a configuration without simulating:** run any entry script with `GUADEX_CONFIG_ONLY=1`; it loads the parameters, prints the provenance line, and exits before data loading or simulation.
* **Per-run overrides:** individual values can still be overridden with `GUADEX_CEDEX_VAR_FILE`, `GUADEX_CEDEX_UTS_FILE`, `GUADEX_OBSTACLES_FILE`, `GUADEX_OBSTACLE_MODE`, `GUADEX_OBSTACLE_TOLERANCE_M`, `GUADEX_OBSTACLE_PASSABILITY`, and `GUADEX_OBSTACLE_DOWNSTREAM_PASSABILITY`.
* Library defaults in `src/data_preparation.jl` are unchanged; the scripts pass the parameter file values explicitly to `prepare_ode_data`.

## What is integrated now

`prepare_ode_data` accepts the optional updated inputs. The scripts pass the values from `parameters.toml`:

```julia
data = prepare_ode_data(
    upstream_cost = 0.05,
    cedex_var_file = "data/2045_CEDEX_GUADALQUIVIR_VAR.csv",
    cedex_esc_uts_file = "data/2045_CEDEX_GUADALQUIVIR_ESC_por_UTS.csv",
    obstacles_file = "data/obstacles_1658_Obstaculos_No_Completamente_Franqueables _2026-02-16_Guadex.csv",
    obstacle_mode = :overlay,
    obstacle_matching_tolerance = 2000.0,
    obstacle_passability = 0.1,
    obstacle_downstream_passability = 0.5,
)
```

* The CEDEX files are loaded, validated, and returned as `data.cedex_var_df` and `data.cedex_esc_uts_df`. The UTS table is tidy and can be queried explicitly:

  ```julia
  row = select_cedex_uts(data.cedex_esc_uts_df;
      uts="ES050_01", scenario="SSP245", measure="percent", season="AMJ")
  ```

* `obstacle_mode = :overlay` maps each obstacle with valid UTM coordinates to the nearest connected site-to-site network segment within the matching tolerance. Matched segments limit movement towards the higher-elevation (upstream) site to `obstacle_passability` (default 0.1) and movement towards the lower-elevation (downstream) site to `obstacle_downstream_passability` (default 0.5): downstream passage is reduced but remains more accessible than upstream passage. Equal-elevation links are not restricted because their flow direction cannot be inferred. The effective matrix is the element-wise minimum of the legacy dam matrix and this directional obstacle overlay.
* Matching results are available in `data.obstacle_mapping_diagnostics`; `matched == false` records obstacles that were not assigned to a network segment. For matched non-flat segments, `restricted_origin` and `restricted_destination` identify the directional movement that receives the upstream (limiting) obstacle passability.
* `data.legacy_dams` is retained for comparison with `data.dams` (the effective matrix). The source `IF` field is not used to infer passability.

With the current repository files, the default 2,000 m matching tolerance matches 973 of 1,658 obstacles. The other 685 records remain in `data.obstacles_df` and in `data.obstacle_mapping_diagnostics`, but are not assigned to a model edge because no connected site segment is within the tolerance. This is a limitation of the spatial crosswalk, not a claim that those obstacles are invalid or absent in the source inventory.

The default library behavior remains `obstacle_mode = :legacy` when no updated files are supplied. This preserves existing callers and tests. The four maintained simulation entry points use the repository files and the `parameters.toml` value (`:overlay` by default). Set `GUADEX_OBSTACLE_MODE=legacy` (or edit `obstacles.mode` in `parameters.toml`) to run a comparison without the new obstacle overlay.

## Important limitation: CEDEX climate is not yet a temperature time series

The `VAR` table contains `PRE`, `ETP`, `ETR`, `REC`, and `ESC`; it does not contain a temperature column. The UTS table contains seasonal runoff changes. Consequently, loading these files does **not** change `params.temperatures` or `params.habitat_suitability`. The existing `temperature_increases` grids in the sensitivity scripts remain the explicit temperature scenarios, and the CEDEX tables are recorded as scenario inputs/provenance until a validated transformation is available.

Do not apply a CEDEX percentage or millimetre value to every site. A UTS-to-`CODIGO`/`CODIGO_S` crosswalk and a hydrological relationship to discharge, habitat, carrying capacity, or mortality are required first.

## Recommended workflow for future runs

1. **Validate the source and crosswalk.** Confirm the meaning of the CEDEX region columns, resolve the `SSP285` versus `SSP585` label discrepancy with the data provider, and add a versioned UTS-to-site or UTS-to-subcatchment crosswalk.
2. **Choose the climate representation.** For a static 2045 sensitivity run, create documented site-level temperature and habitat vectors. For seasonal or year-by-year runs, extend `MetacommunityParams`/`metacommunity_ode!` to read a time-indexed environmental schedule rather than adding one constant delta.
3. **Calibrate hydrological effects.** Define how percentage and millimetre runoff changes modify flow, habitat availability, carrying capacity, or growth. Keep percentage and millimetre measures as separate scenarios and record the selected season.
4. **Validate obstacle matching.** Inspect `obstacle_mapping_diagnostics`, test several tolerances, and compare matched segments against GIS/network geometry. Confirm whether bridges, culverts, weirs, and passages should have different rules.
5. **Replace the provisional passability rule.** Obtain the definition and scale of `IF`, then implement a documented categorical or continuous, potentially species-specific passage function. Until then, `obstacle_passability = 0.1` (upstream) and `obstacle_downstream_passability = 0.5` (downstream) are conservative modelling assumptions, not observations.
6. **Run paired comparisons.** Save legacy and obstacle-overlay runs with the source filenames, tolerance, passability rule, CEDEX scenario, season, and crosswalk version in the output metadata.

For a reproducible baseline comparison, run once with `GUADEX_OBSTACLE_MODE=legacy` and once with the default `overlay`, then compare the resulting dispersal matrices and ecological outputs before interpreting changes as climate effects.
