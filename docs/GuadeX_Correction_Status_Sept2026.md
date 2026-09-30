# GuadeX Corrections — Status and Decision Record

**Status date:** 2026-09-30. **Base commit:** `master` @ `9702169` (all correction work
below is in the **working tree only**; nothing is committed).

## 1. Purpose and how to read

This is the single current status/decision document for the GuadeX review corrections. It
answers three questions:

1. **What is actually done** — which review items (C1–C7, E1–E22, minor issues) are
   implemented and verified in the working tree, and what evidence supports each claim.
2. **What to run next** — the exact runbook for the corrected climate re-run, including the
   endpoint re-selection step and the legacy-configuration trap.
3. **What is still open** — every decision that needs the team and every known residual.

The companion documents are the **assessment** (`docs/GuadeX_Review_Assessment_Sept2026.md`,
which adjudicates the reviewer's claims) and the **plan**
(`docs/GuadeX_Correction_Plan_Sept2026.md`, which orders the work). This document is the
status of that plan. The plan is not rewritten here.

How to read the tables: **Status** is `Implemented` when the change is present and exercised
by the test suite in the working tree. **Evidence / verification** names the code and the
check that supports the claim. Every numeric claim was re-derived from the tree on
2026-09-30 (see §6); where a claim could not be re-derived it is explicitly marked
**unverified** in §5 and §6 rather than silently repeated.

---

## 2. Implemented fixes (working tree, uncommitted)

| ID | What changed | Evidence / verification | Status |
|----|--------------|-------------------------|--------|
| **C2** | Interaction term moved **inside** the intrinsic-growth bracket (competitive Lotka–Volterra, `c = 1 − α`); the legacy outside-growth form is opt-in via `interaction_inside_growth = false`. | `src/ode_model.jl` (`interaction_inside_growth`, `clamp(1 − total/K + interaction, −1, 2)` at ~L144); default `true` in the constructor and docs. `test/test_ode.jl` testset *interaction form: competitive Lotka–Volterra (C2)* (13 tests) checks the two-species equilibrium `K(1−c₁₂)/(1−c₁₂c₂₁)` and per-capita growth scaling. | Implemented |
| **C1** | Interaction direction parsed from the cell text (named target, so `α[target, source]`); undirected competition made symmetric; `No coexist` → no interaction; Spanish `depreda`/`interfiere` handled; mixed-mechanism cells flagged ambiguous; the missing *Lepomis* (`LG`) row/column is now present. New long table `data/BIOTIC/interaction_matrix_long.csv`. | `src/data_preparation.jl` (`parse_interaction_cell`, `build_interaction_long_table`, long-table loader); `scripts/build_interaction_long_table.jl`. Verified: the long table has **227 data rows** and **22 ambiguous source cells** (44 ambiguous rows); the loaded matrix is 24×24 including `LG` and `AAL`. `test/test_data_preparation.jl` (parsing helpers, data-loading testsets). | Implemented |
| **E5** | Reporting offsets are `1:20`, so label `2045` is the **end** of year 2045 (`t = 7300`), not 1 Jan 2045. | `run_climate_scenarios.jl` (`YEAR_OFFSETS = collect(1:SIMULATION_YEARS)`, `SIMULATION_YEARS = END_YEAR − START_YEAR + 1 = 20`). `test/test_outputs.jl` *annual snapshot alignment* (30 tests). | Implemented |
| **E6** | Quasi-extinction is evaluated only over sites where the species was **baseline-established**; a species absent everywhere in the baseline yields fraction `0.0`; the established denominator is exposed. | `src/outputs.jl` (`baseline_established`, `established` argument, `n_sites_established`); tests *quasi_extinction_flags never flags non-established sites (E6)*, *quasi-extinction fraction over established sites (E6)*, *species absent everywhere in the baseline → fraction 0 (E6)*, *trout-style fraction is not 1 − occupancy (E6)*. | Implemented |
| **E9 / E10** | Parameter **and code digest** folded into burn-in cache keys and completed-run fingerprints; an existing run is skipped **only when its stored digest matches**; corrected scenario config and legacy reproduction config written and distinguished. | `sensitivity_core.jl`, `src/run_identity.jl`, `run_climate_scenarios.jl`; `run_metadata.json`/`run_digest.txt` sidecars. Tests *parameter digest*, *burn-in cache key includes the parameter digest*, *run fingerprint includes the parameter digest*, *committed k1x burn-in configuration*, *corrected climate re-run configuration*. | Implemented |
| **E22** | Per-month held-out-basin (LOBO) bias correction applied in the `guadex_tw` temperature pipeline; the 48 model-facing temperature products were rewritten. | `guadex_tw/scripts/13b_apply_bias_correction.py`, `guadex_tw/models/lobo_monthly_bias_correction.json`, `guadex_tw/logs/bias_correction_13b_run.log`. Verified: pooled annual **+1.1059 °C**, Feb **+2.2237**, Mar **+2.2867**; log reports **48 files processed in 699.4 s** (44 per-GCM wide files 8.63 GB + ensemble-median wide + baseline + future-stats + historical ≈ 9 GB). Anomalies unchanged, as designed. `scripts/13_project_guadex_sites.py` also fixed/plumbed (the correction itself lives in `13b`). | Implemented |
| **C5 / E1** | Model site temperature level = corrected `tw_baseline_mean` (1986–2005) selected by **exact column name**, joined by site code, with a loud fallback; `level + anomaly` reconstructs the corrected daily series; k-calibration and exposure use the same array. | `src/data_preparation.jl` (`WATER_TEMPERATURE_COLUMN = "tw_baseline_mean"`, `extract_site_temperatures`), `src/temperature_forcing.jl`. Tests *corrected site water-temperature baseline (C5/E1)*, *model temperature equals corrected daily series (C5/E1)* (16 tests), *driver daily schedules reconstruct the corrected series (ITEM 1)* (18 tests). Corrected **trout-site mean = 12.8305 °C** (`test/test_temperature_forcing.jl`, `isapprox(..., 12.8304986; atol=0.01)`); corrected basin range 9.785–18.064 °C. | Implemented |
| **C7** | Realised applied warming persisted per run (2026–2045 and 2036–2045 means over the sites) and used as the dose–response regressor; exposure restricted to baseline-established sites; GCM-aware dose–response (fixed-effect GCM, no new dependency); endpoint re-selection function `select_climate_endpoints`. | `src/temperature_forcing.jl` (`realised_warming_anomaly`), `src/climate_diagnostics.jl`, `run_climate_scenarios.jl`, `docs/climate_scenarios.md`. Tests *realised warming anomaly (C7)*, *established_exposure_summary (C7)*, *dose_response handles GCM structure (C7)*, *select_climate_endpoints / realised_forcing_axis (C7)*. | Implemented |
| **E4** | Robust composition stop rule: 95th percentile of `|Δ log N|` over **active** site×species cells (density floor 0.1) in at least one compared year, combined with the basin-total rule (`spin_up_criterion = "both"`, `spin_up_composition_tol = 1e-2`, `spin_up_min_years = 10`); the raw all-cells statistic is recorded as a diagnostic; **actual matrix-specific burn-in years** are recorded in metadata. | `src/spin_up.jl` (`SPIN_UP_COMPOSITION_ACTIVE_FLOOR = 0.1`, `composition_q95_change`, `composition_q95_change_all_cells`, `composition_diagnostics`); tests *spin_up convergence rule (E4)* (29 tests) and *burn-in diagnostics flow into run metadata and index*. On real data the robust statistic crosses 1e-2 at ≈ year 423 (< 600 cap) per the config comment; the raw all-cells statistic plateaus near 1.4e-2. | Implemented |
| **C3** | **Interim projection route**: `projection_route = "interim_observed"` starts from the observed community and integrates exactly `interim_spin_up_years = 3` baseline years (`interim_observed_spin_up`), reporting every scenario as per-metric `scenario − control` deltas (`scenario_minus_control` / `scenario_minus_control_table`, written by the figures path). The route is an explicit opt-in; `"full_burnin"` remains the default so the E4 machinery is untouched elsewhere, and the route mandates `control_run = true`. Route and years are recorded in `export/run_metadata.json`. | `src/spin_up.jl` (`projection_route`, `interim_observed_spin_up`), `run_climate_scenarios.jl`, `src/climate_diagnostics.jl`; tests *C3 interim projection route* (spin-up) and *scenario_minus_control (C3 interim route)*, plus corrected-config route assertions. | Implemented (provisional pending calibration) |
| **Consistency** | All three drivers reconstruct the corrected daily series; the standalone exporter uses 1-based offsets; the serializer uses the shared (corrected) extractor; corrected config vs legacy reproduction config are explicitly separated. | `run_climate_scenarios.jl`, `run_climate_experiments.jl`, `run_alt_interactions.jl`, `scripts/export_run_outputs.jl`, `scripts/serialize_data.jl` (imports `Guadex` and calls the shared extractor). Tests *driver daily schedules…*, *export_run_outputs*. | Implemented |
| **Housekeeping** | Initial state aligned by site code; seasonal rate option; vacuous tests replaced and mass-conservation/topology tests added; viewer species-level export + aggregate rendering; ML (`reml=False`) LRT refit in the air–water model. | `run_climate_scenarios.jl`, `src/spin_up.jl`, `test/test_graph_construction.jl`, `test/test_ode.jl` (*dispersal operator: mass conservation and structure*), `test/test_data_preparation.jl`, viewer export under `results/.../export/viewer/species_<sp>_*_timeseries.json`, `guadex_tw/scripts/08_models.py`. | Implemented |
| **E11 / report text** | `IET` renamed to the **bank-stability** index; water-body metadata corrected (**289** reporting units = 288 assigned + 1 unassigned, 774 network links); report and decision brief corrected and marked as pending re-run where numbers will change. | `src/data_preparation.jl` docstring; `legacy/FinalReportSeptember2026.tex`, `legacy/GuadalquivirDecisionBrief.tex` (both also flag `\pendingrerun`). Both TeX documents compiled (29 + 7 pages) in the working tree (`legacy/*.log` present). | Implemented |
| **C4** | New **on-path river-network graph** (fallback, because `SW_Line_4C_.shp` carries no flow-direction attribute or elevation): on-path parent rule plus explicit MST root-join; deterministic tie-breaking. | `src/data_preparation.jl` (`build_on_path_distance_matrix`, defaults `rtol=0.02`, `atol=200 m`); `test/test_on_path_graph.jl` (testset **2647 tests**). Verified on real data: **774 = n−1 edges**, 267 roots joined by 266 MST edges, **total 13,871.5 km** (vs legacy 31,206.3 km), **4** parent elevation inversions vs **238** in the legacy graph, acyclic and deterministic. Default stays `legacy`; the corrected config selects `connectivity.method = "on_path"`. | Implemented (fallback) |
| **E12** | Obstacles snapped to the selected (corrected) tree; passabilities **multiply** along a link (0.1·0.1 = 0.01) instead of `min()`; non-operational/ambiguous structures classified and excluded/annotated; legacy de-duplication; deterministic. `IF` index deliberately **not** used as passability. | `src/data_preparation.jl` (`build_obstacle_passability_matrix`, status classification); tests *param digest* + `test/test_data_preparation.jl` (`n_status_active = 1274`, `n_status_nonoperational = 185`, `n_status_ambiguous = 199`) and `test/test_on_path_graph.jl` overlay counts. Verified overlay line: matched 973/1658 (legacy) or 1036/1658 (on-path), applied 727/796, excluded 387. `IF` non-numeric rows confirmed at **547/1658** in `src/data_preparation.jl` comments — this differs from the review's 445 and is stated as such. | Implemented |
| **E13–E16 / E18** | Six explicit, opt-in biological options added with **defaults unchanged** (byte-identical default regression): `absence_growth_fraction`, `pool_capacity_mode`, `fishless_dificil_capacity`, `nonreproducing_local_growth`, `exclude_fishfarm_eel_records`, `salinity_envelope`. Metadata + digest wired; documented. | `parameters.toml [biological_options]`, `SimulationParameters.biological_options`, `prepare_ode_data`; `test/test_data_preparation.jl` *Explicit biological options (E13-E16, E18)* (46 tests); `test/test_parameters.jl` *biological-assumption option loader*. Documented in `docs/climate_scenarios.md`. | Implemented; the corrected config now selects the team's six values (E13 `1.0`, E14 `"cap"`, E15 `"near_zero"`, E16 `"zero"` + fish-farm drop, E18 `true`), while `parameters.toml` keeps the legacy defaults. |
| **Validation** | Full Julia suite and a real corrected smoke run. | Full suite: **3,693 / 3,693 pass, 67 testsets, 0 failures** (re-run 2026-09-30; see §6). Smoke run present at `results/_smoke_corrected/` covering **4 SSPs + control** (ACCESS-CM2) and exercising the corrected pipeline. Digest-skip behaviour tested (see E9/E10). | Implemented |

---

## 3. Corrected re-run runbook

The corrected re-run **has been executed**: **45 climate runs** (44 GCM x SSP scenarios + a
matched no-warming control), **32 obstacle runs** and **48 interaction runs**. Everything below
is the exact procedure that was followed and is retained as the run record.

### 3.0 The one rule that matters most

> **Do NOT use `legacy/parameters_climate_scenarios_k1x_burnin.toml` for corrected results.**
> That file is the **legacy September reproduction**: it pins `spin_up_criterion = "basin"`,
> `spin_up_min_years = 2`, the pre-E4 stop rule and the pinned September heat-stress slope
> (`k = 3.7817692640400324e-5`, `calibrate = false`). Use it only to reproduce
> `results/climate_scenarios_k1x_burnin/`.

> **Do not circulate the decision brief until the corrected re-run has completed.** The brief
> currently carries a `\pendingrerun` notice; its quantitative claims are provisional.

### 3.1 Configuration

`parameters_climate_scenarios_corrected.toml` (working tree, uncommitted). It differs from
the legacy file only in the corrected settings: E4 `spin_up_criterion = "both"` (in both
`[run_climate_scenarios]` and `[obstacle_sensitivity]`), heat stress derived from the
corrected baseline (`[temperature_stress].calibrate = true`), `connectivity.method =
"on_path"`, and the new output root. The corrected `tw_baseline_mean` source and
`interaction_inside_growth = true` are **code defaults with no config opt-in**, so the file
uses them by construction and cannot select the legacy alternatives. It additionally selects
the team's decisions (§4): the C3 interim projection route
(`projection_route = "interim_observed"`, `interim_spin_up_years = 3`), the six biological
values in `[biological_options]`, and the uniform `IF`/passability assumption noted in
`[obstacles]`.

### 3.2 Command (PowerShell, this platform)

```powershell
$env:GUADEX_PARAMETERS_FILE = 'parameters_climate_scenarios_corrected.toml'
julia --project=. run_climate_scenarios.jl
```

Bash equivalent:

```bash
GUADEX_PARAMETERS_FILE=parameters_climate_scenarios_corrected.toml \
  julia --project=. run_climate_scenarios.jl
```

Matching climate-consistent obstacle/interaction sweeps:

```powershell
$env:GUADEX_PARAMETERS_FILE = 'parameters_climate_scenarios_corrected.toml'
julia --project=. run_sensitivity_report.jl
```

**Before running the sweeps**, override the output directory (see §5):
`$env:GUADEX_SENSITIVITY_OUTPUT_DIR = 'results/sensitivity_obstacles_corrected'`, otherwise
`results/sensitivity_obstacles` is overwritten.

### 3.3 Outputs, time and directory

- **Output root:** `results/climate_scenarios_corrected` (a **new** directory, deliberately
  distinct from `results/climate_scenarios` and `results/climate_scenarios_k1x_burnin`, so
  the pre-correction September artifacts are preserved for comparison).
- **Contents:** 11 GCMs × 4 SSPs + a no-warming control (44 runs + control). Each run writes
  `simulation_output.jld2`, `run_digest.txt`, `export/run_metadata.json` and the four-level
  + viewer exports; `runs_index.csv` is rewritten after every run, so an interrupted ensemble
  resumes without losing finished runs.
- **Expected time:** **~35–45 min** for the 44-run ensemble + control on this machine — this
  is an estimate from the brief, not a measured figure (the real measured smoke run is much
  shorter; see §6). Budget extra for the objective/level figure generation
  (`make_figures = true` by default; disable with `GUADEX_CLIMATE_PLOT=0`).

### 3.4 Digest behaviour

A completed run is **reused only when its stored digest matches** the current one. The digest
covers the resolved model state (interaction matrix, growth/temperature/K/thermal arrays,
heat-stress slope, dispersal matrix and scaling, `interaction_inside_growth`), the resolved
run options (forcing mode, baseline period, burn-in years/tolerance/criterion, K scaling,
upstream cost, obstacle/passability configuration, scenario, GCM, warming curve) and the
**code version** (`git rev-parse --short HEAD`). Outputs written before this change carry no
digest and are re-run once. Use `$env:GUADEX_CLIMATE_FORCE = '1'` to recompute regardless.

### 3.5 Endpoint re-selection (do this **after** the run)

The two climate endpoints used as the background climate for the obstacle/upstream-cost
sweeps must be re-selected from the **realised** forcing of the corrected ensemble, not
hard-coded:

```julia
using DataFrames, CSV, Guadex
idx = CSV.read(joinpath("results", "climate_scenarios_corrected", "runs_index.csv"), DataFrame)
endpoints = Guadex.select_climate_endpoints(idx)
println(endpoints.source)              # "realised" once the columns are populated
println(endpoints.cool.scenario, "/", endpoints.cool.gcm)
println(endpoints.warm.scenario, "/", endpoints.warm.gcm)
```

Then place those `scenario`/`gcm` pairs into `[obstacle_sensitivity].climate_models` (or the
`GUADEX_SENSITIVITY_CLIMATE_MODELS` override; use `scenario_col = :climate_scenario` for a
sensitivity index). Run the sweeps only after this step.

---

## 4. Team decisions

The review items that needed a team choice are listed below. The decisions taken for the
**corrected re-run** (2026-09-30) are marked **DECIDED** and implemented in
`parameters_climate_scenarios_corrected.toml`; the remaining items are still open. The
corrected results under the interim route and the chosen biological options remain
**provisional pending full calibration**.

### Decision record for the corrected re-run (2026-09-30)

| decision | chosen value | where |
| :--- | :--- | :--- |
| **C3** projection route | **interim observed start + short spin-up, scenario − control** (`projection_route = "interim_observed"`, `interim_spin_up_years = 3`) | `[run_climate_scenarios]`; `Guadex.interim_observed_spin_up` (`src/spin_up.jl`); deltas via `Guadex.scenario_minus_control` (`src/climate_diagnostics.jl`) |
| **E13** `absence_growth_fraction` | `1.0` (temperature/habitat decide establishment) | `[biological_options]` |
| **E14** `pool_capacity_mode` | `"cap"` (84 isolated pools clipped at the non-pool median) | `[biological_options]` |
| **E15** `fishless_dificil_capacity` | `"near_zero"` (56 `SIN_PECES = DIFICIL` sites transit-only) | `[biological_options]` |
| **E16** `nonreproducing_local_growth` | `"zero"` (`AA`, `LR`, `MC` local r = 0) | `[biological_options]` |
| **E16** `exclude_fishfarm_eel_records` | `true` (drop the four Guadiato fish-farm eel records) | `[biological_options]` |
| **E18** `salinity_envelope` | `true` (strictly brackish `AB`/`LR`/`MC` below the conductivity threshold; no elevation envelope) | `[biological_options]` |
| **IF passability** | keep the uniform structural assumption (`0.1` upstream / `0.5` downstream); do not decode `IF` | `[obstacles]`; note in `docs/climate_scenarios.md` |

`parameters.toml` **keeps the safe legacy defaults** for all six biological options and does
not select the C3 route, so the legacy reproduction and the default test suite are unchanged;
the route is an explicit opt-in and the E4 convergence machinery still applies to every
configuration that does not select it.

> **Consolidated October-report limitations.** The standing caveats of the team decisions and
> the other current model limitations are consolidated for the forthcoming October 2026 report
> in **`docs/GuadeX_Limitations_and_Scope_Oct2026.md`**. That document records the reflected
> team decisions on **C6** (thermal niches; item 1), **C3** (interim route; item 2), **E21**
> (air–water slope uncertainty band; item 3), the obstacle **IF**/passability assumption
> (item 4) and **E8** (deterministic realised richness loss; item 5). The report text is drawn
> from that document rather than from this status record.

1. **C3 — Calibrate before projecting (ABC-SMC). DECIDED: interim route.**
   *Options:* (a) run simulation-based inference on `r`, `K`, interaction coefficients and
   dispersal with **pre-registered** TSS / abundance-rank / site-richness criteria, spatial
   block CV and posterior predictive checks; (b) the interim route the review allows — start
   from the observed state with a short spin-up and report scenario − control.
   *Decision:* **(b)**, implemented as `projection_route = "interim_observed"` with a 3-year
   fixed baseline spin-up from the observed community and a **mandatory** matched control
   (`control_run = true`). Results are reported as per-metric `scenario − control` deltas. (a)
   is deferred: the **corrected run is therefore provisional, not a calibrated forecast**; the
   acceptance criteria and compute for (a) remain open.
2. **E8 — Stochasticity / extinction probability.**
   *Options:* add demographic/environmental stochasticity with replicates and
   species-specific quasi-extinction thresholds, and rename the current
   richness-loss metric to `realised_richness_loss`; or keep the deterministic model and stop calling the
   realised-richness-loss metric an "extinction probability".
   *Consequence:* adding replicates multiplies the run cost; keeping determinism requires a
   naming/documentation fix so the metric is not over-claimed. Design decision.
3. **C6 — Physiological thermal niches.**
   *Options:* replace range midpoints with per-species optimum-for-growth and CTmax/UILT
   from the literature, on an asymmetric (Sharp–Schoolfield-style) curve, propagating ±2 °C;
   or keep the symmetric empirical-range Gaussian.
   *Consequence:* needs **literature values per species** (not in the repo). Matters mainly
   under stronger warming (all baseline sites ≤ 18.37 °C).
4. **E21 — Lagged air–water refit.**
   *Options:* refit with lagged/distributed-lag air temperature and propagate the slope
   range; or treat the lagged model as an uncertainty analysis only.
   *Consequence:* the repo has only **3 in-sample lag-eligible sites**, so a lagged refit
   cannot become a new primary calibration; the honest use is an uncertainty band.
5. **E7 — Global Sobol sensitivity.**
   *Options:* replace the cross-design factor ranking with within-design sums of squares
   **plus** a global design (Latin hypercube + Sobol indices) over agreed ranges; or keep the
   current design-limited ranking and label it as such.
   *Consequence:* needs **agreed parameter ranges and dependencies**; otherwise the indices
   are uninterpretable.
6. **E13–E16 / E18 — the six biological options. DECIDED for the corrected re-run.** The
   team chose `absence_growth_fraction = 1.0`, `pool_capacity_mode = "cap"`,
   `fishless_dificil_capacity = "near_zero"`, `nonreproducing_local_growth = "zero"`,
   `exclude_fishfarm_eel_records = true`, `salinity_envelope = true` (decision table above;
   set in `[biological_options]` of the corrected config). `parameters.toml` keeps the legacy
   defaults. The consequences below remain the standing caveats of those choices:
   - **E13 `absence_growth_fraction`** — `0.1` (legacy) vs `1.0` (temperature/habitat decide)
     vs `0.0` (forbid local growth when absent). *Consequence:* raising it makes rare/invasive
     spread easier; it **interacts with E15**.
   - **E14 `pool_capacity_mode`** — `"legacy"` vs `"cap"` vs `"exclude"`. *Consequence:* the
     84 isolated-pool sites hold ~56% of density and drive the K tail (max K ≈ 95,333 → ≈
     8,057 under either alternative). The `cap` rule uses the **non-pool median** (a choice).
   - **E15 `fishless_dificil_capacity`** — `"legacy"` vs `"near_zero"` (K = 1e-6 at the 56
     `SIN_PECES = DIFICIL` sites). *Consequence:* makes them transit-only nodes; **interacts
     with E13**.
   - **E16 `nonreproducing_local_growth`** — `"legacy"` vs `"zero"` for `AA`, `LR`, `MC`.
     *Consequence:* `zero` removes local recruitment; without an estuarine-immigration term
     the species only decays locally.
   - **E16 `exclude_fishfarm_eel_records`** — `false` vs `true` (drops the four documented
     Guadiato fish-farm eel records; AA presences 9 → 5).
   - **E18 `salinity_envelope`** — `false` vs `true` (zeros strictly brackish species at sites
     below 1000 µS/cm `CONDUCTIVIDAD`; 1,584 cells). A static **elevation** envelope is
     deliberately **not** offered (it would block upslope range shifts).
7. **Alien species introductions** — which species to include and **where/when** they are
   introduced. *Consequence:* the current model has no explicit introduction events; adding
   them changes invasive-richness projections.
8. **Obstacle `IF` passability-index semantics. DECIDED: keep the uniform assumption.**
   Barrier passability is currently modelled using a uniform structural assumption (0.1
   upstream, 0.5 downstream) pending full standardisation of the regional obstacle
   inventory's qualitative passability index (IF). The `IF` index is therefore **not**
   decoded for the corrected re-run; the `547/1658` non-numeric rows (the repo's count; the
   review said 445) are not yet standardised. *Consequence:* using `IF` prematurely would
   silently mis-scale every barrier, so the uniform `0.1`/`0.5` assumption is retained. This
   note is recorded in `[obstacles]` of the corrected config and in
   `docs/climate_scenarios.md`.

---

## 5. Known residuals and limitations

1. **C4 is a documented fallback, not the GIS build.** The repo's `SW_Line_4C_.shp` has no
   flow-direction attribute or elevation, so the corrected graph uses **chord snapping**
   (on-path parent rule + MST root join), not flow-directed geometry. Root-join edges are
   **flow-undirected** (41 inversions in the root-join set; only 4 parent-edge inversions).
   The `rtol`/`atol` slack is the key modelling choice — the root count is sensitive to it
   (defaults `rtol = 0.02`, `atol = 200 m`, 267 roots).
2. **Obstacle `IF` semantics unresolved; the `ESTADO` codebook is absent.** Structures are
   classified active (1274) / non-operational (185) / ambiguous (199) from the codes that are
   present; unknown/blank codes are reported as ambiguous rather than dropped.
3. **E14 pool cap rule uses the non-pool median** — a defensible but arbitrary choice.
4. **The corrected re-run has been executed** (see §6): the corrected outputs are in
   `results/climate_scenarios_corrected` (45 rows: 44 GCM x SSP scenarios + control),
   `results/sensitivity_obstacles_corrected` (32 obstacle rows) and
   `results/alt_interactions_corrected` (48 interaction rows). **Endpoints must be re-selected
   on the realised forcing after the run** (§3.5). Under the C3 interim route the run starts
   from the observed community with a 3-year baseline spin-up, so the results are
   **provisional and uncalibrated**; they are reported as `scenario − control` deltas and must
   not be read as a calibrated forecast.
5. **Pre-correction documents moved to `legacy/`:** the stale-statement documents
   `legacy/integrated_climate_obstacle_report.md`,
   `legacy/publication_manuscript_climate_fragmentation.md` and
   `legacy/ReportingMay2026.tex` carry pre-correction statements and have been moved to
   `legacy/` (together with `legacy/FinalReportSeptember2026.tex` and
   `legacy/GuadalquivirDecisionBrief.tex`), so they cannot be mistaken for current results.
   The legacy `k1x` result directories remain in place, marked by a `LEGACY_README.md`.
6. **`[obstacle_sensitivity].output_dir` in the default and legacy configs still points at
   `results/sensitivity_obstacles`** (`parameters.toml`,
   `legacy/parameters_climate_scenarios_k1x_burnin.toml`), so a sweep run with those configs
   would **overwrite** the pre-correction sweep. The corrected config
   (`parameters_climate_scenarios_corrected.toml`) writes to the corrected directories; a
   default run can still be redirected with `GUADEX_SENSITIVITY_OUTPUT_DIR`.
7. **`scripts/serialize_data.jl` still has a stale `SPECIES_CODES`** — a duplicated `"AA"`
   and no `"AAL"` (L31–34). The corrected species list is resolved elsewhere; this dead
   constant should be deleted or fixed.
8. **The seasonal rate option is not exposed via `parameters.toml`** (code-only).
9. **Ranking of what the corrected run will change is provisional.** The report and brief
   correctly mark the affected numbers as pending; the qualitative directions (warming
   effect, trout decline, invasive-richness value) may still move.

---

## 6. Verification summary

- **Working tree:** `git status` shows no staged changes and no new commits; 47 tracked files
  modified and the new files untracked (including this document, the assessment, the plan,
  `parameters_climate_scenarios_corrected.toml`, `legacy/parameters_climate_scenarios_k1x_burnin.toml`,
  `data/BIOTIC/interaction_matrix_long.csv`, `src/run_identity.jl`, the new tests, the
  `guadex_tw` bias-correction script/model/logs, and the `.biascorr.json` sidecars).
- **Full Julia test suite (re-run 2026-09-30):** **3,693 / 3,693 pass, 67 testsets, 0
  failures.** The run includes the C4 on-path graph testset (2,647 tests), the C2
  interaction-form testset, the E4 spin-up-rule testset, the C5/E1 reconstruction testsets,
  the C7 testsets, the E6 quasi-extinction testsets, the E13–E16/E18 option testsets, and the
  new C3 testsets (*C3 interim projection route*, *scenario_minus_control (C3 interim
  route)*) plus the corrected-config route/biological-value assertions.
- **Real smoke run:** `results/_smoke_corrected/` contains a completed run for **4 SSPs +
  control** (ACCESS-CM2) with full exports, exercising the corrected daily-forcing /
  `tw_baseline_mean` / interaction / digest pipeline. (A **218 s** wall time was quoted for
  the brief; no timing is persisted in the tree, so that figure is **unverified**.)
- **Invariants checked:**
  - `level + anomaly` reconstructs the corrected daily series (tests assert to `atol = 1e-9`
    and, for the driver schedules, exactly). A specific **8.9e-16** residual was quoted but is
    **not recorded** in the tree; the invariant itself is tested.
  - Corrected trout-site mean temperature = **12.8305 °C**.
  - On-path tree: **774 edges, 13,871.5 km, 4 parent inversions vs 238** in the legacy graph;
    acyclic and deterministic.
  - Obstacle overlay: passabilities multiply; excluded-structure counts consistent.
  - Digest-skip: an existing run with a mismatched digest is recomputed (tested).
  - Default biological options reproduce the historical arrays (regression test).

### Facts that could not be verified against the tree

- **ML (`reml=False`) LRT figures for the air–water `ta:z` interaction.** The code change is
  present (`guadex_tw/scripts/08_models.py` refits both models by maximum likelihood when
  `reml=False`), and the qualitative direction is consistent — but the quoted values
  (Guadalquivir ≈ **16.08, p ≈ 6e-5**; Spain ≈ **41.45, p ≈ 1e-10**) are **not persisted**.
  The committed `guadex_tw/models/model_coefficients.json` stores the **REML** statistics
  (Guadalquivir **11.96, p = 5.42e-4**; Spain **34.94, p = 3.40e-9**), and
  `guadex_tw/logs/acceptance.json` quotes the Spain REML value. Treat the ML figures as
  unconfirmed until `08_models.py` is re-run and its output committed.
- **Smoke-run wall time (218 s)** and **temperature reconstruction residual (8.9e-16)** — not
  recorded anywhere in the tree.
- **E4 real-data convergence year (~423)** — asserted in the corrected config and plan
  comments; the test suite exercises the rule with synthetic data, not the real burn-in.

---

## 7. Pointers and commit status

- Consolidated limitations and scope of inference for the October 2026 report:
  `docs/GuadeX_Limitations_and_Scope_Oct2026.md` (reflects the C6, C3, E21, IF and E8 team
  decisions).
- Assessment (adjudication of the review): `docs/GuadeX_Review_Assessment_Sept2026.md`.
- Plan (ordered work: Waves 0–3, decisions, gates): `docs/GuadeX_Correction_Plan_Sept2026.md`.
- This status document: `docs/GuadeX_Correction_Status_Sept2026.md`.
- Corrected config: `parameters_climate_scenarios_corrected.toml`.
  Legacy reproduction (do not use for corrected results):
  `legacy/parameters_climate_scenarios_k1x_burnin.toml`.
- Runbook detail and environment overrides: `docs/climate_scenarios.md`.

**Nothing is committed.** Base commit is `master` @ `9702169`; all corrections, new tests,
configs, data products and this document are in the working tree only.
