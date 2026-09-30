# GuadeX Correction Plan

> **Status (2026-09-30):** current implementation state, the corrected re-run runbook, open team decisions and known residuals are recorded in **`docs/GuadeX_Correction_Status_Sept2026.md`** — read that for status; this document remains the work plan.

**Companion to `docs/GuadeX_Review_Assessment_Sept2026.md`.** This plan addresses only the
review items judged **FAIR** or **PARTLY FAIR** in that assessment, and explicitly corrects
the reviewer's proposed remedies where they are wrong. Work is grouped so that each wave can
be implemented and tested together, and so that expensive model re-runs happen as few times
as possible.

Repository state: `master`, commit `9702169`. Working tree additionally contains the two
new documents (`GuadeX_Review_and_Corrections_Sept2026.md`, this plan, and the assessment).

---

## 0. Ground rules and constraints

1. **One combined re-run per wave is cheaper than many.** The fixes interact: fixing the
   interaction direction (C1) and scaling (C2) without fixing the temperature level (C5) and
   the basin bias (E22) would produce a new but still uninterpretable result. Group them.
2. **The 44-run ensemble forcing is not in Git.** The daily per-GCM wide files
   (`guadex_tw/outputs/tables/water_temp_daily_guadex_sites_wide_ssp*_*.csv`, 11 GCMs × 4 SSPs)
   are present on this checkout but gitignored (`guadex_tw/.gitignore:15`), so Wave 1 is
   re-runnable here; a fresh clone would require the team's forcing archive. Compute remains
   the binding constraint for the calibration and replicate runs.
3. **Do not apply the review literally.** The specific corrections are listed in §1
   ("rejected/replaced").
4. **Add a test with every behavioural fix.** The repository currently has no test for the
   sign of α, network topology, mass conservation, or the temperature column (minor #13);
   these tests are part of each fix, not an afterthought.
5. **Record everything in `run_metadata.json`** (resolved settings, parameter digest, code
   commit) so the E9/E10 reproducibility failures cannot recur.

### 0.1 Reviewer recommendations that must be rejected or replaced

| Reviewer proposal | Replace with |
|---|---|
| C2: `c_sj = −α_sj` | `c_sj = 1 − α_sj`, `c_ss = 1` (matches the reviewer's own equilibrium check) |
| C1: "iterate over all rows" to include *Lepomis* | Add an `Lg` column / recode to long format; iterating rows cannot supply the missing column |
| C4: dendritic-parents tree "≈ 8,640 km, close to MST" | Add an explicit root-joining rule; the complete tree is ≈ 13,871 km on the current data, so set the correctness criterion on topology and on-path consistency, not on matching the MST weight |
| Part II: obstacle passability from `IF/10` | Do not implement until IF semantics are confirmed; 445/1,658 rows are non-numeric |
| E18: add an elevation envelope | Only add a salinity/brackish envelope; the reviewer's own Part II notes that a static elevation envelope would block upslope range shifts |
| E20: make `06a` outlier flag two-sided | `06a` is already two-sided (`abs(resid) > 4·MAD`); fix the real problem — apply the existing flags in `08_models.prep()` |
| C2 Part I wording: interaction "1–3 orders of magnitude" | Keep (verified; up to ~3.5 orders at the low-`r` end), but describe as a documentation/calibration issue rather than a hidden bug, since Equation 1 documents the form |

---

## 1. Wave 0 — Safe, local fixes (no new model science, can ship with unit tests)

These change code/labels/tests without changing any parameter estimate. Implement together
and validate with the Julia test suite.

### 1.1 C1 — Interaction direction, one-way loading, and the missing *Lepomis* column
- **Files:** `data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv`,
  `src/data_preparation.jl:load_interaction_matrix` (564), `parse_interaction_string` (631),
  `test/test_data_preparation.jl:187`.
- **Do:** recode the matrix as an explicit long table `(source, target, mechanism, α,
  ambiguous, reference)`; parse the named target from the text and store `α[target, source]`;
  make undirected competition symmetric; treat `No coexist` as **absence of interaction**
  (allopatry is already handled by the thermal/habitat filters) — this is a team-confirmable
  biological decision, flagged in §5; add Spanish `depreda`/`interfiere`; add the missing
  `Lg` column; flag mixed-mechanism cells (e.g. `"interfere through competition and Ms
  depreda Gh"`) rather than silently symmetrising them.
- **Do not:** rely on "iterating over all rows" (does not fix *Lepomis*), and use one single
  regex consistently (the review shows two different `TARGET_RX` patterns).
- **Tests:** `α[SA, MS] < 0` and `α[MS, SA] == 0`; summed invasive→native effect < 0; every
  species has both a row and a column; no `Lg` row/column is silently dropped; update the
  test at line 187 that currently pins the inverted value.
- **Effort:** M. **Blocking:** team validation of the long table (§5).

### 1.2 C2 — Put the interaction term inside the growth term
- **File:** `src/ode_model.jl:104, 113–115, 128`.
- **Do:** `dU[i,s] = N_is * (r_eff * clamp(1 - total_i/K_i + interaction_term, -1, 2) -
  heat_stress)`, where `interaction_term` remains `Σ_j α_sj N_j /K_i` (it is already divided
  by `K`). Keep the old form behind a flag for comparison.
- **Tests:** one site, two species, no dispersal → equilibrium matches
  `K(1 − c₁₂)/(1 − c₁₂c₂₁)` with `c = 1 − α`; per-capita growth doubles when `r` doubles.
- **Effort:** S. Note: the clamp now also bounds the interaction; document the clipping
  (minor #10).

### 1.3 E5 — Report end-of-year states
- **Files:** `run_climate_scenarios.jl:144`, `run_climate_experiments.jl:43`,
  `sensitivity_core.jl:135`.
- **Do:** offsets `1:20` (labels stay `2026:2045`). Verify `_snapshot_index` /
  `searchsortedfirst` alignment (`src/outputs.jl:249`) before trusting the month.
- **Tests:** first reported row ≠ initial state; exposure/quasi-extinction series shift
  consistently.
- **Effort:** S.

### 1.4 E6 — Quasi-extinction only where a species was established
- **File:** `src/outputs.jl:267, 425, 546`.
- **Do:** never flag a site with baseline density ≤ threshold; add a `baseline_established`
  mask/argument to `quasi_extinction_summary` and the site table, and divide by the
  established count. (The summary function currently has no such argument — this is the
  missing piece the review omits.)
- **Tests:** trout fraction no longer equals `1 − occupancy`.
- **Effort:** S.

### 1.5 E9/E10 — Reproducible runs and configuration
- **Files:** `sensitivity_core.jl:307–348, 506–525`; `run_climate_scenarios.jl:358–372`;
  `parameters.toml`.
- **Do:** extend the existing FNV-1a digest to include interaction/growth/temperature/`K`/
  thermal arrays and the git commit; add a `[run_climate_scenarios]` block reproducing the
  September settings (`forcing_mode="daily"`, `spin_up=true`, 373 yr,
  `carrying_capacity_base_scaling=1`, heat stress on, control run); write all resolved
  settings to `run_metadata.json`; make the existing-file skip compare digests.
- **Effort:** M.

### 1.6 E11 / report text — Labels and statements
- **Files:** `src/data_preparation.jl:962–989` (IET docstring), `legacy/FinalReportSeptember2026.tex`.
- **Do:** rename `IET` to the bank-stability index throughout code/comments/report; correct
  the report statements in the assessment §"Report statements" (water bodies 289; temperature
  level; habitat label; 16.0 °C; same-anomaly text; dispersal multiples; README counts), and
  rename the `crosswalk_water_bodies` field so it counts sites, not water bodies.
- **Effort:** S.

### 1.7 Minor housekeeping
- Align initial states by site code in **all three** run scripts
  (`run_climate_scenarios.jl:201`, `run_climate_experiments.jl:102`,
  `run_alt_interactions.jl:152`) and use `innerjoin(...; order=:left)`.
- Document/repair the rate conversion (minor #4): either correct the *Gambusia* comment or
  offer `log1p(r)/365`.
- Replace vacuous tests (`test_graph_construction.jl:9,16`, `test_data_preparation.jl:402`,
  `test_visualization.jl`) and add the missing mass-conservation, topology, α-sign and
  temperature-column tests.
- Remove the stale `SPECIES_CODES` duplicate `AA` / missing `AAL` or delete the dead constant.
- Viewer: export one species-level series (data already computed in `species_metrics`) and
  render sub-catchment/water-body/basin aggregates.
- **Effort:** S–M.

---

## 2. Wave 1 — One coordinated re-run of the climate ensemble

Run after Wave 0, and **only** with all of these in place, because each depends on the
others.

### 2.1 C5 + E1 + E22 — A single, corrected water-temperature series
- **Files:** `src/data_preparation.jl:extract_site_temperatures` (922) and the run scripts
  that build the parameter set (`run_climate_scenarios.jl:389–401`);
  `guadex_tw/scripts/13_project_guadex_sites.py` (206–209).
- **Do:** (a) add the monthly held-out-basin residual correction from E22 to script 13;
  (b) plumb the corrected `tw_baseline_mean` (1986–2005) into the model as the site level,
  selected by exact column name, joining on site code and handling the 776-row/775-site
  mismatch and `1.30.20`; (c) fail loudly if the expected column is missing; (d) ensure the
  heat-stress calibration `k` and the exposure diagnostics use the same corrected series.
- **Test:** with this level, model temperature equals the daily guadex_tw series; mean trout
  site temperature ≈ 11.7 °C; exposure and `k` use the same array.
- **Effort:** M. **Note:** this is the dependency the reviewer's C5 fix omits.

### 2.2 C7 — Realised forcing and exposure over occupied sites
- **Files:** `run_climate_scenarios.jl` (run index, lines 366/489), `src/outputs.jl:973` and
  `exposure_table` (~407).
- **Do:** persist the realised mean applied anomaly (2036–2045 and 2026–2045 over the 775
  sites) for every run; use it as the dose–response regressor; fit with a random intercept
  per GCM and report model spread; re-select the cool/warm endpoints on the realised forcing;
  compute exposure only over baseline-established sites for each species (mean/median, not
  basin maximum).
- **Tests:** run index contains the realised anomaly; exposure figure names its site set.
- **Effort:** M.

### 2.3 E2/E3 — Like-for-like controls and burn-ins
- **Files:** `run_climate_scenarios.jl`, `sensitivity_core.jl`, `run_sensitivity_report.jl:88–141`.
- **Do:** add a no-warming control per GCM using that GCM's own historical series, and report
  scenario − control per GCM; for the obstacle sweep, either run a burn-in per
  (upstream-cost, passability) combination or pair each run with a matched control.
- **Effort:** M (compute-heavy). **Dependency:** per-GCM historical series must exist in the
  forcing files.

### 2.4 E4 — Stronger convergence criterion
- **File:** `src/spin_up.jl:121–134`.
- **Do:** stop on both total biomass and a species×site criterion (e.g. 95th percentile of
  |Δ log N|) plus a minimum number of years, and record both in metadata. Record the
  **actual** matrix-specific burn-in years (the report and review both wrongly say 21).
- **Effort:** S–M.

### 2.5 Report corrections implied by the re-run
- Update all dependent numbers (trout slope/R², exposure days, invasive richness, richness
  loss) and state the survey dates (2006–2009).

---

## 3. Wave 2 — Network and obstacles (moderate; separable)

### 3.1 C4 — Rebuild the river network from the river lines
- **New routine**, using `data/Version_02-01-2026-Masas/SW_Line_4C_.shp` (in-repo) with
  geopandas: snap each site to its segment, derive flow direction, connect each site to its
  nearest downstream neighbour along the channel; break ties deterministically.
- **Fallback** only if the GIS route is not taken: implement the on-path parent rule **and**
  an explicit root-joining rule, and validate by topology (`n−1` edges, acyclic, on-path,
  no elevation inversions beyond documented cases). Do **not** claim "close to the MST":
  on current data the complete fallback tree is ≈ 13,871 km.
- **Tests:** `n−1` edges; every edge on a flow path; deterministic tie-breaking; link length
  reported but not used as the pass criterion.
- **Effort:** L.

### 3.2 E12 — Obstacles on the corrected network
- **File:** `src/data_preparation.jl:build_obstacle_passability_matrix` (456),
  `build_dam_passability_matrix` (849).
- **Do:** re-snap obstacles to the corrected network; accumulate passabilities along a link
  (currently `min()`, so multiple barriers do not compound); exclude/annotate demolished and
  out-of-basin structures; de-duplicate against the legacy dam layer.
- **Do not:** use the `IF` index as passability until its scale/direction is confirmed and
  the 445 non-numeric rows are handled.
- **Effort:** M–L.

---

## 4. Wave 3 — Calibration and inputs (research / team decisions)

These are the fixes the review correctly marks as depending on the team's decisions, data,
or compute. They should be scoped but not started until §5 is settled.

### 4.1 C3 — Calibrate before projecting
- Implement simulation-based inference (ABC-SMC) on `r`, `K`, interaction coefficients and
  dispersal, with pre-registered acceptance criteria (per-species TSS, abundance rank
  correlation, site-richness distribution), spatial block cross-validation and posterior
  predictive checks.
- **Interim:** start projections from the observed state with a short spin-up and report
  scenario-minus-control, as the review suggests.
- **Effort:** L. **Dependency:** compute and team decisions; forcing files are present locally.

### 4.2 E8 — Stochasticity and extinction probability
- Add demographic/environmental stochasticity with replicates and species-specific
  quasi-extinction thresholds; rename the current richness-loss metric to `realised_richness_loss`.
- **Effort:** L.

### 4.3 C6 — Physiological thermal niches
- Replace range midpoints with per-species optimum-for-growth and CTmax/UILT from the
  literature, with an asymmetric (Sharp–Schoolfield) curve; propagate uncertainty (±2 °C)
  and sensitivity-test.
- **Note:** with all baseline sites ≤ 18.37 °C, asymmetry matters mainly under stronger
  warming. **Effort:** M–L (literature-dependent).

### 4.4 E13–E19 — Input and biological options
- E13: make the 10% absence-growth factor an explicit parameter (review suggests 1.0 — team
  decision; note it interacts with E15).
- E14: exclude/cap isolated pools when estimating `K`; standardise densities by season;
  document units.
- E15: allow `K ≈ 0` at `SIN_PECES = DIFÍCIL` sites.
- E16: set local `r = 0` for non-reproducing species and model them as estuarine immigration;
  remove the four Guadiato fish-farm eel records.
- E17: rebuild growth/dispersal from verifiable Iberian sources; single table including
  *Alburnus*.
- E18: add a **salinity** envelope only (not elevation).
- E19: represent dry reaches as seasonal barriers/transit nodes and scale `K`/connectivity
  with the CEDEX runoff projections.
- **Effort:** L, mostly data/decision-bound.

### 4.5 E21 — Air–water model
- Refit with lagged/distributed-lag air temperature and propagate the slope range. The repo
  has only 3 in-sample lag-eligible sites, so this must be treated as an uncertainty
  analysis, not a new primary calibration.
- **Effort:** M.

### 4.6 E7 — Global sensitivity
- Replace the cross-design factor ranking with within-design sums of squares plus a global
  design (Latin hypercube + Sobol indices) over agreed parameter ranges.
- **Effort:** L.

---

## 5. Decisions required from the team

1. **Interaction table (C1):** validate every pairwise entry, especially predation direction
   and the mixed-mechanism cells; confirm whether `No coexist` should be "no interaction"
   (recommended) or retained as a negative effect.
2. **Biological options (E13–E16):** default value of the absence-growth factor; whether to
   zero `K` at `DIFÍCIL` sites; whether to model non-reproducing species as estuarine
   immigrants; whether to remove the fish-farm eel records.
3. **Thermal niches (C6) and growth/dispersal rates (E17):** adopted values and uncertainty
   ranges, with references.
4. **Calibration acceptance criteria (C3):** thresholds fixed before calibrating.
5. **Alien species not in the 2006–2009 survey:** which to include, and where/when
   introduced.
6. **Resources:** access to the per-GCM daily forcing files and compute for the re-runs,
   replicates and calibration.

---

## 6. Suggested order and verification gates

| Wave | Contents | Gate before proceeding |
|------|----------|------------------------|
| 0 | C1, C2, E5, E6, E9/E10, E11, minor | `Pkg.test()` green; new tests for α sign, equilibrium, quasi-extinction denominator, cache invalidation; no script writes a result whose digest is absent |
| 1 | C5+E1+E22, C7, E2/E3, E4 | model temperature equals corrected daily series at the 775 sites; run index records realised forcing; control has no trend; convergence metadata present |
| 2 | C4, E12 | graph has `n−1` edges, on-path edges, deterministic ties; obstacle accumulation documented |
| 3 | C3, E8, C6, E13–E19, E21, E7 | pre-registered criteria met; posterior predictive checks pass; team sign-off on parameter values |

**Definition of done for the corrected report:** every number in the report either traces to
a committed output with a matching configuration digest, or is explicitly labelled as a team
decision / unverified inference. The "brief" should not be circulated until Wave 1 has been
re-run on the corrected series.
