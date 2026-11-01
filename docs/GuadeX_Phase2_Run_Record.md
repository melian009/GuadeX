# GuadeX Phase 2 — Run Record

**Date:** 2026-11-01.
**Environment:** `E:\projects\GuadeX`, Windows, PowerShell 7, Julia 1.12.6.
**Code version:** `c340398` (`git rev-parse --short HEAD`). The Phase 1 parser fix and the Phase 2
edits are uncommitted working-tree changes.
**Scope:** every Phase 2 experiment, the validation run, and the full test suite. Exact commands,
output roots, row counts, elapsed time and digest evidence are recorded below. The pre-fix
`*_corrected` directories were read only and are retained as the audit baseline.

---

## 1. Run summary

| Item | Output root | Runs (index lines) | Elapsed | Identity evidence | Pre-fix baseline |
|---|---|---|---|---|---|
| **P2.1a** climate ensemble | `results/climate_scenarios_phase2` | **45** (46) | **≈23.1 min** (log 13:25:56 → 13:49:03; 1,386.7 s span) | `code_version = c340398` in 45/45; `run_digest.txt` = metadata `parameter_digest` in 45/45; example `ssp126/ACCESS-CM2` `f5f95873698951df` → `96a0c5621104ebb9` | `results/climate_scenarios_corrected` (45; `code_version = 9702169`) |
| **P2.1b** obstacle / upstream-cost sweep | `results/sensitivity_obstacles_phase2` | **32** (33) | **831.8 s** (13 min 52 s) | `parameter_digest` in 32/32, all changed; `.sensitivity_complete.jld2` in 32/32; **no `run_digest.txt`**; example `ssp126/IITM-ESM uc=0.01 baseline` `6b97f30399ade16c` → `580fba8364bccfe4` | `results/sensitivity_obstacles_corrected` (32; code version documented as `9702169`) |
| **P2.1c** alternative-interaction sweep | `results/alt_interactions_phase2` | **48** (49) | **1,234.3 s** script (20 min 34 s); 21 min 47 s wall | `parameter_digest` in 48/48 (24 distinct); `.sensitivity_complete.jld2` in 48/48; **no `run_digest.txt`**; example `original/sig_0.3/ssp126__IITM-ESM/uc_0.01/baseline` `314e1de67e715462` → `c071ab6d7ada9a9a` | `results/alt_interactions_corrected` (48; code version documented as `9702169`) |
| **P2.3** competition-structure sensitivity | `results/climate_scenarios_competition_alt` | **45** (46) | **30 min 08 s** (log 14:50:18 → 15:20:26; 1,808.1 s span) | `code_version = c340398` in 45/45; `resolved_config.competition_structure = species_specific` in 45/45; `run_digest.txt` = metadata digests; example `ssp126/ACCESS-CM2` `96a0c5621104ebb9` → `3c48bfabfa8b2ed5` | `results/climate_scenarios_phase2` (primary, same 45 keys) |
| **P2.4a** passability floor `1e-3` | `results/sensitivity_obstacles_floor1e3` | **32** (33) | **789.3 s** (13 min 09 s) | `passability_floor = 0.001` in 32/32; `parameter_digest` in 32/32, all changed from P2.1b; `.sensitivity_complete.jld2`; representative digest `f2d862ed4f4eeea2` | `results/sensitivity_obstacles_phase2` (32) |
| **P2.4b** parent edges only | `results/sensitivity_obstacles_parentonly` | **32** (33) | **806.2 s** (13 min 26 s) | `obstacles_on_path_only = true` in 32/32; `parameter_digest` in 32/32, all changed; `.sensitivity_complete.jld2`; representative digest `39b9faa914d04dc7` | `results/sensitivity_obstacles_phase2` (32) |
| **Validation** | `results/climate_scenarios_phase2/validation/initial_state_validation.csv` | **24 species × 775 sites** (25 lines incl. header) | not reported (completes in seconds; no timing persisted) | Reads `results/climate_scenarios_phase2/control/baseline`; 24 rows | `results/climate_scenarios_corrected/validation/initial_state_validation.csv` |
| **Test suite** (full Julia suite) | `results/phase2_test_suite.log`; final post-figure-fix re-run in `results/phase2_test_suite_final.log` (same 3,941/3,941 result) | **3,941 / 3,941 pass; 72 testsets; 0 failures/errors/broken** | 822.2 s log span (13:10:35 → 13:24:17) | 72 `Test Summary:` blocks; 28 more tests than the previous 3,913 | previous 3,913 / 3,913, 70 testsets |

All closing log lines report normal completion:
`Climate scenarios complete: 45 run(s)`, `Obstacle sensitivity complete: 32 run(s)`,
`Alternative-interaction sensitivity complete: 48 run(s)`; the logs contain no `ERROR`,
`LoadError`, `failed` or `Stacktrace`.

---

## 2. Exact commands

All commands were run from `E:\projects\GuadeX`. Logs are under `results/`.

### 2.1 P2.1a — climate ensemble (post-parser-fix primary)

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_phase2'
$env:GUADEX_CLIMATE_PLOT='0'
julia --project=. run_climate_scenarios.jl *> results/phase2_climate_run.log
```

### 2.2 P2.1b — obstacle / upstream-cost sweep (post-parser-fix primary)

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_phase2'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_run.log
```

### 2.3 P2.1c — alternative-interaction sweep (post-parser-fix primary)

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_ALT_OUTPUT_DIR='results/alt_interactions_phase2'
julia --project=. run_alt_interactions.jl *> results/phase2_alt_interactions_run.log
```

### 2.4 P2.3 — competition-structure sensitivity (`c = −α`)

The config-only smoke test that resolved the switch:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_COMPETITION_STRUCTURE='species_specific'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_competition_alt'
$env:GUADEX_CLIMATE_PLOT='0'
$env:GUADEX_CONFIG_ONLY='1'
julia --project=. run_climate_scenarios.jl
```

It printed `competition structure: :species_specific` and
`projection route: :interim_observed (interim_spin_up_years=3; control_run=true)`. The 45-run
ensemble:

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_COMPETITION_STRUCTURE='species_specific'
$env:GUADEX_CLIMATE_OUTPUT_DIR='results/climate_scenarios_competition_alt'
$env:GUADEX_CLIMATE_PLOT='0'
julia --project=. run_climate_scenarios.jl *> results/phase2_competition_alt_run.log
```

### 2.5 P2.4a — compounded-passability floor `1e-3`

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_OBSTACLE_FLOOR='0.001'
$env:GUADEX_OBSTACLES_PARENT_ONLY='0'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_floor1e3'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_floor1e3_run.log
```

### 2.6 P2.4b — obstacles on parent edges only

```powershell
$env:GUADEX_PARAMETERS_FILE='parameters_climate_scenarios_corrected.toml'
$env:GUADEX_OBSTACLE_FLOOR='0'
$env:GUADEX_OBSTACLES_PARENT_ONLY='1'
$env:GUADEX_SENSITIVITY_OUTPUT_DIR='results/sensitivity_obstacles_parentonly'
julia --project=. run_sensitivity_report.jl *> results/phase2_obstacle_parentonly_run.log
```

### 2.7 Validation against the Phase 2 control

```powershell
$env:GUADEX_VALIDATION_RUN='results/climate_scenarios_phase2/control/baseline'
$env:GUADEX_VALIDATION_OUT='results/climate_scenarios_phase2/validation'
julia --project=. scripts/validate_initial_state.jl
```

The script writes `initial_state_validation.csv` under the output directory. Result: pooled
sensitivity **0.975**, pooled specificity **0.920**, LS TSS **0.417** (observed occupancy 45.3 %,
modelled 77.2 %). The values are identical at reported precision to the pre-fix validation.

### 2.8 Full test suite

```powershell
julia --project=. test/runtests.jl *> results/phase2_test_suite.log
```

The log records **72** `Test Summary:` blocks and **3,941 / 3,941** passing tests, with **0**
matches for `Fail`, `Error` or `Broken`. This is **28** more assertions than the previous
3,913 / 3,913, 70-testset run. The command itself is not echoed in the log; the invocation above is
the reproducible equivalent against `test/runtests.jl`.

**Clean-checkout caveat:** on a checkout without the gitignored daily forcing file, the C5/E1
testset *model temperature equals corrected daily series* skips with one `@test_skip` (shown as
*Broken*) instead of erroring. The Phase 2 log was produced with the forcing file present; no exact
clean-checkout denominator was recorded for this run.

---

## 3. Digest evidence and identity convention

- **Climate runs** (`run_climate_scenarios.jl`, both the primary and the competition alternative)
  write `run_digest.txt` **and** `export/run_metadata.json`; the two digests are identical in every
  run. The metadata also records `code_version = c340398` and
  `resolved_config.competition_structure`.
- **Sweeps** (`run_sensitivity_report.jl`, `run_alt_interactions.jl`) write the
  `parameter_digest` in `export/run_metadata.json` **and** a `.sensitivity_complete.jld2`
  completion marker/fingerprint. They do **not** write `run_digest.txt`, and their metadata
  carries no raw `code_version` key. The code version still enters every digest via
  `Guadex.code_version() = c340398` (`src/run_identity.jl`), and all post-fix digests differ from
  the pre-fix sweeps.
- The interaction matrix is a digest input. All post-fix digests differ from the pre-fix runs, so
  the new directories are provably not reuses of the audit baseline.
- The realised forcing is bit-identical between the pre-fix and post-fix climate, obstacle and
  interaction runs (0 differing rows in the five forcing columns checked); the comparison
  therefore isolates the interaction-matrix change.

---

## 4. Pre-fix audit baseline

| Baseline | Rows | Code version evidence | Digest evidence |
|---|---|---|---|
| `results/climate_scenarios_corrected` | 45 (46) | `code_version = "9702169"` in metadata; example `ssp126/ACCESS-CM2` | `parameter_digest = f5f95873698951df`; `run_digest.txt` present |
| `results/sensitivity_obstacles_corrected` | 32 (33) | no `code_version` key; documented as `9702169` | example `parameter_digest = 6b97f30399ade16c`; `.sensitivity_complete.jld2`, no `run_digest.txt` |
| `results/alt_interactions_corrected` | 48 (49) | no `code_version` key; documented as `9702169` | 24 distinct `parameter_digest`; example `314e1de67e715462`; `.sensitivity_complete.jld2`, no `run_digest.txt` |
| `data/BIOTIC/interaction_matrix_long_prefix.csv` | 227 directed rows | pre-fix parser output | retained byte-for-byte (SHA-256 `2716DC41…C55098`) |

The `*_corrected` directories were not overwritten and retain their original timestamps (run
directories 30 Sep 2026; indices and report plots 29 Oct 2026).

---

## 5. Derived Phase 2 evidence files

These are build artefacts under the gitignored `results/phase2/` tree; they are the source of the
numbers in the Phase 2 evidence documents.

| File | Rows (incl. header where CSV) | Content |
|---|---|---|
| `climate_prefix_vs_phase2.csv` | 271 | 45 runs × 6 metrics, pre/post values and deltas |
| `parser_fix_effect_summary.csv` | 7 | per-metric parser-fix effect summary |
| `basin_aggregate_prefix_vs_phase2.csv` | 136 | 45 runs × 3 final-year basin aggregates |
| `interaction_changed_cells.csv` | 35 | 34 changed directed cells (16 reverse + 18 ambiguity flags) |
| `obstacle_prefix_vs_phase2.csv` | 257 | 32 runs × 8 metrics, pre/post and deltas |
| `obstacle_effect_summary.csv` | 9 | per-metric obstacle-sweep delta summary |
| `obstacle_rank_comparison.csv` | 25 | factor ranking pre vs post (8 metrics) |
| `obstacle_factor_contrasts.csv` | 41 | 40 signed factor contrasts with sign-flip flag |
| `obstacle_digest_comparison.csv` | 33 | 32 run-directory digest pairs |
| `alt_interactions_prefix_vs_phase2.csv` | 385 | 48 runs × 8 metrics, pre/post and deltas |
| `alt_interactions_effect_summary.csv` | 9 | per-metric interaction-sweep delta summary |
| `alt_interactions_matrix_order_comparison.csv` | 9 | matrix ordering pre vs post (8 metrics) |
| `alt_interactions_matrix_block_order.csv` | 9 | per-block ordering checks |
| `alt_interactions_digest_comparison.csv` | 49 | 48 run-directory digest pairs |
| `competition_structure_deltas.csv` | 89 | 44 scenarios × 2 structures, deltas and sign-flip flags |
| `competition_endpoint_summary.csv` | 25 | 6 metrics × 2 endpoints × 2 structures |
| `competition_sign_flips.csv` | 103 | every sign-flip run per metric |
| `obstacle_floor_vs_baseline.csv` | 257 | 32 runs × 8 metrics, floor variant − baseline |
| `obstacle_parentonly_vs_baseline.csv` | 257 | 32 runs × 8 metrics, parent-only − baseline |
| `passability_variant_summary.csv` | 17 | per-variant per-metric summary and ranking |
| `passability_link_diagnostics.csv` | 342 | per-link obstacle counts and compounded factors |
| `passability_matching_diagnostics.csv` | 21 | matched/applied/root-join change counts |
| `passability_floor_scenario_counts.csv` | 5 | floor burden and 33-link share by scenario |
| `passability_variant_by_scenario.csv` | 65 | paired deltas by passability scenario |
| `passability_rank_detail.csv` | 73 | factor scores for baseline and both variants |
| `passability_variant_verification.csv` | 7 | index/digest/metadata checks |

---

## 6. Release-checklist status (plan §6)

| Gate | Status | Evidence |
|---|---|---|
| **Evidence** | **Met** | `docs/GuadeX_SecondReview_Assessment_Oct2026.md` (objection-by-objection dossier) and the five Phase 2 evidence documents (`ParserFix_Effect`, `ObstacleSweep_Update`, `InteractionSweep_Update`, `Competition_Sensitivity`, `PassabilityFloor_Sensitivity`). |
| **Documentation** | **Met** | The decisions, limitations, status record and response are updated, and `docs/FinalReportOctober2026.tex` / `.pdf` has been regenerated from the Phase 2 outputs in P2.6 (2026-11-01; 39 pages). |
| **Code** | **Met** | Full Julia suite **3,941 / 3,941**, **72 testsets**, **0 failures/errors/broken** (`results/phase2_test_suite.log`; the final re-run after the figure-regressor fix, `results/phase2_test_suite_final.log`, gives the same 3,941/3,941); parser regression tests pass; the two-species equilibrium/stability test
pins the expected stability. |
| **Data** | **Met** | Regenerated long table is **227 directed rows**, **24×24** matrix, **16** reverse-direction cells corrected; **31 ambiguous source cells / 62 rows** (22/44 before); pre-fix table retained as `data/BIOTIC/interaction_matrix_long_prefix.csv`. |
| **Validation** | **Met** | Per-species table at `results/climate_scenarios_phase2/validation/initial_state_validation.csv` (24 species × 775 sites); pooled sensitivity 0.975, specificity 0.920, LS TSS 0.417; control change per species reported. |
| **Phase 2 gate** | **Met** | Decisions recorded in `docs/GuadeX_Phase2_Decisions.md` on 2026-11-01 before the runs. P2.1, P2.3, P2.4 and P2.6 executed; P2.2 and P2.5 deferred by decision. |

---

## 7. What the tree does not record

- The test-suite log does not echo the command; the invocation in §2.8 is the reproducible
  equivalent, not a captured shell line.
- The validation run's wall time is not persisted.
- The sweep metadata does not contain a raw `code_version` key; identity is the
  `parameter_digest` (which folds `c340398`) plus the completion-marker fingerprint. This is a
  writer limitation, not a missing run.
- `results/` is gitignored; the new run directories are present on this machine and are not part
  of the commit. The shared archive `https://k00.fr/eqnsvc2y` is the external reproduction route.
