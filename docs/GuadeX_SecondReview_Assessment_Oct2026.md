# GuadeX — Second Review: Evidence Dossier and Dispositions

**Status date:** 2026-11-01.
**Review under analysis:** `docs/GuadeX_Reply_to_Response_Oct2026.md` (Carlos, October 2026).
**Snapshot checked:** `master` at `c34039809cf1e402d7deff78c50063cbf58513b9` (base `9702169`,
correction-series tip `61fcd4c`, AAL reclassification `099a87d`, 29 Oct 2026).
**Scope:** Phase 1 only — evidence, documentation, cheap code/test fixes and post-processing.
No experiment was re-integrated; the Phase 2 re-runs stay behind the team gate (§6).

This dossier records, for every objection in the second review, what the current tree shows, the
value we recomputed, the verdict, and the action taken. Every number below has a file path or a
script that reproduces it. The companion deliverables are the revised response
(`docs/GuadeX_Response_to_Reviewers_Oct2026.md`), the report
(`docs/FinalReportOctober2026.tex`), the limitations note
(`docs/GuadeX_Limitations_and_Scope_Oct2026.md`) and the status record
(`docs/GuadeX_Correction_Status_Sept2026.md`).

---

## 1. Frozen snapshot and hashes

| Artefact | Value |
|---|---|
| `git rev-parse HEAD` | `c34039809cf1e402d7deff78c50063cbf58513b9` |
| Review document SHA-256 | `06D6DEEF536A1E38912B97D9D4E7D01FCD02A101BC1C735D66457C553A9BC667` |
| `docs/FinalReportOctober2026.pdf` at `c340398` (blob) | `5280AC915292C6E9A423E557D1D0F9EB006469E0688F560C70D6D66ACF4247C7` |
| Regenerated Phase 1 PDF (36 pages) | `6D4B0418F898CC1DD539F35FD3BFD5D4140EA9AE52C596A27B0B6C34FA2D265B` |
| Pre-fix `interaction_matrix_long.csv` SHA-256 | `2716DC41FEAD2464F712869B6264D5B0E2428598699E4DB790FCFAFBA0C55098` |
| Post-fix `interaction_matrix_long.csv` SHA-256 | `31A83BB16E1CFEBB9191D0497607B3CDC01CF37D59EFA4F0340B0F120418AA8E` |
| Pre-fix copy | `data/BIOTIC/interaction_matrix_long_prefix.csv` (unchanged, same hash) |

**Which snapshot each objection was checked against.** Objections 4.1–4.7, 4.9 and 4.10 were
checked against the current tree. Objection 4.8 was written against a pre-`099a87d` snapshot:
`git show 099a87d^:parameters_climate_scenarios_corrected.toml` line 132 and
`git show 099a87d^:parameters.toml` line 85 both hold
`migratory = ["AA", "AAL", "LR", "MC"]`, while in the current files those lines are the
`[species]` headers. The reviewer's line numbers therefore identify the old file exactly.

---

## 2. Dispositions

### 4.1 The control is shared, not matched by GCM — **Accept (text)**

| Evidence | Value |
|---|---|
| Control GCM assignment | `run_climate_scenarios.jl:493` — `scenario_gcms = is_control ? ["baseline"] : CLIMATE_GCMS` |
| Delta function | `src/climate_diagnostics.jl:235-255` — `control = first(controls)` for every scenario |
| Control forcing | `src/temperature_forcing.jl:251-273` — day-of-year climatology of the ensemble-median historical series |
| Run index | `results/climate_scenarios_corrected/runs_index.csv` — 44 scenarios + 1 control |

The review is right: there is one control, and it is not driven by each GCM's own daily series.
The control's anomalies are zero by construction (the `control` row of `runs_index.csv` reports
`-2.18e-17 °C`), so "matched" language cannot be recovered by pointing at the GCM pairing.

**Action:** the response, report and limitations now say "shared control" and state the residual
GCM-specific offset. Per-GCM controls (11 runs) are a Phase 2 option.

### 4.2 The quoted slope drives the projections — **Accept (text); withdraw the rebuttal**

| Evidence | Value |
|---|---|
| Stage-B model | `guadex_tw/scripts/10_project.py:4,46-55` — all-Spain fit |
| Loader | `guadex_tw/scripts/13_project_guadex_sites.py:370-374` loads `10_project.py` dynamically |
| Stored coefficients | `guadex_tw/models/projection_model_coefficients.json` — `ta = 0.5097248`, `ta_z = −0.1629927`, `n = 39944`, `n_sites = 980` |

The national (all-Spain) model produces the GuadeX water-temperature forcing, so the quoted slope
is the one that matters and the 257 coded Waterbase values do reach the projections. Two claims
must be separated: (a) the residual-SD improvement 2.65 → 2.34 °C belongs to the Spain model and
the Guadalquivir-only SD is 2.54 — this narrows *which model's SD improves*; (b) the slope that
drives the projections is also the Spain model's — this part of the old rebuttal is invalid and is
withdrawn.

**Action:** E20/E21 and §3 items 11–12 rewritten. Keeping the national model is now an explicit
choice; a coded-value refit that regenerates the forcing is a Phase 2 decision.

### 4.3 The 21/494/532 burn-in years do not belong to the corrected runs — **Accept (text)**

| Evidence | Value |
|---|---|
| Route selection | `parameters_climate_scenarios_corrected.toml:259-260` (obstacle sweep) and `:374-378` (climate runs) — `projection_route = "interim_observed"`, `interim_spin_up_years = 3` |
| Rule bypass | `sensitivity_core.jl:481-494` — interim branch integrates exactly 3 years |
| Corrected summary | `results/alt_interactions_corrected/burnin_summary.csv` — 3 rows, `years=3`, `criteria=basin`, `converged=false` for `original`, `random`, `invasive_favoring` |

The E4 rule is implemented and tested but inactive in every corrected run. The 21/494/532-year
figures are September full-burn-in values; the first response wrongly attached them to the
corrected file.

**Action:** E4 and §3 item 7 corrected; report and status record now state 3 years / interim route.

### 4.4 Residual parsing error in the interaction table — **Accept (bug fix + gated re-run)**

The undirected clause was cut as `lower[m.offset:end]`, so a cell of the form
`interfere through competition and  Gh through predation` gave the reverse direction the
predation value. The parser now cuts the undirected clause at the next named clause
(`src/data_preparation.jl:906-932`).

| Item | Before | After |
|---|---|---|
| Reverse-direction rows changed −0.5 → −0.3 | — | **16** (table below) |
| Ambiguous source cells / rows | 22 / 44 | **31 / 62** |
| Pure "competition and predation" cells flagged | 0 of 9 | **9 of 9** (18 rows) |
| Long-table rows | 227 | 227 |
| 24×24 matrix | 24×24 | 24×24 |

The 16 corrected rows (all −0.5 predation → −0.3 competition):

`GH→AAL`, `GH→CC`, `GH→AM`, `GH→LG`, `AAL→CC`, `AAL→AM`, `AAL→OM`,
`CG→EL`, `CG→AAL`, `CG→CC`, `CG→AM`, `CG→LG`, `TT→CC`, `TT→AM`, `TT→LG`, `LG→OM`.

The nine pure mixed cells — MS–EL, EL–CC, EL–AM, EL–OM, EL–LG, CC–AM, CC–LG, AM–LG, LG–AAL —
keep −0.5 in both directions as a documented judgement and are now flagged `ambiguous = true`.
Directed multi-mechanism text such as "competition and sperm parasitism" is preserved.

**Action:** parser fixed; long table regenerated; pre-fix copy retained
(`interaction_matrix_long_prefix.csv`); regression testset added (105 tests) covering the 16
rows, the named direction, the nine cells, competition-only cells and the sperm-parasitism text.
The affected climate/obstacle/interaction experiments were **not** re-run (Phase 2, §6).

### 4.5 Competition coefficients have strong ecological consequences — **Accept in substance; keep the specification**

| Evidence | Value |
|---|---|
| Interaction block | `src/ode_model.jl:168-181` — `c_sj = 1 − α_sj`, `c_ss = 1` |
| Observed α | −0.3, −0.5, −0.8 → `c ∈ {1.3, 1.5, 1.8}` |
| Consequence | `c₁₂c₂₁ > 1` for every interacting pair; the interior equilibrium `K(1−c₁₂)/(1−c₁₂c₂₁)` is a saddle |
| Test coverage | `test/test_ode.jl` C2 testset, **35 tests** (13 before) |

The reviewer is right on all counts: the two-species equilibrium test claimed in the first
response did **not** exist; the chosen mapping puts every pair in the exclusion/priority-effect
regime; and `c = −α` is a valid alternative LV structure, not a mathematical error.

**Action:** the false test claim is withdrawn and the missing test added. The new tests verify the
coexistence formula for `c₁₂c₂₁ < 1`, the positive Jacobian eigenvalue for the observed regime,
bistability under two initial conditions, and the `c = 1 − α` range. The report and limitations
state the ecological consequence. `c = −α` / species-specific competition is a Phase 2
structural sensitivity.

### 4.6 The starting state is not validated against the survey — **Accept (validation added)**

New script: `scripts/validate_initial_state.jl`. It loads the control run's saved post-spin-up
state (`baseline_density` in `export/levels/species_timeseries.csv`), so no re-integration is
needed, and compares it with the 2006–2009 survey over 775 sites and 24 species.

| Quantity | Value |
|---|---|
| Observed presence | survey density > 0 |
| Modelled presence | spun-up density > 0.1 |
| Pooled sensitivity | 0.975 |
| Pooled specificity | 0.920 |
| Widespread native LS | TSS 0.42, observed occupancy 45.3 %, modelled 77.2 % |
| Rare species | TSS 0.49–1.00, driven by correctly empty sites (the low values are AA 0.55 and TT 0.49) |
| Abundance rank ρ | 0.48–1.00 where ≥3 sites were occupied |

Raw transient in the same control (2026 → 2045): total biomass −35.5 %, native biomass −40.8 %.
The validation is a finding, not a success: the model over-predicts occupancy, especially for the
dominant native cyprinid, and several invasive species move strongly in the control (for example
GH +241 %, MS −87.9 %, LG −84.5 %).

Full per-species table:
`results/climate_scenarios_corrected/validation/initial_state_validation.csv` (24 rows); the
report reproduces it as Table~\ref{tab:validation} and the response quotes the headline.

**Action:** validation table added to the report and the evidence dossier; the TSS/occupancy
criteria remain the acceptance condition for calibration (§6).

### 4.7 Parts of the network are reconstructed rather than observed — **Accept (reporting + sensitivity)**

New script: `scripts/network_root_join_diagnostics.jl`. Recomputed on the corrected config:

| Quantity | Value |
|---|---|
| Parent edges | 508 (8,636.9 km, 62.3 %) |
| Root-join edges | 266 (5,234.5 km, **37.7 %**) |
| Roots | 267 of 775 sites (34.5 %) |
| Root-join length | median 15.8 km, max 206.6 km |
| Elevation inversions | 4 parent, **41 root** |
| Obstacles on root-join edges | **341 matched / 245 applied** |
| Obstacles incident to a root site | 531 matched / 392 applied |
| Links with obstacles | 330 of 774 |
| Max obstacles on one link | **33** (upstream factor `1.0e-33`) |
| Links with ≥5 obstacles | **28** |
| Links with factor < 1e-6 | 16 |
| Roots at (0.01, 100 m) / (0.02, 200 m) / (0.05, 500 m) | **440 / 267 / 119** |

The review's 38 %, 33 obstacles / 1e-33, 28 links, 206.6 km and **341 matched / 245 applied**
root-join obstacle count all reproduce exactly. An earlier version of the diagnostic reported
297 / 215 because it canonicalised root pairs by site index rather than by site code; 45 root
joins have different index and string ordering, which dropped 44 mapping rows. The script now
canonicalises by site code (`scripts/network_root_join_diagnostics.jl` line 67) and agrees with
the review. A broader definition that also counts parent edges entering a root site gives
531 matched / 392 applied; that is a different quantity and is reported only as a cross-check.

**Action:** figures added to report §2.2 and limitation 8; a passability floor and on-path-only
obstacle assignment are Phase 2 structural sensitivities.

### 4.8 *Alburnus alburnus* still classified as migratory — **Already fixed at `099a87d`**

| Evidence | Value |
|---|---|
| Fix commit | `099a87d` (29 Oct 2026), "Refactor species classification and metrics in response to reviewer feedback" |
| Current `parameters.toml` | `:86-88` — invasive list includes `AAL`; migratory `["AA","LR","MC"]` |
| Current corrected config | `:133-135` — same groups |
| Tests | `test/test_parameters.jl` species-group assertions |
| Pre-fix evidence | `git show 099a87d^` puts AAL in migratory at exactly lines 85 and 132 |

AAL-dependent metrics were recomputed from saved states without re-integrating
(`results/climate_scenarios_corrected/invasive_reclassified_AAL.csv`, 45 rows: 44 scenarios plus the control, which the means below include):

| Metric | Before reclassification | After |
|---|---|---|
| Invasive richness (mean) | 0.7130 | **0.7488** (0.7432–0.7523) |
| Invasive biomass (mean) | 4.892 | **4.936** |
| Change | — | +0.036 species/site, +0.044 units |

The response's "≈0.71" is the pre-reclassification value; the report's 0.749 is current. Species
counts are 10 native, 11 invasive, 3 migratory.

**Action:** response E17 and the §2.4 rows corrected; report already used the current values.

### 4.9 Results cannot be reproduced from the repository — **Resolved by a shared archive**

The corrected outputs and daily forcing were shared with the reviewer at
`https://k00.fr/eqnsvc2y`. The repository still excludes `results/` (`.gitignore:32`) and the
daily wide forcing (`guadex_tw/.gitignore:15`), so a clone cannot rebuild the numbers; the shared
archive is the reproduction route for this release. The clean-checkout failure is also fixed
cheaply: the testset *model temperature equals corrected daily series (C5/E1)* now skips with an
explicit message when the forcing file is absent instead of erroring.

**Action:** link repeated in the response, report and status record; report availability sentence
now distinguishes the demonstration viewer from the corrected archive.

### 4.10 Smaller inconsistencies — **Accept (documentation)**

| Claim | Current value / action |
|---|---|
| Limitations still say 18.37 °C and ~1 °C | Fixed to **18.06 °C** and **~0.6 °C** (`GuadeX_Limitations_and_Scope_Oct2026.md`, item 1) |
| Report headlines span-based mean ranks but response says ranking replaced | Balanced within-design sums of squares are now the headline in the abstract and the hydraulic-controls subsection (Results 4.4); the mean rank 1.88 is retained as a descriptive cross-metric summary |
| E18 says 1,584 cells | Measured **528** AB cells with the corrected biological options; `prepare_ode_data(...; salinity_envelope = true, nonreproducing_local_growth = :zero)` returns `salinity_envelope_filtered_cells = 528` once E16 has zeroed LR/MC, while the legacy E16 default returns 1,584 |
| Endpoints quoted without the window | Endpoints are selected on **2036–2045**: 0.180 and 1.292 °C; the 2026–2045 means of the same runs are 0.208 and 1.091 °C. Both windows are now named |
| Alien introductions absent from report limitations | Added as report limitation 13 and limitations-note item 12 |
| Response/status say work is uncommitted | Fixed and made precise: the correction work through `099a87d` is committed on `master` (series tip `61fcd4c`), and the second-review edits are explicitly labelled uncommitted working-tree changes; HEAD `c340398` |

---

## 3. Section-2 corrections (the reviewer's own corrections) — acknowledged

The review's §2 corrects its own earlier claims; we accept all of them as stated:
C1 (79 pairs, Σα = −58.1 without the *Lepomis* row, as the pre-correction loader computed it;
the current 24-species long table gives 82 invasive→native rows summing to −56.6 with the
report's 10-native / 11-invasive groups), C2 (`c = 1 − α` is what the code does),
C4 (8,640 km is the parent-edge sum; the tree is 13,871 km), C7 (r ≈ 0.77, not 0.81),
E10 (settings were documented, the gap was a committed config), E20 (`resid` is already an
absolute value) and the minor-number corrections. Two differences are reference-set effects, not
errors: the 294-vs-293 fishless sites (776 environmental rows vs 775 model sites) and the
21-year burn-in being correct for one of the September matrices.

## 4. Section-3 confirmations — accepted

The review independently reproduced the corrected test suite, the interaction direction and
24×24 matrix, the inside-growth LV form, the interim route, the network totals, the corrected
temperature range (9.79–18.06 °C; basin mean 15.82 °C; trout sites 12.83 °C), offsets 1–20,
established-only quasi-extinction, the obstacle status counts (1,274 / 185 / 199),
E13–E16/E22 and the viewer items A1–A8. Nothing in this dossier re-litigates those.

## 5. Phase 1 validation outputs

| Gate | result |
|---|---|
| Parser regression tests | 105 tests pass; 16 reverse rows corrected; 9 cells flagged |
| C2 stability tests | 35-test testset passes in both regimes |
| Full Julia suite | **3,913 / 3,913 pass, 70 testsets, 0 failures, 0 errors** (sum of the per-testset totals, 2026-11-01; an outer-`@testset` run of the same files reports the same 3,913). Clean checkout without the daily forcing file: 3,898 total, 3,897 pass + 1 skipped (`@test_skip`) |
| Post-spin-up validation | 24 species × 775 sites; pooled sensitivity 0.975, specificity 0.920 |
| Network diagnostics | 508 parent / 266 root joins; 37.7 % length; 341 matched / 245 applied root-join obstacles (531/392 incident-to-root); 33 obstacles / 1e-33; 28 links ≥5 |
| Document consistency | the 4.10 items are closed in the response, report, limitations and status record; `docs/climate_scenarios.md:566` now gives both the 1,584 legacy-default and the 528 corrected counts, and `docs/global_sensitivity_methods_note.md:65` uses 18.06 °C. The old values survive only as explanatory contrasts |

## 6. Decisions and the Phase 2 gate

Still required from the team (no Phase 2 run has started):

1. **Competition structure.** Keep `c = 1 − α` as primary and run `c = −α` as a sensitivity
   (recommended), or switch the primary specification.
2. **Per-GCM controls.** Keep the shared control for this release (recommended) or run the 11
   matched controls now.
3. **Coded-value refit.** Refit the all-Spain model and regenerate the forcing, or document the
   national-model choice and defer.
4. **Validation acceptance criteria.** Fix the TSS/occupancy thresholds that will define a
   successful calibration. Proposed starting point: per-species TSS ≥ 0.4 for every species with
   ≥10 occupied sites, and |modelled − observed occupancy| ≤ 0.10 for the widespread natives
   (LS, SA, SP, PW, CP); the current LS TSS of 0.42 and occupancy bias of +31.9 points fail the
   second condition and would therefore be a genuine calibration target.

Phase 2 items are gated in `.kilo/plans/1793526105792-second-review-triage-plan.md` §4:
P2.1 parser-affected re-runs, P2.2 per-GCM controls, P2.3 competition-structure sensitivity,
P2.4 passability floor, P2.5 coded-value refit (conditional), P2.6 report regeneration.

## 7. Residual risks and honest limits

- **The corrected run numbers still reflect the pre-fix interaction table.** The parser fix
  changes 16 α entries; the climate, obstacle and interaction results were not re-run. Every
  number in the report is therefore a pre-4.4 result until Phase 2. The saved runs record
  `code_version = "9702169"`, so a current re-run will (correctly) refuse to reuse them under the
  digest rule.
- **The shared control carries a GCM-specific offset.** It is now stated, not removed.
- **The start-state validation shows over-prediction.** TSS is flattered by rare species; the LS
  result is the honest one.
- **The first network diagnostic undercounted root-join obstacles.** An index-vs-code ordering
  bug gave 297/215; the fixed script gives the review's 341/245. The broader incident-to-root
  definition is 531/392. This is a fixed defect, but any earlier note that quoted 297/215 is
  superseded.
- **Three parser gaps are latent and currently untriggered.** The directed-verb branch evaluates
  the whole cell rather than the clause (no cell contains two directed verbs); `dispersa` is not
  in the verb list (the word does not occur in the data); and when a bare-target and a `depreda`
  clause name the same ordered pair, the first specific clause wins (no current cell has that
  shape). They are hardening items for Phase 2, not defects in the committed 24×24 matrix.
- **`results/` is not in the repository.** The shared archive is the only reproduction route until
  a persistent accession exists.
- **The `c > 1` regime may drive part of the weak-warming conclusion.** The Phase 2 sensitivity
  must be reported, not buried.
