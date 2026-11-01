# GuadeX Phase 2 — Decision Record

**Date:** 2026-11-01.
**Basis:** §5 of `.kilo/plans/1793526105792-second-review-triage-plan.md`; answers recorded before any
Phase 2 re-run.
**Status:** decisions taken and Phase 2 executed at `c340398`. The primary post-fix results are the
`*_phase2` directories; the `*_corrected` directories remain the pre-fix audit baseline.

---

## 1. Decisions

| # | Decision | Chosen value | Where implemented / why deferred | Consequence |
|---|---|---|---|---|
| 1 | **Competition structure** | Keep **`c = 1 − α`** (shared per-site carrying capacity) as the primary specification; run **`c = −α`** (species-specific competition) as a sensitivity. | Implemented: `competition_structure = "shared_k"` is the default in `run_climate_scenarios.jl`; the sensitivity is selected with `GUADEX_COMPETITION_STRUCTURE='species_specific'`. The 45-run alternative is in `results/climate_scenarios_competition_alt` (evidence: `docs/GuadeX_Phase2_Competition_Sensitivity.md`). | The weak-warming magnitude is robust (basin native-richness delta ≤ 0.0129 species, ≤ 10 of 775 sites, under either structure). The **sign of the ensemble-mean richness response depends on the structure**: shared-K mean +0.0020 (26 positive / 4 zero / 14 negative), species-specific mean −0.0071 (44/44 negative); 30/44 runs flip sign. The shared-K richness-loss signal collapses to ≈0 under the alternative (mean −0.00155 → −0.000037; 31/44 exactly zero), while total-biomass and trout responses amplify (~4.8× and ~7.8–8.3×). Primary results are unchanged. |
| 2 | **Per-GCM controls** | Keep the **shared no-warming control** as an interim counterfactual and document it; **defer** the 11 matched per-GCM controls (P2.2). | Deferred by decision. The shared control is `control/baseline` (`run_climate_scenarios.jl:493`); `scenario_minus_control` uses `first(controls)` for every run (`src/climate_diagnostics.jl:235-255`). Documented in Limitations 2 and the response. | Deltas are scenario-minus-shared-control, not like-for-like by GCM. The control transient cancels, but a **residual GCM-specific offset remains**; it is stated in the report and limitations, not removed. |
| 3 | **Coded-value refit** | Keep the **national all-Spain air–water model** and document that choice; **defer** the coded-value refit and forcing regeneration (P2.5). | Deferred by decision. The projections use the Stage-B all-Spain model (`guadex_tw/scripts/10_project.py`; `13_project_guadex_sites.py`), with `ta = 0.5097`, `n = 39,944`, 980 sites in `projection_model_coefficients.json`. Documented in Limitations 3 and E20/E21 of the response. | The 257 coded Waterbase values do reach the GuadeX forcing. The forcing is conditional on the national calibration; **E20 remains a stated limitation**, together with the E21 same-day slope and lagged uncertainty band. |
| 4 | **Validation acceptance criteria** | Fix the thresholds that will define a successful **future calibration campaign** (they are **not** gates for this release): **per-species TSS ≥ 0.40 where ≥ 10 sites are occupied**; **\|modelled − observed occupancy\| ≤ 0.10 for LS, SA, SP, PW, CP**. | Recorded here as the future-calibration target; no code gate is added to this release. The Phase 2 validation run against the post-spin-up control is `results/climate_scenarios_phase2/validation/initial_state_validation.csv`. | The current control **passes the TSS criterion where it applies** (LS TSS 0.417) but **fails the occupancy criterion** for the widespread natives (LS observed 45.3 %, modelled 77.2 %, bias +31.9 points; similar over-prediction for SA, SP, PW and CP). These are genuine calibration targets, not release blockers. |

---

## 2. Phase 2 scope executed

All Phase 2 outputs were produced on 2026-11-01 at code version `c340398`. Exact commands, row
counts, elapsed times and digest evidence are in `docs/GuadeX_Phase2_Run_Record.md`.

| Item | Output root | Runs / rows |
|---|---|---|
| **P2.1a** climate ensemble re-run | `results/climate_scenarios_phase2` | 45 runs (44 GCM × SSP + 1 shared control) |
| **P2.1b** obstacle / upstream-cost sweep re-run | `results/sensitivity_obstacles_phase2` | 32 runs |
| **P2.1c** alternative-interaction sweep re-run | `results/alt_interactions_phase2` | 48 runs |
| **P2.3** competition-structure sensitivity | `results/climate_scenarios_competition_alt` | 45 runs |
| **P2.4a** passability floor `1e-3` | `results/sensitivity_obstacles_floor1e3` | 32 runs |
| **P2.4b** obstacles on parent edges only | `results/sensitivity_obstacles_parentonly` | 32 runs |
| **Validation** against the Phase 2 control | `results/climate_scenarios_phase2/validation/initial_state_validation.csv` | 24 species × 775 sites |
| **Test suite** (full Julia suite) | `results/phase2_test_suite.log` | 3,941 / 3,941 pass, 72 testsets, 0 failures |

**Pre-fix audit baseline (read only):** `results/climate_scenarios_corrected` (45 runs),
`results/sensitivity_obstacles_corrected` (32), `results/alt_interactions_corrected` (48), and
`data/BIOTIC/interaction_matrix_long_prefix.csv`.

**Deferred by decision:** P2.2 per-GCM controls (decision 2) and P2.5 coded-value refit /
forcing regeneration (decision 3). The report and its figures are regenerated from the Phase 2
outputs in the P2.6 report step (completed 2026-11-01); this decision record does not itself
change the `.tex` or `.pdf`.
