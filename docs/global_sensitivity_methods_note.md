# Global sensitivity (E7): methods note

Status: implementation + estimator validation + **pilot** run. The pilot is *not*
the production index. This note is written so the October report can cite it and
so the team can confirm the parameter ranges before the production index is
quoted.

This addresses correction **E7** ("replace the cross-design factor ranking with
within-design sums of squares **plus** a global design (Latin hypercube + Sobol
indices) over agreed parameter ranges"). The within-design ranking is **kept
unchanged**; the global design is added alongside it.

## 1. Which analysis answers which question

| Analysis | Files | Question it answers | Valid scope |
|---|---|---|---|
| Within-design factor ranking | `scripts/plot_sensitivity_effects.jl` → `parameter_effect_rank.csv`, `variance_decomposition.csv` | Which swept factor moves the response most **at the chosen factor levels**, inside one balanced full-factorial sweep? | Only inside that balanced sweep; the thermal-sigma axis is confounded across sweeps |
| Global Sobol design (new) | `src/global_sensitivity.jl`, `scripts/global_sensitivity_sobol.jl` → `sobol_indices.csv` | How much of the response **variance** is attributable to each parameter and its **interactions** over its **range**? | The declared parameter box `[global_sensitivity]`; assumes independent uniform inputs |

Both are legitimate and complementary; neither replaces the other.

## 2. Method

Implemented in `src/global_sensitivity.jl` in pure Julia (no new dependency; the
SplitMix64 PRNG and the Saltelli design are written in-repo).

1. **Latin hypercube** (`latin_hypercube`, `lhs_samples`): each of the `D`
   columns is a random permutation of the `n` equiprobable strata with one
   uniform draw per stratum, giving marginally uniform samples with exact
   one-per-stratum coverage in every dimension. A fixed `seed` makes the design
   deterministic.
2. **Saltelli A/B/AB/BA design** (`saltelli_design`): `A` and `B` are independent
   uniform samples; for each parameter `i`, `AB[i]` is `A` with column `i` taken
   from `B`, and `BA[i]` the reverse. Total model evaluations
   **`n * (2D + 2)`**.
3. **Estimators** (`sobol_indices_from_outputs`), averaged over the AB and BA
   orientations (Saltelli et al. 2010; Jansen 1999), with `V_k` the variance of
   the pooled A∪B outputs for metric `k`:

   ```
   S1_i = mean( Y_B[:,k] .* (Y_AB[i][:,k] .- Y_A[:,k]) ) / V_k
   ST_i = mean( (Y_A[:,k] .- Y_AB[i][:,k]).^2 ) / (2 V_k)
   ```

4. **Bootstrap confidence intervals**: `n` rows are resampled **with
   replacement** `bootstrap` times and the percentile interval is reported. (A
   permutation resample would preserve every sample mean and collapse the
   intervals; there is a regression test for this.)
5. **Driver** (`sobol_analyze`): maps every design row onto the parameter box
   and an `evaluate(x)` closure returning the response metrics.

## 3. Parameter set and agreed-range default proposal (team to confirm)

Defined in `parameters.toml` under `[global_sensitivity]`. The defaults are the
**proposal** assembled from the correction plan and the data's own uncertainty;
they are an interpretation of uncertainty ranges, **not calibrated posteriors**.
Each maps to an existing model parameter through an existing `set_*` helper (no
model equation changes, no run-offset changes):

| Parameter | Range | Default | Transform | Rationale |
|---|---|---|---|---|
| `growth_multiplier` | 0.75–1.25 | 1.0 | `scale_intrinsic_growth_rate` | Literature local `r` span source populations by roughly ±25% |
| `carrying_capacity_scaling` | 0.5–2.0 | 1.0 | `scale_carrying_capacity` | The E14 isolated-pool distortion (84 pools, ~56% of basin density, max K ≈ 95,333) makes the capacity level the least certain state variable |
| `interaction_strength_multiplier` | 0.5–1.5 | 1.0 | `set_interaction_matrix` | The interaction matrix is qualitative (displaces −0.8, predation −0.5, competition −0.3); the factor propagates that coding uncertainty without changing signs |
| `thermal_optimum_shift_c` | −2.0–2.0 | 0.0 | `set_thermal_optima` | Trait-table optima are range midpoints; C6 notes all baseline sites are ≤ 18.37 °C |
| `thermal_breadth_multiplier` | 0.6–1.4 | 1.0 | `set_thermal_sigma_multiplier` | `sigma ≈ range/6` is an approximation (C6) |
| `dispersal_multiplier` | 0.5–2.0 | 1.0 | `scale_dispersal` | Species dispersal rates carry order-of-magnitude literature ranges (E17); the base matrix uses the median species rate |
| `upstream_cost` | 0.01–0.5 | 0.05 | `precompute_dispersal_matrix` | Exactly the existing obstacle-sweep grid `{0.01, 0.05, 0.1, 0.5}`, treated continuously |
| `heat_stress_k_multiplier` | 0.0–2.0 | 1.0 | `set_heat_stress_rate` | `k` is calibrated to a `max_annual_loss = 0.05` assumption (WP3); 0 disables the term |

The base heat-stress slope used by the pilot is the corrected obstacle-sensitivity
calibrated value `k = 3.81183578844573e-05` (recorded in
`results/sensitivity_obstacles_corrected/.../run_metadata.json`).

**Decision needed:** the team should confirm these ranges (and the `n`,
warming forcing and burn-in) before the production index is cited.

## 4. Estimator validation against the analytic Ishigami function

The Ishigami function `f = sin(x1) + 7 sin²(x2) + 0.1 x3⁴ sin(x1)` on
`[-π, π]³` has known indices. With `n = 2¹³` base samples (`n(2D+2) = 65,536`
evaluations) and 300 bootstrap replicates
(`test/test_global_sensitivity.jl`):

| Parameter | S1 analytic | S1 estimate [95% CI] | ST analytic | ST estimate [95% CI] |
|---|---|---|---|---|
| x1 | 0.3139 | **0.3104** [0.2889, 0.3301] | 0.5576 | **0.5395** [0.5146, 0.5638] |
| x2 | 0.4424 | **0.4492** [0.4326, 0.4650] | 0.4424 | **0.4467** [0.4359, 0.4576] |
| x3 | 0.0000 | **0.0129** [−0.0018, 0.0285] | 0.2437 | **0.2396** [0.2310, 0.2485] |

Analytic variance 13.8446 vs estimated 13.5994 (1.8% low). All estimates are
within 0.02 of the analytic values and the intervals bracket them; x3 is correctly
identified as having zero first-order and non-zero total-order effect. The unit
tests also validate LHS marginal stratification/determinism and the Saltelli
design structure.

## 5. Pilot

Command (defaults from `[global_sensitivity]`):

```bash
julia --project=. scripts/global_sensitivity_sobol.jl
```

Pilot configuration: `N = 128`, `D = 8` → **2,304 model runs**, 80 sites (evenly
spread over model order), 1 simulation year, a uniform **+2 °C** warming stand-in,
bootstrap 200, seed 20261001. It ran in **14.7 s** after data loading (≈1.5 min
total including `prepare_ode_data`). Outputs in
`results/global_sensitivity_pilot/` (`sobol_indices.csv`, `parameter_ranges.csv`,
`saltelli_samples.csv`, `run_summary.txt`), clearly marked PILOT.

Selected pilot results (`S1` [95% CI] / `ST`):

| Response | `dispersal_multiplier` | `upstream_cost` | `carrying_capacity_scaling` |
|---|---|---|---|
| realised richness loss | S1 **0.554** [0.329, 0.781]; ST 0.616 | S1 **0.483** [0.318, 0.657]; ST 0.560 | S1 0.062 [0, 0.134]; ST 0.041 |
| native biomass | S1 ≈ 0.001; ST ≈ 0 | S1 ≈ 0; ST ≈ 0 | S1 1.087 [0.674, 1.507]; ST 1.107 |
| native richness | S1 −5.98 (unstable); ST 0.437 | S1 −1.54; ST 0.619 | S1 −2.79; ST 0.038 |

**Pilot caveats (must not be quoted as the production index):**

* `N = 128` is far below the production `N = 512`. Sobol estimators are
  high-variance at small `N`; first-order indices can fall outside `[0, 1]` and
  total-order indices can go negative. The continuous `realised richness loss`
  response is stable enough to read qualitatively (dispersal rate and upstream
  cost dominate); the discrete `native richness` response is not.
* The pilot uses a short static integration with a uniform +2 °C stand-in rather
  than the per-site daily 2026–2045 schedule, starts from the observed snapshot
  (no long burn-in) and a site subset. `heat_stress_k_multiplier` is therefore
  inert (warmed temperatures barely cross the empirical upper limits at this `k`).
* The `realised_richness_loss` baseline is the matched **unwarmed** control under
  the central parameter values (not the observed snapshot), matching the
  matched-control convention of the climate runs.

## 6. Production run and compute estimate

```bash
GUADEX_GS_MODE=production julia --project=. scripts/global_sensitivity_sobol.jl
```

Production uses all 775 sites and the configured 2026–2045 daily warming
schedule (warmest `obstacle_sensitivity.climate_models` pair), `N = 512`,
`D = 8` → **N(2D+2) = 9,216 model evaluations**, bootstrap 1000.

Wall-time estimate: the pilot evaluated 2,304 runs at 80 sites × 1 yr in 14.7 s
(≈6.4 ms/run). Scaling to 775 sites (~9.7×) and 20 years (~20×) gives ≈1.2 s/run,
so ≈ **3 h** (single core; budget **1–4 h**), plus one burn-in (cached) and the
~527 MB daily-forcing load. Reduce `GUADEX_GS_N` for a faster run with wider
confidence intervals. The production driver currently starts from the observed
state; the team should decide whether to substitute the spun-up initial condition.

Override knobs: `GUADEX_GS_N`, `GUADEX_GS_BOOTSTRAP`, `GUADEX_GS_SEED`,
`GUADEX_GS_MAX_SITES`, `GUADEX_GS_SIM_YEARS`, `GUADEX_GS_WARMING_C`,
`GUADEX_GS_HEAT_K`, `GUADEX_GS_CALIBRATE_HEAT`, `GUADEX_GS_OUTPUT_DIR`.

## 7. Files

* Added `src/global_sensitivity.jl` — LHS, Saltelli design, estimators, bootstrap,
  parameter transforms, basin response metrics.
* Added `scripts/global_sensitivity_sobol.jl` — thin runner.
* Added `test/test_global_sensitivity.jl` (registered in `test/runtests.jl`).
* Added `[global_sensitivity]` to `parameters.toml` (ranges + rationale).
* Modified `src/Guadex.jl` — include + exports.
* Modified `docs/climate_scenarios.md` — pointers to this note.
* Added `results/global_sensitivity_pilot/` — pilot artifacts (PILOT).
* Unchanged: the within-design ranking scripts and every existing default.

## 8. Residual uncertainty

* The ranges are a defensible **proposal**, not calibrated posteriors; they need
  team sign-off (see `docs/GuadeX_Correction_Plan_Sept2026.md`, E7 gate).
* Sobol indices assume independent inputs; correlated parameters (e.g. growth and
  dispersal trait co-variation) would need a different design.
* The response set is basin-level and deterministic; no stochasticity or
  extinction probability is modelled here (E8).
* `native richness` is a poor Sobol target at small `N`; production should report
  it with the same caution or use a smoothed/aggregated form.
* Production still needs the burn-in route decided (observed vs spun-up).
