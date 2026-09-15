# Climate-scenario results: why richness keeps rising

**Scope:** review of `run_climate_scenarios.jl` and its outputs under
`results/climate_scenarios` (44 runs = 4 SSPs x 11 GCMs, 2026-2045, daily step,
monthly saves).

**Short answer:** the rising richness is **not a warming response**. It is the
model filling in from a very sparse initial community toward its coexistence
equilibrium. The climate signal over 2026-2045 is small, scenario-insensitive,
and (for the species that make up the richness metric) points *toward* their
thermal optima, so warming slightly *increases* richness instead of reducing it.
No implementation bug produces this; the pattern follows directly from the
model's structure and initial conditions.

> **Update (modelling-improvement plan).** The runs reviewed here were produced
> by the pre-plan model: annual-mean forcing anchored at 2026, no mortality
> term, `K = 10x` observed, and a native list that omitted the cold-water
> keystone `ST`. The reasons above are now addressable, and opt-in by design:
> WP0 adds `ST` to the native metric, WP1/WP2 deliver per-site daily
> seasonality against a fixed baseline, WP3 adds heat-stress mortality above the
> empirical limits, WP4 starts from a spun-up equilibrium and rebases the loss
> metrics, and WP5 adds abundance/occupancy/quasi-extinction and exposure
> diagnostics. See `docs/climate_scenarios.md` ("Modelling-improvement features")
> and the staged driver `run_climate_experiments.jl`. The default settings still
> reproduce the runs analysed below.

---

## 1. What the outputs show

Mean basin-level values (runs_index + `level_basin.csv`):

| metric | 2026 (initial) | 2045 | change |
| :--- | ---: | ---: | ---: |
| native richness / site (of 9) | 1.27 | 2.73 | +115% |
| invasive richness / site (of 10) | 0.26 | 0.32 | +22% |
| total richness / site (of 24) | 1.62 | 3.15 | +94% |
| total biomass / site | 84.4 | 647.2 | +667% |
| native extinction-risk metric | 0.0 | 0.0126 | flat after year 1 |

Per site (2026 -> 2045, representative run):

* **546 sites gained** native species, 191 were unchanged, only **38 lost**
  species.
* **345 sites started with zero** native species; this fell to **57**.
* The metric is monotone and still rising at the end of the run (total biomass
  is roughly 3/4 of the mean carrying capacity implied by the observed
  densities; it had not equilibrated by 2045).

Warming forcing delivered over the horizon (median across GCMs):

| SSP | min | median | max | final native richness |
| :--- | ---: | ---: | ---: | ---: |
| ssp126 | 0.28 | 0.71 | 1.09 °C | 2.719 |
| ssp245 | 0.55 | 0.76 | 1.37 °C | 2.720 |
| ssp370 | 0.59 | 0.79 | 1.22 °C | 2.723 |
| ssp585 | 0.65 | 0.91 | 1.31 °C | 2.723 |

So even ssp585 warms less than ~1.4 °C by 2045, the scenario separation is
tiny, and the final richness differs only in the third decimal.

![forcing and response](../results/climate_scenarios/figures/diagnostics/forcing_and_response.png)

The top row (left to right) shows that the temperature forcing *is* separated by
SSP while the richness and biomass responses are visually identical (the four
scenario lines overlap). The bottom row shows the counterintuitive result: more
warming is associated with slightly *higher* richness (Pearson r = 0.91) and
biomass (r = 0.98), and the "extinction risk" metric stays near zero.

---

## 2. Why richness rises

There are three independent reasons, in order of importance.

### 2.1 The community starts far below carrying capacity and fills in

`u0` is the observed density matrix, and carrying capacity is set to **10x the
observed total density per site** (`build_carrying_capacity`). The integrated
model is logistic growth with dispersal:

```
dU[i,s] = N[i,s] * ( r_effective * clamp(1 - ΣU/K, -1, 2) + interaction ) + dispersal
```

* Every species always has a **strictly positive** intrinsic growth rate: full
  rate where observed, **10% of it where absent**
  (`build_intrinsic_growth_rates`).
* There is **no explicit mortality / extinction term**. Warming can only scale
  growth down; it cannot drive a species to zero.
* Presence/richness counts a species when density > 0.1, so as soon as a
  dispersal event plus growth lifts a colonist above 0.1 it is counted.

Starting from a mean of only 1.27 of 9 native species per site and 84 biomass
vs a mean K of ~842, the system spends the whole 20-year window growing into its
attractor. This is the dominant term in the trend.

![community filling](../results/climate_scenarios/figures/diagnostics/community_filling.png)

### 2.2 The warming forcing is small and barely scenario-dependent

The projection table has **three windows only** (2021-2040, 2036-2055,
2041-2070). `basin_warming_curve` anchors zero warming at `start_year` (2026)
and the window medians at their midpoints (2030.5, 2045.5, 2055.5). Over a
2026-2045 horizon this means:

* The 2021-2040 window contributes very little (its midpoint is only 4.5 years
  after the forced zero), and the 2041-2070 window is only reached at 2045.5.
* The simulated interval therefore samples just the *early, slowly diverging*
  part of the projections, so all four SSPs look similar (0.3-1.4 °C).
* Because warming is anchored to zero in 2026, warming already realised before
  the start year is discarded.

This is a deliberate, documented choice (`docs/climate_scenarios.md`, steps 2-3)
and not a bug, but it fundamentally limits what the 2026-2045 ensemble can show.
To expose scenario differences you would need a later `end_year` (e.g. 2070) or
the daily per-GCM series rather than the window table.

### 2.3 Warming moves the modelled species toward their optima

The richness metric counts the 9 species in `[species].native`. Their thermal
optima (from `caracteristicas_peces_Guadalquivir_03-04-2018.csv`, sigma =
range/6) are:

* AB, AH, PW, LS, SA, IL -> optimum **19.0 °C**
* IO -> **19.5 °C**, SP -> **16.5 °C**, CP -> **20.0 °C**

Mean baseline site temperature is **16.0 °C** (range 11.2-18.4 °C), and the
median 2045 temperature is about **17.2 °C**. The Gaussian thermal filter
`exp(-(T-opt)^2 / 2σ^2)` therefore *increases* for the dominant opt-19 group as
sites warm, which is exactly why warmer runs end with marginally more richness
and biomass. The only native species moving away from its optimum is the rare
cold-water specialist *Salmo trutta* (optimum 12 °C) — which is **not** in the
native richness list (see caveats). This is the mechanism behind the positive
r = 0.91 / 0.98 correlations.

![thermal niches](../results/climate_scenarios/figures/diagnostics/thermal_niches.png)

---

## 3. Implementation review

**Correct / working as intended**

* Metrics are defined consistently in `src/outputs.jl:345-374`:
  `native_richness = count(density > 0.1)`, and
  `native_extinction_risk = max(0, 1 - richness / richness(t=0))`.
  The risk metric is a *loss relative to the initial observed state*, floored at
  zero, so a basin that gains species can never register risk.
* The scheduled temperature actually reaches the ODE: projected
  `temperature_c` in the outputs rises from 16.0 to 17.2 °C and matches the
  warming matrix; `test_ode.jl` also covers the zero-anomaly and warming cases.
* The warming curve, warming matrix and file exports match the documented
  behaviour; the resume/skip logic and the run index are consistent.
* Ensemble spread across GCMs is genuine; the bands in the diagnostics are
  10-90% across GCMs.

**Caveats / things to decide**

1. **`ST` is missing from the native list.** `parameters.toml` lists
   `native = ["AB","AH","SP","PW","LS","SA","IL","CP","IO"]`, but *Salmo trutta*
   (ST) is a native cold-water keystone and is also not in the invasive list.
   `AA`, `AAL`, `LR`, `MC` are likewise unclassified. Native richness therefore
   undercounts the actual native assemblage, and the 38 sites that lose richness
   are plausibly the cold headwater sites where ST/GL/OM dynamics dominate.
   Worth confirming whether this is intentional.
2. **No extinction process.** With strictly positive growth and no explicit
   mortality, richness loss can only be produced indirectly by competition
   pushing a species below 0.1. Under 1-1.4 °C of warming this never happens at
   basin scale. Any claim about climate-driven extinctions needs either a longer
   horizon, stronger warming, or a temperature-dependent mortality term.
3. **Horizon vs forcing.** As in 2.2, the 2026-2045 window cannot separate SSPs.
   Consider `end_year = 2070` (or a daily projection series) for a scenario
   comparison that is actually informative.
4. **Carrying capacity from a single snapshot.** K = 10x observed density means
   the initial state is far from equilibrium by construction, which is the main
   driver of the trend. If the intent is to start at equilibrium, initialise at
   (a fraction of) K instead of the observed snapshot.
5. **Biomass is still rising at 2045**, so the final year is a transient, not an
   equilibrium; do not read long-term projections from it.

---

## 4. New diagnostic figures

Added `src/climate_diagnostics.jl` (exported as `plot_climate_diagnostics`) and
`scripts/plot_climate_diagnostics.jl`. They read the exported CSVs only and
never re-run the model:

```bash
julia --project=. scripts/plot_climate_diagnostics.jl
```

Outputs under `results/climate_scenarios/figures/diagnostics/`:

| file | message |
| :--- | :--- |
| `forcing_and_response.png` | separated warming forcing vs overlapping response; final-year warming-response scatter; flat risk metric |
| `thermal_niches.png` | native thermal optima vs site-temperature distribution and the end-of-horizon warming |
| `community_filling.png` | initial vs final per-site richness; 546 sites gained, 38 lost |

The diagnostics are also generated automatically at the end of
`run_climate_scenarios.jl` when `make_figures = true`.

---

## 5. Recommended next steps

1. Confirm the intended native/invasive classification of `ST`, `AA`, `AAL`,
   `LR`, `MC`.
2. Re-run with `end_year = 2070` (and/or the daily ensemble series) to obtain a
   forcing signal that separates SSPs.
3. If the goal is climate *impact*, add a temperature-dependent mortality or an
   extinction threshold, otherwise richness cannot decline under warming.
4. Decide whether to start at (or near) carrying capacity so the reported trend
   is a climate response rather than a colonisation transient.
