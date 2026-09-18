# Warming and the fish community of the Guadalquivir basin

### A plain-language report on the 2026–2045 climate-scenario simulations (`climate_scenarios_k1x_burnin`)

**Who this is for.** Biologists, ecologists and managers who want to know what this
simulation ensemble found and what it did **not** find. No programming background is
assumed. Companion documents: `docs/climate_scenarios.md` (run setup) and
`docs/climate_scenarios_results_report.md` (technical review). All results are in
`results/climate_scenarios_k1x_burnin/`.

---

## 1. Summary in plain language

- **The question.** If water temperature rises as projected by 11 climate models
  (GCMs) under four emission scenarios (SSP1-2.6 to SSP5-8.5), how does the
  basin-wide fish community change between 2026 and 2045?
- **The starting point that governs everything.** Most native species are not in
  their thermal optimum today: the dominant cyprinids have optima ≈19–20 °C while
  the basin averages ≈16 °C, so warming *helps* them. Only brown trout has its
  optimum below current temperature, so warming *hurts* it (§2).
- **The basin-wide indicators barely move.** Mean native richness stays at ≈2.50
  species per site in every scenario, total biomass changes by <0.7 %, and the
  "richness-loss" indicator stays below 1 %. Scenario differences are smaller than
  differences between climate models, so no scenario ranking can be defended.
- **One clear biological signal:** the native cold-water specialist **brown trout
  (*Salmo trutta*, code ST)** declines. Its basin biomass falls by ≈**9 % even under
  the mildest scenario** and ≈**11 % under the warmest**, while the no-warming
  control stays flat (−0.2 %). At occupied sites the median population is 68–92 % of
  its 2026 level by 2045.
- **Why the averages miss this:** trout occurs at only ≈4 % of sites (≈34 of 775).
  A species that rare shifts basin *mean* richness by at most ~0.04 species even if
  it disappears.
- **Directionally robust, magnitude-uncertain.** The *sign* of each response is fixed
  by the position of the optimum relative to temperature; the *size* depends on the
  thermal niche width (σ), set to one value per species and never varied. The trout
  warming penalty roughly doubles across a plausible range of σ (§2).
- **Not shown:** climate-driven extinction, richness loss, invasive expansion, or any
  selection mode other than stabilizing selection toward a fixed optimum (§7).
- **Why this ensemble is interpretable at all:** unlike earlier runs it starts from a
  converged 373-year burn-in at observed density, so the control is stationary and the
  colonisation transient no longer hides the climate signal.

---

## 2. The starting point: most native species are not in their thermal optimum

The single most important feature of this system — and of the model — is that the
available temperature is **not** where most species perform best. Each species is
assigned a thermal optimum (the midpoint of its empirical temperature range) and a
thermal niche width σ (σ = range/6). Growth is maximal at the optimum and decays away
from it as a Gaussian; warming therefore shifts the whole community relative to fixed
optima rather than uniformly helping or hurting it.

Using the trait table and a 16.0 °C basin mean (T̄), the mismatch is systematic:

**Table 1 — Where the native species sit relative to current temperature.**

| species | thermal optimum T\* (°C) | niche width σ (°C) | maladaptation (T\*−T̄)/σ | suitability at T̄ | at +1.6 °C | dW/dT at T̄ (°C⁻¹) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| *Salmo trutta* (ST) | 12.0 | 2.67 | **−1.50** | 0.32 | 0.11 | **−0.18** |
| *Squalius pyrenaicus* (SP) | 16.5 | 2.83 | +0.18 | 0.98 | 0.92 | +0.06 |
| *Luciobarbus sclateri*, *Squalius alburnoides*, *Pseudochondrostoma willkommii*, *Iberochondrostoma lemmingii*, *Aphanius baeticus*, *Anaecypris hispanica* | 19.0 | 3.67 | +0.82 | 0.72 | 0.93 | +0.16 |
| *Iberochondrostoma oretanum* (IO) | 19.5 | 3.50 | +1.00 | 0.61 | 0.87 | +0.17 |
| *Cobitis paludica* (CP) | 20.0 | 3.33 | +1.20 | 0.49 | 0.78 | +0.18 |

Only brown trout sits on the **warm** side of its optimum (−1.5σ); every other native
is 0.2–1.2σ on the **cold** side. This is why ~1 °C of warming improves growth for the
dominant cyprinids and degrades it only for the cold-water specialist — and why the
basin averages, which they dominate, are insensitive.

![Native thermal niches against the site-temperature distribution.](../results/climate_scenarios_k1x_burnin/figures/diagnostics/thermal_niches.png)

*Figure 1 (entry figure) — Thermal niches versus available temperature. Grey bars: distribution of baseline site temperatures (≈11–19 °C); curves: Gaussian suitability of each native species; dotted lines: the 2026 mean and the mean after ~1 °C of warming. Most curves peak to the right of the site-temperature distribution — the community's maladaptation to current conditions is the mechanism behind the results that follow.*

### What selection mode is being simulated

The model simulates exactly **one** thermal selection mode: **stabilizing selection
toward a fixed, species-specific optimum**, with no evolution, no plasticity and no
acclimation. Thermal suitability is a symmetric Gaussian that multiplies the intrinsic
growth rate only; density dependence enters through the logistic total-biomass term
and a fixed interaction matrix. There is therefore no disruptive selection, no
frequency-dependent selection on the thermal trait, and no fluctuating selection
other than the seasonal cycle and the imposed warming trend. When the environment
shifts away from the optimum, the Gaussian becomes a **directional selection
gradient**, `∂logW/∂T = −(T − T*)/σ²`, whose sign is fixed by the position of the
optimum relative to ambient temperature and whose magnitude is fixed by σ. The
strength of stabilizing selection, `1/σ²`, is therefore not a nuisance parameter: it
sets both how strongly the population is selected toward its optimum *and* how
sensitive it is to climate change. Real populations may experience other selection
modes at the local level; quantifying that is beyond this analysis. What these runs
report is the community response conditional on this one, explicitly stated mode.

### How robust is the result to that assumption?

Only partially tested, and this is the main caveat of the report. **σ was a single
deterministic value per species (σ = range/6) and was not varied** in this ensemble.
The only structural thermal sensitivity implemented is an optimum-position sweep
(the staged E4 experiment), and these runs used the midpoint (fraction 0.5) only. The
qualitative result is nonetheless robust, because the **sign** of each response is set
by (T\* − T) and not by σ: warming favours the cold-marginal cyprinids and penalises
the only warm-marginal native, trout. The **magnitude** is strongly σ-dependent.
Holding the trout optimum at 12 °C and changing only σ:

**Table 2 — Sensitivity of the trout warming response to niche width σ.**

| σ (°C) | interpretation | suitability at +1.6 °C / at T̄ |
| ---: | :--- | ---: |
| 2.67 | range/6 — the value used | 0.33 |
| 3.00 | +12 % wider | 0.42 |
| 4.00 | range/4 | 0.61 |

The warming penalty roughly doubles across a 1.5-fold range of σ. A related assumption
to verify: the code comment for σ = range/6 states that the empirical limits then sit
at ±2σ (≈95 % of the niche), but they in fact sit at ±3σ (≈99.7 %). If the reported
thermal range is meant to span ±2σ, σ should be range/4, which would weaken every
modelled warming response. **Until a σ sweep is run, the trout decline should be
reported as a directionally robust but magnitude-uncertain result**, and the selection
gradient `−(T̄ − T\*)/σ²` in Table 1 is the quantity to report alongside the
abundance change, because it makes the assumed selection strength explicit.

---

## 3. What was simulated

| item | setting |
| :--- | :--- |
| System | Guadalquivir basin: 775 sites, 774 water bodies, one basin (`ES050`) |
| Species | 24 modelled: 10 native residents, 10 invasive, 4 migratory (reported separately) |
| Horizon | 2026–2045, daily integration, yearly reported values |
| Climate forcing | per-site **daily** water temperature, 11 GCMs × 4 SSPs, anomalies against 1986–2005 |
| Runs | 44 scenario runs (11 × 4) + 1 no-warming control = **45** |
| Initial state | converged 373-year burn-in under the seasonal baseline climate |
| Carrying capacity | 1 × observed density (mean K ≈ 86.6 vs mean observed ≈ 84.4) |
| Heat stress | one shared slope, calibrated so the worst baseline species/site loses ≤5 % yr⁻¹ |

Each site holds 24 species whose density changes through logistic growth toward a
site carrying capacity, modified by biotic interactions, dispersal along the river
network (distance, upstream cost, dam passability) and the Gaussian thermal niche
described in §2. Above each species' empirical upper limit an extra heat-stress
mortality applies. This is a directional metacommunity model for comparing scenarios
and mechanisms, not a calibrated abundance forecast.

**Why the burn-in matters.** In earlier GuadeX runs the community started far below
carrying capacity and spent the window filling in, producing a rising-richness trend
unrelated to climate. Here the baseline is at equilibrium first, so change over
2026–2045 can be attributed to warming rather than recovery.

---

## 4. How to read the outputs

Every run reports **five indicators at four nested scales** (sampling point →
sub-basin → water body → basin), each a mean over the units within a level.

| indicator | what it means |
| :--- | :--- |
| Native richness | number of the 10 native species above density 0.1 at a site, averaged over sites |
| Invasive richness | same for the 10 invasive species |
| Native richness loss (relative) | `max(0, 1 − current richness / 2026 equilibrium richness)`; labelled "extinction risk" in some figures. It can only exceed zero if basin-average richness falls below its starting value |
| Total biomass | summed density of all 24 species per site |
| Projected water temperature | basin-mean water temperature |

The standard figure manifest is
`results/climate_scenarios_k1x_burnin/figures/figure_inventory.csv`.

![Final-year summary across all four spatial levels and five indicators. Markers are GCM means per scenario; bars show 25–75 % and 10–90 % spread across GCMs.](../results/climate_scenarios_k1x_burnin/figures/summary_2045.png)

*Figure 2 — 2045 summary. The temperature panel separates scenarios; richness, biomass and richness-loss do not.*

---

## 5. Results

### 5.1 How much warming the runs delivered

| scenario | basin ΔT in 2045 vs 1986–2005 | range across GCMs | added over 2026–2045 |
| :--- | ---: | ---: | ---: |
| control | 0.00 °C | — | 0.00 °C |
| SSP1-2.6 | +0.73 °C | +0.41 to +1.10 | +0.74 °C |
| SSP2-4.5 | +0.91 °C | +0.53 to +1.60 | +0.87 °C |
| SSP3-7.0 | +0.70 °C | +0.43 to +1.16 | +0.88 °C |
| SSP5-8.5 | +0.97 °C | +0.15 to +1.63 | +0.97 °C |

The two temperature columns are **not the same quantity** and should not be
conflated. "Added over 2026–2045" is the imposed warming curve, anchored at zero in
2026; it is monotonic across emission scenarios (SSP1-2.6 +0.74 °C → SSP5-8.5
+0.97 °C). "Basin ΔT in 2045 vs 1986–2005" is the realised daily forcing expressed
against the fixed historical baseline; because each SSP's downscaled product already
carries a different anomaly in 2026, this column is **not directly comparable across
emission scenarios** and is not monotonic (SSP3-7.0 sits below SSP2-4.5 for this
model set). Use the added-warming column for scenario comparisons.

Even there, this is **early-warming territory**: the most aggressive emission
scenario adds only ~1 °C by 2045, because the 20-year window samples just the start of
the century-scale rise. The spread across the 11 climate models overlaps the
differences between scenarios, so at this horizon the climate model matters as much
as the emission scenario.

### 5.2 The basin-wide indicators do not respond

Values for 2045 relative to the 2026 equilibrium (ensemble mean across GCMs):

| indicator at basin scale | control | SSP1-2.6 | SSP2-4.5 | SSP3-7.0 | SSP5-8.5 |
| :--- | ---: | ---: | ---: | ---: | ---: |
| native richness (relative) | 0.9993 | 0.9993 | 0.9990 | 0.9991 | 0.9993 |
| native biomass (relative) | 1.0078 | 1.0052 | 1.0050 | 1.0051 | 1.0046 |
| total biomass (relative) | 1.0065 | 1.0038 | 1.0036 | 1.0038 | 1.0032 |
| richness-loss indicator | 0.0074 | 0.0077 | 0.0078 | 0.0078 | 0.0077 |

All changes are well below 1 %, and the control sits inside the spread of every
scenario. The residual drift in the control (+0.65 % total biomass over 20 years) is
a slow internal cycle of the model, not a warming effect; climate effects must be
read as a contrast *against the control*, not against zero.

![Basin-scale time series, one line per SSP with a 10–90 % band across GCMs. Only the temperature panel separates the scenarios.](../results/climate_scenarios_k1x_burnin/figures/comparison/across_scenarios_basin.png)

*Figure 3 — Across-scenario comparison at basin scale. The four SSP lines overlap for richness, richness loss and biomass; the black line is the control.*

![SSP5-8.5 ensemble at basin scale; line is the mean across GCMs, bands are 10–90 % and 25–75 % across GCMs.](../results/climate_scenarios_k1x_burnin/figures/ensemble/ensemble_ssp585_basin.png)

*Figure 4 — The warmest scenario ensemble. The GCM band is wide relative to the 20-year trend.*

### 5.3 Why the community does not lose richness as it warms

The explanation is §2: the mean richness and biomass are dominated by warm-adapted
cyprinids whose optima sit *above* present temperature, so warming raises their growth
suitability while reducing trout's. Because those species hold most of the sites and
most of the biomass, the community-wide averages are effectively blind to the one
species that is genuinely at risk.

![Per-site richness distribution and site-by-site change, control run.](../results/climate_scenarios_k1x_burnin/figures/diagnostics/community_filling.png)

*Figure 5 — Community turnover at site level (control). The community is at equilibrium: 7 sites gained native species, 751 were unchanged, 17 lost species. Every warming scenario shows the same picture (7–8 gained, 15–18 lost), so per-site richness turnover is model noise, not a climate response.*

### 5.4 The signal that is there: brown trout

Basin-scale change in ST biomass between 2026 and 2045:

| scenario | mean change across GCMs | range across GCMs |
| :--- | ---: | ---: |
| control | −0.2 % | — |
| SSP1-2.6 | −9.0 % | −4.2 to −14.6 % |
| SSP2-4.5 | −9.7 % | −5.8 to −13.9 % |
| SSP3-7.0 | −9.5 % | −6.2 to −14.2 % |
| SSP5-8.5 | −10.9 % | −6.8 to −18.5 % |

At the 50 sites with non-zero trout density at the start, the **median population
falls to 92 % (SSP1-2.6), 83 % (SSP2-4.5), 75 % (SSP3-7.0) and 68–81 % (SSP5-8.5)**
of its 2026 value. Unlike the basin averages, every scenario is separated from the
control. The exposure diagnostics show why:

| ST thermal exposure | baseline 1986–2005 | SSP1-2.6 | SSP5-8.5 (MRI) | SSP5-8.5 (UKESM) |
| :--- | ---: | ---: | ---: | ---: |
| max days/yr above 20 °C | 131 | 145 | 162 | 167 |
| max exceedance energy | 1 582 | 2 145 | 3 977 | 3 433 |
| mean days/yr above 20 °C (all site-years) | 75 | 79 | 101 | 97 |

Trout is already stressed under the baseline (up to 131 days/yr above its limit at
the most exposed site), and warming lengthens that exposure and more than doubles its
intensity. This is physiologically coherent: a cold-water specialist already living
1.5σ on the warm side of its optimum is squeezed further. Note that exposure is a
**habitat thermal diagnostic** reported at every site, including sites without trout,
so the mean is not a trout-population average.

![Forcing versus community response, and the final-year warming–response relationship.](../results/climate_scenarios_k1x_burnin/figures/diagnostics/forcing_and_response.png)

*Figure 6 — Warming forcing is separated by SSP (top left) while richness and biomass responses overlap. Bottom row: richness is essentially uncorrelated with warming (r = −0.08), biomass is weakly negative (r = −0.43), and the richness-loss indicator rises in every run including the control.*

### 5.5 Individual runs and level detail

Each `per_run/` figure shows all five indicators at all four scales and is mainly a
quality-control view for checking that a run is coherent; at nested levels the
indicators look even flatter than at basin scale because averaging removes
site-to-site variability.

![One run in full: SSP5-8.5 / MRI-ESM2-0. Rows are indicators, columns are levels.](../results/climate_scenarios_k1x_burnin/figures/per_run/ssp585__MRI-ESM2-0.png)

*Figure 7 — A single warm run. The clear trend is the temperature rise; richness and biomass are flat or very slightly negative.*

---

## 6. What these findings show

1. **A converged-baseline ensemble now exists and behaves as intended** — the
   control is stationary, so the design can detect a climate signal if one exists.
2. **The community is systematically maladapted to present temperature** (most natives
   are below their optimum), which is the mechanism that determines who wins and who
   loses under warming (§2, Table 1).
3. **Under 0.7–1.6 °C of warming to 2045 the basin's aggregate community composition
   is essentially unchanged**, and scenario differences are within model spread.
4. **A real, directionally reproducible negative signal exists for *Salmo trutta*:**
   present in all 44 scenario runs, absent from the control, stronger with warming,
   and corroborated by independent thermal-exposure diagnostics.
5. **Aggregate richness is the wrong indicator for this question.** A rare
   specialist's decline is diluted to invisibility by abundant warm-adapted species;
   species-level abundance, exposure and selection gradients are the informative
   metrics.
6. **The method separates "no signal" from "no data":** the flat aggregates are
   explained mechanistically by the thermal niches of the dominant species, not by a
   failed run.

## 7. What these findings do not show

1. **No climate-driven extinction or richness loss.** The model has no general
   mortality process; heat stress is the only warming-driven mortality and it only
   activates above a species' upper limit, while richness can only fall below the 0.1
   density threshold. That happens in 15–18 sites per run — including the control.
2. **No scenario ranking.** The GCM spread exceeds emission-scenario differences by
   2045, so these runs cannot say SSP5-8.5 is worse than SSP1-2.6 in any metric. The
   result is "warming harms trout", not "this emission scenario costs X %". A 2070
   horizon or stronger forcing separation is needed. (Absolute ΔT against 1986–2005
   is also not comparable across SSP products, because they carry different 2026
   offsets; see §5.1.)
3. **The "extinction-risk" metric is not a risk measure.** It is `1 − current / initial
   richness` floored at zero, so it is ≤1 % here by construction: a realised loss, not
   a probability or danger of future extinction.
4. **The quasi-extinction fraction is unusable and is not a result.** It reaches
   ≈0.75 even in the control because the threshold is a fraction of each species' own
   baseline and flags naturally rare species. Ignore it until the definition is revised.
5. **Only brown trout responds measurably.** Warm-adapted natives and invasive species
   are unchanged; no invasive-expansion signal was detected. Other species (e.g.
   *Squalius pyrenaicus*, *Gobio lozanoi*, *Oncorhynchus mykiss*, *Liza ramada*,
   *Mugil cephalus*) have upper limits of 25–30 °C and essentially never exceed them,
   so the heat-stress mechanism is inert for them.
6. **The thermal selection mode is assumed, not tested, and no other mode is
   represented.** The model has stabilizing selection toward a fixed optimum and no
   evolution, plasticity or acclimation; disruptive, frequency-dependent and
   (beyond seasonality) fluctuating selection are absent. Warming in this model
   therefore cannot drive adaptation or evolutionary rescue.
7. **Selection strength (σ) was not explored.** σ is a single value per species
   (range/6) and only its sign-level consequence is robust; the magnitude of the trout
   decline changes substantially with σ (Table 2). The range/6 convention also places
   the empirical limits at ±3σ even though the code comment describes ±2σ — a
   discrepancy to resolve before publication.
8. **Nothing here is a spatial projection.** Every site receives the same basin
   warming anomaly (elevation scaling is off), so the model does not resolve which
   reaches warm most.
9. **Key ecological processes are outside the model:** flow and drought regime, water
   abstraction, land-use/habitat deterioration (the habitat index is a static
   placeholder), age structure, fishing, and any organisms beyond the 24 modelled
   fish — including Mediterranean drought and flow alteration, plausibly first-order
   in this basin.
10. **The trout decline includes pre-existing stress.** The baseline already puts trout
    above 20 °C on up to 131 days/yr; warming worsens an already marginal position,
    it is not a pristine-baseline climate effect.
11. **One classification issue to verify.** `parameters.toml` lists `AAL` as
    "migratory", but the source trait tables identify it as *Alburnus alburnus*
    (Alburno), an exotic/invasive cyprinid. If that is an error, the invasive metrics
    currently omit one invasive species. It does not affect the trout result.
12. **A residual control drift remains** (total biomass +0.65 % over 20 years). The
    burn-in converged to 5 × 10⁻⁵ yr⁻¹ but a slow internal cycle persists, so
    differences smaller than ~1 % are not interpretable as climate effects.

---

## 8. Figure and output guide

All figures are under `results/climate_scenarios_k1x_burnin/figures/`. The manifest
with one row per standard figure is `figures/figure_inventory.csv`; the three
`diagnostics/` figures are additional.

| category | files | content |
| :--- | ---: | :--- |
| `diagnostics/thermal_niches.png` | 1 | native optima vs site-temperature distribution — **entry figure** (Fig. 1) |
| `summary_2045.png` | 1 | final-year summary, 4 levels × 5 indicators, scenario mean ± GCM spread (Fig. 2) |
| `comparison/across_scenarios_<level>.png` | 4 | one level, one line per SSP with 10–90 % GCM band (Fig. 3) |
| `ensemble/ensemble_<scenario>_<level>.png` | 20 | one scenario × one level, 5 indicators, GCM bands (Fig. 4) |
| `per_run/<scenario>__<gcm>.png` | 45 | one run, 5 indicators × 4 levels (Fig. 7) |
| `diagnostics/forcing_and_response.png` | 1 | forcing vs response, warming–response scatter, richness-loss trajectory (Fig. 6) |
| `diagnostics/community_filling.png` | 1 | initial vs final per-site richness and site-by-site change (Fig. 5) |

Per-run numeric outputs are in `<run>/export/levels/`:
`level_basin.csv`, `level_subcatchment.csv`, `level_water_body.csv`,
`level_sampling_point.csv`, `species_timeseries.csv` (site × species × year),
`quasi_extinction_summary.csv` and `exposure_sites.csv`. Each run's configuration is
in `<run>/export/run_metadata.json`; `runs_index.csv` summarises all 45 runs.

## 9. What we recommend next

1. **Explore the strength of stabilizing selection directly.** Add a σ-multiplier
   sweep (e.g. 0.5, 0.75, 1.0, 1.25, 1.5, 2.0× per species) to the climate driver,
   and re-run the optimum-position sweep (E4) under this converged-baseline design.
   A small σ × optimum-fraction grid for brown trout is the highest-value experiment.
2. **Calibrate σ and the optimum from empirical thermal-performance curves** (T_opt,
   CT_max) instead of the midpoint/range-6 heuristic, so the assumed selection
   strength is auditable; resolve the range/6 vs ±2σ discrepancy in the code comment.
3. **Report brown trout (and other cold-water specialists) with species-level
   abundance, exposure and the selection gradient** `−(T̄ − T\*)/σ²`, not basin
   richness — that is where the climate signal lives.
4. **Extend the horizon to 2070** (or use the daily per-GCM series to its full extent)
   so emission scenarios separate; the current window cannot rank scenarios.
5. **Add the processes most likely to dominate Mediterranean fish vulnerability:** flow
   regime and drought, water abstraction, and a real habitat-degradation input.
6. **Revise or retire the quasi-extinction threshold** so rare-but-stable species are
   not flagged at baseline.
7. **Verify the `AAL` species classification** before publishing invasive-species
   metrics.

---

*Generated from `results/climate_scenarios_k1x_burnin` (45 runs; 775 sites;
2026–2045). Run setup: `docs/climate_scenarios.md`. Engineering review of earlier
ensembles and the burn-in / carrying-capacity change:
`docs/climate_scenarios_results_report.md`.*
