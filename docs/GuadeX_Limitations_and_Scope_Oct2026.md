# Limitations and Scope of Inference

**GuadeX — October 2026 report.** This note consolidates the standing limitations of the
current model and its projections into a single numbered statement, written to be read as the
Limitations section of the forthcoming October 2026 report. It supersedes the scattered
caveats in `docs/climate_scenarios.md`, `docs/GuadeX_Correction_Status_Sept2026.md` and
`docs/GuadeX_Correction_Plan_Sept2026.md` for the items below, and it records how the team's
**C6, C3, E21, IF and E8** decisions are represented. No model code or `.tex` report is changed
by this document.

---

## Limitations

1. **Thermal niches are distributional, not physiological (C6).** The thermal response of each
   species is represented by a symmetric Gaussian whose optimum is the midpoint of the reported
   thermal range and whose breadth (σ) is fixed at one sixth of that range. This construction
   is a first-order sensitivity baseline, as the review explicitly accepts, and it is retained
   here rather than treated as a calibrated physiological niche. Because every baseline site in
   the basin lies at or below 18.37 °C, and because the warming delivered over the simulation
   horizon is mild (approximately 1 °C in the central projections), the bluntness of the
   symmetric curves has little leverage on the present thermal results. That insensitivity is
   not a general property of the model: it is a consequence of the mild forcing and the cool
   baseline. Under stronger warming, and for warm-adapted species whose true growth optima lie
   above the range midpoint, the symmetric-Gaussian approximation becomes a materially limiting
   assumption. Replacement of the range midpoint with literature-based optimum-for-growth and
   critical thermal maximum (CTmax / UILT) values, and replacement of the symmetric Gaussian
   with an asymmetric physiological performance curve (for example the Sharpe–Schoolfield
   formulation), are therefore reserved as future work and are required before the thermal
   results can be interpreted under strong-warming scenarios.

2. **The projections are interim and uncalibrated (C3).** The projections follow the review's
   interim route: they are initialised from the observed 2006–2009 community and advanced
   through a short three-year baseline spin-up, rather than from a long model-generated burn-in.
   Every scenario is reported as a scenario-minus-control difference matched to the no-warming
   control for the same general-circulation model (GCM), so that any residual transient common
   to the scenario and its control cancels. The intrinsic growth rates (r), carrying capacities
   (K), interaction coefficients and dispersal parameters are **not** calibrated, and the model
   does **not** reproduce the observed community to a pre-registered standard. The results are
   consequently conditional and relative, not absolute forecasts: they should be read as
   internally consistent differences between a warming pathway and its matched control, and not
   as predictions of the absolute state of the community. Full simulation-based calibration
   (approximate Bayesian computation with sequential Monte Carlo, ABC-SMC) against
   pre-registered true-skill-statistic (TSS) and occupancy criteria, with spatial block
   cross-validation and posterior predictive checks, remains future work.

3. **The air–water warming response is a same-day calibration with a lagged uncertainty band
   (E21).** The primary air–water warming relationship is calibrated on same-day air
   temperature, giving a slope of approximately 0.510 across Spain and 0.436 for the
   Guadalquivir. Lagged air fits (3-, 7- and 30-day means) yield higher slopes, up to
   approximately 0.65 Spain-wide, but they rest on only about three in-sample lag-eligible sites
   (n = 104, one year of data); a lagged refit therefore cannot serve as a valid primary
   calibration. The lagged relationship is instead expressed as an uncertainty band on the
   projected warming, spanning ×0.63 to ×1.28 of the central value
   (`guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv`). Within the Guadalquivir
   subset the lagged slopes are lower than the same-day slope, so the upper end of the band may
   not apply to the basin, and the band should be treated as a conservative sensitivity
   envelope rather than as a calibrated range for this system.

4. **Barrier passability rests on a uniform structural assumption (IF).** Barrier passability is
   currently modelled using a uniform structural assumption (0.1 upstream, 0.5 downstream)
   pending full standardisation of the regional obstacle inventory's qualitative passability
   index (IF). The IF index is not decoded in the present runs; using it prematurely would
   silently mis-scale every barrier. This assumption is a deliberate, documented stand-in for
   the qualitative index and is a first-order control on connectivity that should be revisited
   once the inventory's passability semantics are standardised.

5. **The richness-loss metric is a deterministic, realised contraction, not an extinction
   probability (E8).** The richness-loss metric reports a realised, threshold-based contraction
   of the community under the deterministic trajectory; it is not a probability of extinction
   for any species. The model contains no demographic or environmental stochasticity and is
   therefore run without replicates. Consequently, the metric describes the single deterministic
   trajectory implied by the parameters and forcing, and it cannot be interpreted as, or
   converted into, a quasi-extinction or extinction probability. Adding demographic and
   environmental stochasticity with replicate simulations and species-specific quasi-extinction
   thresholds remains future work, and until it is added the reported metric should be described
   only as realised richness loss.

6. **Absence of hydrological dynamics, drought and water abstraction.** The model does not
   represent river discharge, flow intermittency, drought or water abstraction. Habitat
   suitability enters through a static environmental index and through the site temperature
   series, but there is no dynamic hydrological driver and no representation of dry reaches,
   reduced runoff or abstraction-driven habitat loss. Warming is therefore the only
   climate-linked process that varies over the projection, and the results should not be read as
   projections of the combined effect of climate and water management.

7. **No age, stage or body-size structure and no evolution or acclimation.** Individuals are
   represented by a single aggregated abundance or biomass variable per species and site; there
   is no age, stage or body-size structure and no representation of ontogenetic shifts in
   habitat or thermal tolerance. There is likewise no evolution or acclimation: thermal traits
   are fixed, so populations cannot adapt to warming, and the model cannot represent
   range shifts that depend on adaptation. Both are plausible buffering mechanisms under
   sustained warming and their omission means the projections may overstate thermal stress
   relative to a model that allowed adaptation.

8. **The river network is an on-path reconstruction, not a flow-directed channel layer.** The
   connectivity graph is reconstructed by an on-path parent rule with an explicit
   minimum-spanning-tree root join, because the regional channel layer carries no flow-direction
   attribute or elevation. It is a topology defined by on-path proximity rather than a
   flow-directed channel network derived from the hydrography. The reconstruction is acyclic,
   deterministic and validated by topology, but its edges do not encode true flow direction, and
   the root-joining step produces a small number of flow-undirected connections. Connectivity
   results should be interpreted relative to this reconstructed topology.

9. **Biomass is reported in model units without an absolute scale.** Abundance and biomass are
   expressed in model units. They are internally consistent across sites, species and scenarios
   and are suitable for relative comparison, but they carry no absolute physical scale and
   cannot be read as real standing stocks or densities without an external calibration.

10. **Short projection horizon (2026–2045).** The projections span 2026–2045, a twenty-year
    horizon. This is short relative to the generation times of the long-lived species and to the
    timescales of community reorganisation, and it truncates the analysis before the
    differentiating effects of the higher-emission pathways are fully expressed. Statements
    about scenario divergence therefore refer to differences that are emerging by 2045, not to
    long-run outcomes.

11. **Climate-model spread overlaps pathway differences.** The scenario differences are small
    relative to the spread among the individual GCMs within any single pathway, and the two
    sources of variation overlap. The central warming delivered by the lower pathways cannot be
    separated from pathway differences on the strength of a single GCM, which is why the results
    are reported as ensemble statistics and as scenario-minus-control differences matched by
    GCM, and why no ordering of pathways should be inferred where the GCM spread overlaps it.

---

## Scope of inference

Taken together, these limitations define a narrow and explicit scope. The present results are a
conditional, internally consistent, deterministic sensitivity assessment of the effect of a
mild, prescribed thermal forcing on a lightly structured and uncalibrated metacommunity model,
reported against matched no-warming controls. They are informative about the sign, relative
magnitude and spatial pattern of the modelled thermal response, and about how that response
behaves within the modelled assumptions. They are not absolute forecasts, not extinction
probabilities, and not integrated assessments of climate together with hydrology, management or
adaptation. Each of the limitations above names the specific mechanism that would have to be
added, calibrated or standardised before a wider interpretation could be supported.

---

## Source artifacts and decision record

* Thermal niches (C6), the interim projection route (C3), the air–water uncertainty band (E21),
  barrier passability (IF) and the richness-loss metric (E8) are decided and documented in
  `docs/GuadeX_Correction_Status_Sept2026.md` (§4 Team decisions, §5 Known residuals) and
  `docs/climate_scenarios.md`.
* The air–water slope uncertainty band is tabulated in
  `guadex_tw/outputs/tables/air_water_slope_uncertainty_band.csv`.
* The corrected-run configuration that selects the C3 interim route and the team's biological
  options is `parameters_climate_scenarios_corrected.toml`.
