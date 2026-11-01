# 1. Summary

- **Most of the response is confirmed.** The corrections are implemented, the test suite passes and nearly every number we checked reproduces exactly (section 3).
- **We accept the errors you identified in our review** (section 2). None of them changes a diagnosis, as you kindly note.
- **A few statements in the response do not match the code**, and we would kindly ask you to revisit them: the control is shared, not matched by GCM (4.1); the slope we quoted is the one that drives the projections (4.2); and the long burn-in years do not belong to the corrected runs (4.3).
- **Some points deserve further discussion**: a residual parsing error in the corrected interaction table (4.4), the ecological consequences of the chosen competition coefficients (4.5), the absence of a validation of the starting state (4.6), the reconstructed parts of the river network (4.7), the classification of *Alburnus alburnus* (4.8) and the reproducibility of the reported results (4.9).

# 2. Where our review was wrong

We are grateful for these corrections and accept them.

| **Item** | **Our statement** | **Correct reading** |
| --- | --- | --- |
| C1 | 89 invasive→native pairs with Σα = −65.5; iterating over all rows would recover *Lepomis* | Without the *Lepomis* row the figures are 79 pairs and −58.1. The CSV has an Lg row but no Lg column, so a column or a long-format table is needed. |
| C2 | Equation with c = −α | Our equation and the code we sketched were not consistent with each other: the code kept the shared carrying capacity, which implies c = 1 − α. See also 4.5. |
| C4 | An on-path tree of about 8,640 km | That was the sum of the parent edges only. With the 267 roots joined, the full tree is 13,871 km. |
| C7 | r = 0.81 | Your recomputation gives about 0.77 for the 2045 snapshot. |
| E10 | Settings existed only as undocumented environment variables | They were documented in docs/climate\_scenarios.md and the run metadata. The real gap was a committed configuration. |
| E20 | The flag in 06a is one-sided | It is two-sided: resid is an absolute value (06a\_qc\_waterbase\_direct.py, line 84). |
| Minor numbers | Native richness 2.52; 16 species share 8–30 °C; 200 m coordinate shift; 25 juvenile-only cases | 2.51; 15 species; about 234 m; 23 value-cells over 8 rows. |

On two items, the difference comes from the reference set rather than from an error: the 294 fishless sites in the review count all 776 rows of the environmental matrix (293 among the 775 model sites), and the 21-year burn-in is correct for one of the three matrices (21, 494 and 532 years in the September runs).

# 3. What we could confirm

We recomputed the following from the repository and found them exactly as stated in the response.

| **Topic** | **Confirmed** |
| --- | --- |
| Test suite | We ran all test files: **3,769 tests pass**. One testset errors only because the daily forcing file water\_temp\_daily\_guadex\_sites\_wide.csv is excluded from the repository (see 4.9). |
| C1 interaction table | 227 rows, 22 ambiguous cells and 44 ambiguous rows; 24 × 24 matrix including LG and AAL; correct direction in all cells we spot-checked (e.g. MS→SA −0.8, ST→SA −0.5); the corrected runs load the long table. |
| C2 growth form | Interaction term inside the growth bracket, clamp(1 − ΣN/K + interaction, −1, 2) (src/ode\_model.jl, lines 168–181). |
| C3 interim route | Observed start, exactly three baseline years, control required, scenario-minus-control deltas. |
| C4 network | 774 edges (508 on-path parents + 266 root joins), 13,871.5 km, 4 parent-elevation inversions, acyclic and deterministic. |
| C5/E1 temperature | tw\_baseline\_mean by exact column name; range 9.79–18.06 °C; basin mean 15.82 °C; mean 12.83 °C over the 36 field trout sites; sub-catchment 32.1 has two TEMP\_MEDIA\_SC values. |
| E5, E6, E8 | Offsets 1–20 (2045 = end of year); quasi-extinction over baseline-established sites only; metric renamed realised\_richness\_loss. |
| E12 obstacles | Status classes 1,274 / 185 / 199; 1,036 of 1,658 matched; passabilities multiply along links; IF not used. |
| E13–E16, E22 | Absence factor 1.0; 84 pools holding 56 % of density, capped at the non-pool median; 56 DIFÍCIL sites at K = 1e-6; AA/LR/MC with zero local growth and AA 9 → 5; bias correction +1.106 °C (Feb +2.224, Mar +2.287). |
| Viewer A1–A8 | All implemented and already deployed: cursor-centred zoom, panel-aware framing, arrow-key panning, working layer list, neutral water bodies, "only sites with data" by default, warning for unmatched keys, symmetric diverging scale, and per-species exports. |

# 4. Points we would kindly ask you to revisit

## 4.1 The control is shared, not matched by GCM (E2)

**What we found**

The response (E2 and section 2.5), the limitations document and the report state that every scenario is paired with a matched no-warming control for the same GCM. In the code there is a single control: the control scenario runs only gcm = "baseline" (run\_climate\_scenarios.jl, line 493), and scenario\_minus\_control\_table uses first(controls) for every run (src/climate\_diagnostics.jl, lines 238–240). This agrees with Table 1 of the report (44 runs plus one control). In addition, the control and the spin-up are still driven by the smoothed day-of-year climatology of the ensemble-median historical series (src/temperature\_forcing.jl, lines 251–273), while each scenario uses its own GCM’s daily series, so the original E2 concern remains.

**What we would suggest**

- Either correct the description (one shared control) or run one control per GCM driven by that GCM’s own historical series, which would make the scenario-minus-control differences like for like.

## 4.2 The slope we quoted drives the projections (E20, E21)

**What we found**

The response states that the coefficients we quoted (0.510 and −0.163) belong to the national model, whereas the Guadalquivir model gives 0.436 and −0.299. Both values are correct, but the projections that feed GuadeX are produced with the national model: 10\_project.py fits the "Stage-B (all-Spain)" model (lines 1–6 and 46), 13\_project\_guadex\_sites.py uses it (lines 368–376), and projection\_model\_coefficients.json stores ta = 0.5097. The slope we quoted is therefore the one that determines the forcing of the fish model. For the same reason, the 257 coded Waterbase values in the all-Spain calibration do reach the GuadeX forcing, even though the Guadalquivir-only model does not contain them.

**What we would suggest**

- Revise the wording of E20 and E21 in the response accordingly.
- Either refit the all-Spain model without the coded values, or document why the national rather than the basin model is used for the projections.

## 4.3 The 21, 494 and 532 burn-in years do not belong to the corrected runs (E4)

**What we found**

Under the corrected configuration the alternative-interaction sweep uses the same interim route as the climate runs, with three baseline years and the E4 stop rule bypassed (parameters\_climate\_scenarios\_corrected.toml, lines 250–257; sensitivity\_core.jl, lines 481–494). The 21, 494 and 532 years cited from results/alt\_interactions\_corrected/burnin\_summary.csv can therefore only come from the September runs. The new composition stop rule, although implemented and tested, is inactive in every corrected run.

**What we would suggest**

- Clarify the source of these figures and state that the E4 rule is not used by the corrected route.

## 4.4 A residual parsing error in the corrected interaction table (C1)

**What we found**

In 16 cells of the form «interfere through competition and Xx through predation», the reverse direction, which the text supports only as competition, receives the predation value −0.5 instead of −0.3. The undirected clause is cut as lower\[m.offset:end\] (src/data\_preparation.jl, lines 906–907), so it still contains «through predation», which then takes priority. Examples in data/BIOTIC/interaction\_matrix\_long.csv: GH→LG (line 82), AAL→OM (line 121) and CG→AM (line 163), all −0.5. In addition, the nine cells reading «interfere through competition and predation» receive −0.5 in both directions and are not flagged as ambiguous; this is a reasonable judgement, but it would be worth stating.

**What we would suggest**

- Restrict the undirected clause to the competition part, regenerate the long table and re-run the affected experiments.

## 4.5 The chosen competition coefficients have strong ecological consequences (C2)

**What we found**

With the shared carrying capacity, the implemented form is equivalent to a competitive Lotka–Volterra model with c = 1 − α. Since every non-zero α in the matrix is −0.3, −0.5 or −0.8, every interacting pair has c between 1.3 and 1.8, that is, interspecific competition stronger than intraspecific competition. In that regime, two-species Lotka–Volterra dynamics lead to competitive exclusion or to outcomes that depend on the starting state, never to stable local coexistence; and pairs with no interaction (α = 0) still compete with c = 1 through the shared K. Our own equation, by contrast, dropped the shared capacity and used c = −α (0.3–0.8), which allows coexistence. Both are valid Lotka–Volterra forms; the choice is ecological rather than mathematical, and it probably has a large influence on the results.

We also note that the response mentions a test of the two-species equilibrium K(1 − c₁₂)/(1 − c₁₂c₂₁), but we could not find it in test/test\_ode.jl (the C2 testsets, lines 380–482, test the single-species equilibrium, the per-capita identity, the scaling with r and the legacy switch). With all c \> 1, that interior equilibrium would in any case be unstable.

**What we would suggest**

- Discuss explicitly which competition structure is intended (shared K plus interactions, or species-specific competition only).
- Run both forms as a structural sensitivity, at least for the climate endpoints, and report how much the conclusions depend on the choice.

## 4.6 The starting state is not yet validated against the survey (C3)

**What we found**

The interim route is a sensible temporary measure, and the report is admirably clear that the results are uncalibrated. Still, there is no species-by-species comparison of the post-spin-up state with the 2006–2009 survey, and the no-warming control loses 35.5 % of its total biomass and 40.8 % of its native biomass between 2026 and 2045 (report §4.1). The report itself notes that the indicators depend much more on this transient and on network hydraulics than on the roughly 0.6 °C of realised warming. In a non-linear model the response to warming depends on the state, so subtracting the control does not fully remove the transient.

**What we would suggest**

- Report, for each species, observed versus modelled occupancy (and a skill score such as TSS) after the three-year spin-up, and the change of each species in the control between 2026 and 2045.
- Use this comparison as the acceptance criterion for the future calibration.

## 4.7 Parts of the river network are reconstructed rather than observed (C4, E12)

**What we found**

The on-path tree reproduces exactly, but 267 of the 775 sites (34 %) have no on-path parent and are joined by minimum-spanning-tree edges. These root joins make up 38 % of the total length (median 15.8 km, maximum 206.6 km) and are straight connections between different flow paths. Of the matched obstacles, 341 (245 of them applied) lie on these edges. Because passabilities now multiply, one link carries 33 obstacles and an upstream factor of 1e-33, and 28 links carry five or more. The root joins also add 41 elevation inversions to the 4 reported for the parent edges. Finally, the number of roots is sensitive to the on-path tolerance: 440, 267 and 119 roots for tolerances of (0.01, 100 m), (0.02, 200 m) and (0.05, 500 m).

**What we would suggest**

- Report these figures alongside the network description, and test the sensitivity of the obstacle results to the tolerance.
- Consider a floor for compounded passability, or assign obstacles only to on-path edges, until a flow-directed network is available.

## 4.8 *Alburnus alburnus* is still classified as migratory

**What we found**

The report and the response acknowledge that the bleak (AAL) is an exotic species misclassified as migratory, so invasive richness and biomass are lower bounds. The fix is very small: AAL appears in migratory in parameters\_climate\_scenarios\_corrected.toml (line 132) and parameters.toml (line 85), while scripts/dump\_interaction\_matrices.jl and test/test\_data\_preparation.jl already treat it as invasive. Because the species groups only affect the output metrics, not the dynamics, the invasive and migratory indicators could be recomputed from the saved states without re-integrating the model.

**What we would suggest**

- Move AAL to the invasive group in both configuration files, update test/test\_parameters.jl and recompute the invasive and migratory metrics.

## 4.9 The reported results cannot yet be reproduced from the repository

**What we found**

The results/ directory is excluded by .gitignore (line 32), and the 45 daily forcing files are excluded by guadex\_tw/.gitignore (line 15). As a result, none of the numbers in the October report can be checked or rebuilt from the repository (control transient, invasive richness, warming endpoints, trout exposure and quasi-extinction, factor rankings), the corrected runs cannot be re-executed, and the testset «model temperature equals corrected daily series» fails in a clean checkout. The report also states that the corrected exports can be browsed in the published viewer, but only the demonstration files are hosted there.

**What we would suggest**

- Deposit the corrected results/ and the daily forcing files in a citable archive (for example Zenodo or the project’s Koofr space), with checksums, and reference them from the report.
- Host at least one corrected run in the published viewer, or adjust the sentence in the report.

## 4.10 Smaller inconsistencies between documents

**What we found**

We noticed a few points that are easy to align:

**What we would suggest**

- The limitations document still mentions about 1 °C of central warming and a maximum of 18.37 °C; the report and the data give about 0.6 °C and 18.06 °C.
- The response states that the span-based factor ranking was replaced by within-design sums of squares, but the report still headlines the span-based mean ranks, including in the abstract.
- Under E18, 1,584 cells are said to be zeroed by the salinity envelope; as far as we can see, in the corrected run the envelope acts after the E16 zeroing and therefore changes only the 528 AB cells.
- The endpoints are quoted as 0.208 and 1.091 °C (2026–2045 mean), while the configuration selects them on the 2036–2045 window (0.180 and 1.292 °C); it would help to state the window.
- Alien-species introductions are listed as a deferred team decision in the response, but are not mentioned in the limitations of the report, although they cap invasive richness.
- The response and the correction-status document still say that all work is uncommitted, whereas it has been in the repository since 1 October.

# 5. Closing remarks

We would like to repeat that the October release is a major step forward and that most of the corrections are implemented exactly as described. The points above are offered in the same constructive spirit as the review: a few statements need to be aligned with the code, one parsing error remains, and two modelling questions, the competition structure and the validation of the starting state, seem to us the most important open issues before the projections can support management advice. We would be very happy to discuss any of them, to look at the corrected outputs once they are available, and to help with the biological inputs, in particular the physiological thermal data and the selection of alien species and introduction points.

Carlos
