> # ANNOTATION — REVIEW RESPONDED TO (October 2026)
>
> This review **has been responded to**. The reviewer's original text below is **unchanged and
> retained in full**; this banner is the only addition. See:
>
> - **`docs/GuadeX_Response_to_Reviewers_Oct2026.md`** — the formal point-by-point response.
> - **`docs/GuadeX_Review_Assessment_Sept2026.md`** — the independent adjudication of each claim.
> - **`docs/GuadeX_Correction_Status_Sept2026.md`** — what was implemented, with evidence.
>
> The model was corrected and re-run, and the report reissued as
> **`docs/FinalReportOctober2026.tex`** (superseding `legacy/FinalReportSeptember2026.tex`).
>
> **Status key** (as used in the response): **Accepted** — point valid, corrected; **Partly** —
> valid in substance but a stated figure/remedy is wrong or only partially implemented;
> **Not accepted** — the claim is inaccurate or not reproducible; **Deferred** — valid but
> deliberately out of scope for this iteration.
>
> Section status at a glance: **C1–C7** — Accepted in substance; C1/C2/C4/C7 partly (specific
> figures/remedies corrected); **E1–E24** — Accepted / Partly as detailed in the response;
> **minor issues** — mostly Accepted, with #2, #5, #7, #8, #11, #12 Partly; **report-statements
> table** — Accepted except "K = 1× observed", which is Not accepted as a report error;
> **"What stands"** — some positives marked Not verifiable; **Part II fixes** — adopted, partly
> adopted, or Deferred (notably the `IF/10` obstacle mapping, Not accepted).

---

**GuadeX: Technical Review and Proposed Corrections**

This document has two parts. **Part I** is the technical review: it
describes the issues we found, why they matter and the evidence for
each. **Part II** proposes a complete correction for every issue,
includes a simple check that the fix has worked, and, at the end,
gathers the few items that depend on the team's decisions, data, or
computing resources.

The outset is that **many of these issues are easy to fix**. Most are
local coding or configuration matters, and the critical ones that affect
interactions, site temperatures, and the warming axis can each be
corrected with a few lines of code.

**Part I. Technical review**

**Material reviewed.** GitHub repository melian009/GuadeX (branch
master, commit 9702169, 25 September 2026) including all Git LFS data;
the water-temperature methodology in guadex_tw/ (REPORT.md,
LIMITATIONS.md and scripts 01--13); the simulation outputs shared on
Koofr (climate_scenarios_k1x_burnin, sensitivity_obstacles,
report_viewer_exports); the September 2026 technical report and its
seven-page brief; and the source code of the viewer.

# Critical issues

Each of the following issues is, on its own, sufficient to change the
sign or the magnitude of the main conclusions. We would recommend
addressing all of them before any experiment is re-run.

## C1. The interaction matrix is loaded with inverted direction and only one-way

**What we found**

-   load_interaction_matrix (src/data_preparation.jl, lines 564--612)
    stores the text of cell (row A, column B) in α\[A,B\], whereas the
    model equation (src/ode_model.jl, lines 111--128) reads α\[s,j\] as
    the effect of species *j* on species *s*.

-   114 cells name the column species as the one affected (for example
    row Ms, column Sa: *"displaces Sa through predation"*), so they are
    loaded inverted: in the model, *S. alburnoides* depresses largemouth
    bass. Only 7 cells name the row species and are loaded correctly.

-   The CSV is filled only below the diagonal (253 cells against 0). The
    107 undirected relationships (54 *"No coexist"* and 53 competition
    entries) therefore enter in one direction only.

-   As a result, in the "original" matrix the effect of the ten invasive
    species on native species is **zero in all 100 pairs**, while
    natives depress invasives in 89 pairs (Σα = −65.5).

-   *"No coexist"*, which in the source table describes allopatry (for
    example brown trout and *Aphanius baeticus*), is coded as
    competitive exclusion (α = −1).

-   *Lepomis gibbosus* has a row but no column in the CSV, so the parser
    skips it and it receives no interactions at all. The unit test in
    test/test_data_preparation.jl (lines 187--188) asserts the inverted
    value.

**Why it matters**

The mechanism proposed in the report for the effect of barriers
(restricting competitively superior invasive species) cannot operate in
the model as currently written. Conversely, the "invasive-favouring"
matrix, presented as a deliberately extreme case, is the only one with
the ecologically correct direction of effects, although its magnitudes
remain extreme.

**Recommended solution**

-   **Recode the matrix as an explicit long-format table** with one row
    per directed interaction: source species, target species, mechanism
    (predation, competition, hybridisation, none), qualitative strength,
    confidence and supporting reference. This removes the ambiguity of
    free text, makes the matrix easy to review by fish biologists and
    lets the code read it without parsing. We would be glad to help
    review the pairwise entries.

-   **In the meantime, parse the affected species from the text** rather
    than transposing the matrix (see the snippet below), and treat
    undirected competition as symmetric.

-   **Treat "No coexist" as absence of interaction.** Whether two
    species co-occur is already governed by the thermal and habitat
    filters.

-   **Add the missing *Lepomis* column**, extend the parser to Spanish
    wording (for example *"depreda"*), and correct the unit test.

-   **Use the qualitative categories as prior ranges** rather than fixed
    values, and estimate the coefficients during calibration (C3).

> const TARGET_RX = r\"(?:affects\|displaces\|interfere\[s\]?(?:
> with)?)\\s+(\[A-Z\]\[a-z\]{1,2})\"
>
> \# cell in row A, column B; v = parse_interaction_string(txt)
>
> if occursin(\"no coexist\", lowercase(txt))
>
> \# allopatry: not an interaction
>
> elseif (m = match(TARGET_RX, txt)) === nothing
>
> α\[i(A), i(B)\] = v; α\[i(B), i(A)\] = v \# undirected competition:
> symmetric
>
> else
>
> tgt = m.captures\[1\]; src = tgt == B ? A : B
>
> α\[i(tgt), i(src)\] = v \# α\[s,j\] = effect of j on s
>
> end

**How to verify the fix**

Add unit tests asserting, for example, that α\[SA, MS\] \< 0 and α\[MS,
SA\] == 0 (largemouth bass preys on *S. alburnoides*, not the reverse),
that the summed effect of invasives on natives is negative, and that
every species has both a row and a column.

## C2. The interaction term is not scaled by the growth rate

**What we found**

-   The intrinsic growth rate is converted to a daily rate (r/365, line
    1088; between 0.0003 and 0.011 per day), but the term Σα·N/K enters
    the equation without being multiplied by *r* or by the environmental
    suitability *W*. With α as large as −1, interactions act on a
    per-day scale that is between one and three orders of magnitude
    stronger than growth.

-   Equation 1 of the technical report has the same form, so the issue
    is in the formulation as well as in the code.

**Why it matters**

Equilibria are set by the α coefficients rather than by r·K: any species
receiving a negative α from an abundant species is effectively excluded.
This explains both the 90% collapse of native biomass under the
invasive-favouring matrix and why the interaction matrix ranks first
among factors.

**Recommended solution**

-   **Reformulate the local dynamics as a competitive Lotka--Volterra
    model** in which interactions sit inside the growth term, with
    competition coefficients c derived from the qualitative scale:

dN~i,s~/dt = r~s~ W~i,s~(t) N~i,s~ \[1 − (Σ~j~ c~sj~ N~i,j~) / K~i~\] −
m^heat^~i,s~(t) N~i,s~ + Σ~j∈nb(i)~ (m~ji~ N~j,s~ − m~ij~ N~i,s~)

with c~ss~ = 1 and c~sj~ = −α~sj~ ≥ 0

-   **If predation is to be represented explicitly**, model it as a
    separate consumer--resource term (for example a Holling type II
    functional response with a conversion efficiency) rather than as a
    negative competition coefficient, so that predators benefit from
    prey.

-   **Keep all rates in consistent units** (per day) and document them
    in a parameter table.

**How to verify the fix**

For a single site with two species and no dispersal, the simulated
equilibrium should match the analytical Lotka--Volterra equilibrium of
species 1, K(1 − c₁₂)/(1 − c₁₂c₂₁). A unit test with this case will
catch any future scaling error.

## C3. The initial state after burn-in does not match the observed community

**What we found**

-   After a 373-year burn-in the model's 2026 starting state differs
    profoundly from the field survey it should represent (see the
    species-level table in Section 4, run control/baseline).

-   Of the ten invasive species only *Lepomis* persists, which is the
    only species with no interactions (C1). *Gambusia holbrooki*, the
    most abundant species in the survey (15,074 density units at 65
    sites), falls to zero, as do all four migratory species.

-   Mean native richness rises from 1.32 species per site (observed) to
    2.52 (model). Brown trout exceeds the presence threshold at 34
    sites, but only 13 of them coincide with the 36 sites where it was
    caught, and 58% of its modelled biomass lies at sites with no trout
    records.

-   The constant "invasive richness = 0.2168" reported for every
    scenario corresponds to *Lepomis* being present at 168 of the 775
    sites.

**Why it matters**

All projections start from a community generated by the model rather
than from the observed one, so no statement about invasive species, and
few about natives, can currently be supported by these outputs.

**Recommended solution**

-   **Calibrate before projecting.** Estimate r, K, the interaction
    coefficients and dispersal so that the model reproduces the observed
    community, using simulation-based inference. Approximate Bayesian
    computation with sequential Monte Carlo is well suited to
    deterministic dynamical models of this kind (Toni et al., 2009;
    Hartig et al., 2011), and it fits naturally with the ABC framework
    already outlined in the project's working paper.

-   **Use informative summary statistics:** per-species occupancy
    agreement measured with the true skill statistic (Allouche et al.,
    2006), rank correlation of abundances, and the distribution of site
    richness.

-   **Fix acceptance criteria before calibrating** (for example a
    minimum TSS for species with at least ten occurrences) and validate
    on held-out sub-catchments (spatial block cross-validation),
    followed by posterior predictive checks.

-   **As an interim measure**, start projections from the observed state
    with a short spin-up, and report every scenario relative to a
    matched control rather than as an absolute trajectory.

**How to verify the fix**

Publish, for every species, observed versus modelled occupancy and
abundance at the starting state, as in the table in Section 4. No
projection should be interpreted until the pre-registered criteria are
met.

## C4. The dispersal graph does not follow the river network

**What we found**

-   build_distance_matrix (lines 698--816) does not reconstruct the
    dendritic topology. Within each sub-catchment it chains sites in
    order of their distance to the Guadalquivir, and it then chains the
    77 sub-catchment outlets in order of elevation.

-   79% of the within-sub-catchment links jump between different
    tributaries, and in 34% of them the site labelled "upstream" is
    lower than its "downstream" neighbour.

-   Total link length is 31,215 km, against 8,783 km for the minimum
    spanning tree built from the same network distances. The outlet
    chain links, for example, marsh sites with the Guadiamar by
    elevation rather than by their position on the main stem.

-   Ties in outlet elevation are resolved by the iteration order of a
    Dict, so the number of matched obstacles varies between 966 and 986
    depending on that order.

**Why it matters**

The effect of upstream movement cost, the reported sign reversal of
passability effects and Figures 8--11 all rest on artificial links.
Results are also not strictly reproducible.

**Recommended solution**

-   **Preferred solution: build the directed graph from the river-line
    layer** (SW_Line_4C\_.shp). Snap each site to its river segment,
    derive the network topology and flow direction, and connect each
    site to its nearest downstream neighbour along the channel. Standard
    tools exist for this in R (for example sfnetworks or riverdist) and
    Python (networkx with geopandas).

-   **Fallback: reconstruct the tree from the network-distance matrix.**
    Site *k* lies downstream of site *i* on the same flow path if D(i,k)
    ≈ d(i) − d(k), where d is the distance to the Guadalquivir. A Python
    prototype of this rule produced a tree of 8,637 km with only 4
    elevation inconsistencies, against 238 in the current graph.

-   **Rebuild the dam and obstacle overlays on the corrected graph** and
    break ties deterministically.

-   **For later statistical analyses of network data**, spatial
    stream-network models that account for flow-connected and
    flow-unconnected dependence are a natural option (Ver Hoef &
    Peterson, 2010).

> function dendritic_parents(D, dG; rtol=0.02, atol=200.0) \# D and dG
> in metres
>
> n = size(D, 1); parent = zeros(Int, n)
>
> for i in 1:n, k in 1:n
>
> (k == i \|\| dG\[k\] \>= dG\[i\]) && continue
>
> on_path = abs(D\[i,k\] - (dG\[i\] - dG\[k\])) \<= rtol\*D\[i,k\] +
> atol
>
> if on_path && (parent\[i\] == 0 \|\| D\[i,k\] \< D\[i,parent\[i\]\])
>
> parent\[i\] = k
>
> end
>
> end
>
> return parent \# 0 = joins the main stem; order these roots along the
> main channel
>
> end

**How to verify the fix**

The graph should have n − 1 edges, every edge should satisfy the on-path
condition, elevation should not increase downstream except for a handful
of documented cases, and total link length should be close to the
minimum spanning tree.

## C5. The site temperature level is sub-catchment air temperature for 1961--1990

**What we found**

-   extract_site_temperatures (lines 922--956) takes the first column
    whose name contains "TEMP", which is TEMP_MEDIA_SC. The data README
    defines it as the *"temperatura media anual de la subcuenca (serie
    histórica 1961--1990)"*. It is almost certainly an air-temperature
    climatology and has a single value per sub-catchment (69 distinct
    values across 775 sites).

-   The model adds to this level the seasonal cycle and anomaly of water
    temperature relative to 1986--2005 (temperature_forcing.jl, line
    213), so the level and the anomaly also refer to different periods.

-   In sub-catchment 33 (Genil) all sites between 67 m and 1,997 m
    receive 15.17 °C. The 36 trout sites receive on average 14.35 °C,
    0.88 σ above the trout optimum; with the water-temperature level of
    guadex_tw they would receive 11.74 °C, slightly below it.

-   The "basin mean of 16.0 °C" in the report is the mean of this
    column. A function to set a water-temperature baseline
    (with_temperature_baseline, spin_up.jl, line 258) already exists but
    is not called by any run script.

**Why it matters**

The statement that brown trout lies "1.5 thermal standard deviations on
the warm side of its optimum", and with it much of the reported 9--11%
decline, appears to be an artefact of the temperature level. The model
also cannot represent thermal refugia within a sub-catchment.

**Recommended solution**

-   **Drive the model with the site-level daily absolute water
    temperature produced by guadex_tw**, after the bias and elevation
    corrections described in E22, and select the column by its exact
    name.

-   **Use one and the same temperature series** for the model equation,
    the calibration of the heat-stress constant and the exposure
    diagnostics (E1).

-   **Validate the site temperatures** against the 46 Guadalquivir
    monitoring sites before running experiments.

**How to verify the fix**

A test should fail if the temperature column is not the expected
water-temperature field, and the report should show modelled versus
observed water temperature at the monitoring sites.

## C6. Thermal niches are derived from distributional ranges rather than physiology

**What we found**

-   The thermal optimum is the midpoint of the TEMPERATURE_C range and σ
    is one sixth of that range. Sixteen of the 24 species share the same
    range, "8 to 30 °C" (optimum 19 °C, σ = 3.67 °C), and no site
    exceeds 18.37 °C.

-   *Gambusia*, *Lepomis*, *Cyprinus* and *Micropterus* are given an
    optimum of 19 °C, although their growth optima and thermal
    tolerances are considerably higher. For *S. trutta* the upper limit
    is set at 20 °C, whereas its upper incipient lethal temperature is
    close to 25 °C (Elliott & Elliott, 2010).

**Why it matters**

Warming improves conditions for these sixteen species at every site, and
the model cannot distinguish native from invasive species by their
thermal response.

**Recommended solution**

-   **Parameterise each species with physiological data**: optimum
    temperature for growth and critical thermal maximum (or upper
    incipient lethal temperature), compiled from the literature with
    references. For brown trout, Elliott and Elliott (2010) provide a
    thorough synthesis.

-   **Use an asymmetric thermal performance curve** that rises gradually
    to the optimum and falls steeply towards the upper limit, for
    example the Sharpe--Schoolfield formulation (Schoolfield et al.,
    1981), instead of a symmetric Gaussian.

-   **Propagate the uncertainty**: treat optima and limits as priors in
    the calibration and run a sensitivity analysis on the position of
    the optimum (for example ±2 °C). The current range midpoints can be
    kept as one sensitivity case.

**How to verify the fix**

A table of thermal parameters per species, with sources, should
accompany the model description.

## C7. The warming axis and the trout exposure metric do not reflect the applied forcing

**What we found**

-   warming_end_degc is taken from basin_warming_curve, which is
    computed from water_temp_future_2045.csv: six gauging stations and
    20-year window means. The simulations, however, were forced with
    daily per-GCM series at the 775 sites.

-   We reproduced the report's values exactly (0.737, 0.873, 0.880 and
    0.973 °C; range 0.284--1.365 °C). The warming actually applied by
    2045 differs (r = 0.81 across runs): for example ssp126/IITM-ESM,
    the "cool endpoint", is listed at 0.28 °C but received 0.70 °C, and
    SSP3-7.0 actually lies below SSP1-2.6.

-   The "131 → 139--174 days above 20 °C" is the maximum over all sites,
    reached at site 1.22.2, a lowland reach without trout, and it is
    computed on the guadex_tw series rather than on the temperature the
    model uses. At trout sites the control gives a mean of 5.2 days per
    year and a maximum of 65.

**Why it matters**

The −9.2% per °C slope, the total-biomass regression and the choice of
cool and warm endpoints for the other two experiments are fitted against
a variable that is not the forcing. The exposure mechanism attributed to
trout does not describe the sites where trout lives.

**Recommended solution**

-   **Record the realised forcing for every run**, for example the mean
    daily anomaly over 2026--2045 and over 2036--2045 across the 775
    sites, and use it as the regressor.

-   **Account for the ensemble structure** when fitting dose--response
    relationships, for example with a mixed model including a random
    intercept per GCM (Zuur et al., 2009), and report the spread across
    climate models.

-   **Re-select the cool and warm endpoints** on the realised forcing.

-   **Compute exposure at the sites occupied by each species**,
    reporting the mean or median rather than the basin maximum.

**How to verify the fix**

The run index should contain the realised forcing, and the exposure
figure should state the set of sites over which it is computed.

# Initial model state compared with the observed community

The table compares, for the 775 model sites, the field survey with the
2026 state after burn-in (run control/baseline). A site counts as
occupied in the model when density exceeds 0.1 units, the threshold used
in the report for richness. Densities are in the units of the data file
(500 × individuals caught / area sampled in m²), which are not
documented in the file itself.

  ---------------------------------------------------------------------------------------------
  **Species**           **Group**    **Observed   **Model   **Both**   **Observed   **Model
                                     sites**      sites**              density**    density**
  --------------------- ------------ ------------ --------- ---------- ------------ -----------
  *Luciobarbus          Native       351          576       289        15,532       8,567
  sclateri* (LS)                                                                    

  *Squalius pyrenaicus* Native       161          99        33         8,831        1,735
  (SP)                                                                              

  *Squalius             Native       159          449       136        13,758       15,879
  alburnoides* (SA)                                                                 

  *Cobitis paludica*    Native       122          431       109        1,909        10,208
  (CP)                                                                              

  *Pseudochondrostoma   Native       112          92        30         1,813        1,348
  willkommii* (PW)                                                                  

  *Iberochondrostoma    Native       70           236       65         3,291        15,223
  lemmingii* (IL)                                                                   

  *Salmo trutta* (ST)   Native       36           34        13         680          839

  *Anaecypris           Native       11           18        6          151          520
  hispanica* (AH)                                                                   

  *Aphanius baeticus*   Native       1            5         1          167          7,079
  (AB)                                                                              

  *Iberochondrostoma    Native       1            9         1          1            145
  oretanum* (IO)                                                                    

  *Gambusia holbrooki*  Invasive     65           0         0          15,074       0
  (GH)                                                                              

  *Lepomis gibbosus*    Invasive     45           168       30         3,260        10,351
  (LG)                                                                              

  *Cyprinus carpio*     Invasive     40           0         0          127          0
  (CC)                                                                              

  *Micropterus          Invasive     19           0         0          132          0
  salmoides* (MS)                                                                   

  *Carassius gibelio*   Invasive     14           0         0          138          0
  (CG)                                                                              

  *Oncorhynchus mykiss* Invasive     12           0         0          87           0
  (OM)                                                                              

  *Gobio lozanoi* (GL)  Invasive     4            0         0          160          0

  *Tinca tinca* (TT)    Invasive     4            0         0          6            0

  *Esox lucius* (EL)    Invasive     2            0         0          1            0

  *Ameiurus melas* (AM) Invasive     1            0         0          1            0

  *Alburnus alburnus*   Invasive     11           0         0          134          0
  (AAL)                 (coded as                                                   
                        migratory)                                                  

  *Anguilla anguilla*   Migratory    9            0         0          11           0
  (AA)                                                                              

  *Liza ramada* (LR)    Migratory    5            0         0          39           0

  *Mugil cephalus* (MC) Migratory    4            0         0          69           0
  ---------------------------------------------------------------------------------------------

*Aphanius baeticus*, recorded at a single site, accounts for about 10%
of modelled basin biomass.

# Major issues

The following issues bias the results appreciably, although none of them
on its own reverses the conclusions. They are grouped into model and
experimental design (E1--E10), input data (E11--E19) and the
water-temperature model (E20--E24).

  ---------------------------------------------------------------------------
  **ID**   **Issue and evidence**              **Recommended solution**
  -------- ----------------------------------- ------------------------------
  E1       **Three different temperature       Use a single water-temperature
           series in one simulation.** The     series for the model equation,
           model equation uses TEMP_MEDIA_SC   the calibration of *k* and the
           plus the anomaly, whereas the       exposure diagnostics (see C5).
           heat-stress constant *k* and the    
           exposure table use the absolute     
           guadex_tw series, which is on       
           average 1.3 °C colder (2.6 °C at    
           trout sites).                       

  E2       **The control is not a              Run a control for each GCM
           like-for-like counterfactual.**     with its own historical
           Control and burn-in use a smoothed  series, define a present-day
           daily ensemble median, which halves baseline (for example
           day-to-day variability, whereas     2016--2035), and report every
           scenarios use each GCM's daily      result as scenario minus the
           series. At 2026 each scenario also  control of the same GCM.
           jumps from the 1986--2005 climate   
           to its GCM, including 0.34--0.52 °C 
           of warming already realised.        

  E3       **The obstacle sweep starts from a  Run a separate burn-in for
           single equilibrium.** All 32 runs   each combination of upstream
           start from the burn-in at upstream  cost and passability, or pair
           cost 0.05 and baseline passability  every run with a matched
           (run_sensitivity_report.jl, lines   no-warming control.
           88--141), so the sums of squares of 
           0.96--0.97 for upstream cost        
           measure transient re-equilibration  
           rather than sensitivity.            

  E4       **The burn-in convergence criterion Use a species-by-site
           is weak.** It only requires total   criterion (for example the
           biomass to change by less than      95th percentile of \|Δ log
           5·10⁻⁵ per year, while composition  N\|) together with the
           can keep drifting. The 21-year      residual of the equations and
           burn-in of the interaction          a minimum number of years, and
           experiment is probably premature.   report both values.

  E5       **"2045" is the state on 1 January  Report the state at the end of
           2045.** Snapshots are taken at      each year, or the mean of the
           offset\*365 for offsets 0 to 19     twelve monthly states.
           (outputs.jl, lines 385, 403 and     
           503), so only 19 years are reported 
           and the 2026 row precedes any       
           forcing.                            

  E6       **The quasi-extinct fraction counts Divide only by the sites
           sites that were never occupied.**   occupied in the baseline
           With a baseline density of 0 the    state, species by species.
           threshold is 0.1 and the site is    
           flagged from year 3 (outputs.jl,    
           lines 267--277 and 566--567). For   
           trout this fraction is ≈ 1 −        
           occupancy, which explains the       
           reported "saturation" of 0.956.     

  E7       **The factor ranking is not a valid Report sums of squares within
           comparison.** Each effect is the    each design only. To compare
           span between factor levels, so it   factors, use a global
           depends on how many levels were     sensitivity design, such as
           chosen and how extreme they are,    Latin hypercube sampling with
           and effects are standardised with   Sobol indices over plausible
           standard deviations from different  parameter distributions
           run sets. The interaction matrix is (Saltelli et al., 2010).
           also confounded with the            
           thermal-breadth multiplier and      
           burn-in length.                     

  E8       **The model does not estimate       Add demographic and
           extinction probability**, which is  environmental stochasticity
           the stated aim of the project. It   with replicate simulations and
           is deterministic and has no         species-specific
           replicates; the metric exported as  quasi-extinction thresholds,
           realised_richness_loss is realised  in the spirit of population
           richness loss above a threshold.    viability analysis (Morris &
                                               Doak, 2002; Lande et al.,
                                               2003). Rename the current
                                               metric accordingly.

  E9       **Old results can be silently       Add a hash of all parameters
           reused.** The burn-in cache key and and the git commit to cache
           the run fingerprints do not include keys and fingerprints, and
           the model parameters or the code    require a match before reusing
           version, and                        any result.
           run_climate_scenarios.jl (line 360) 
           skips any run whose output file     
           already exists.                     

  E10      **parameters.toml does not          Version the exact TOML file
           reproduce the report.** It sets     used for each experiment,
           annual-mean forcing, no per-GCM     write the full resolved
           series, no burn-in, K × 10, no heat configuration to
           stress and no control; the settings run_metadata.json, and
           actually used existed only as       document the model following
           undocumented environment variables. the ODD protocol (Grimm et
                                               al., 2020).

  E11      **The habitat index is not a        Correct the label, or build a
           trophic index.** According to the   habitat index from relevant
           data README, IET is the             variables (conductivity,
           bank-stability index (*índice de    pollution indicators, land
           estabilidad del talud*); the code   use, channelisation) with
           and report describe it as a         robust scaling, and test its
           trophic-state index. One outlier    sensitivity.
           (27.67) sets the scale.             

  E12      **Obstacles are matched poorly and  Snap obstacles to the real
           given uniform passability.** They   river network, map each
           are assigned to straight lines      obstacle's own fields to a
           between sites of the incorrect      passability value, multiply
           graph (797 of 968 fall on           passabilities along a link,
           artificial links). All receive 0.1  exclude demolished and
           upstream and 0.5 downstream,        out-of-basin structures, and
           although the inventory provides a   de-duplicate against the
           passability index, the presence and legacy dam layer. Connectivity
           condition of fish passes and the    indices such as the dendritic
           status of each structure (99        connectivity index (Cote et
           matched structures are demolished   al., 2009) and established
           or abandoned). Multiple obstacles   prioritisation procedures
           on one link do not accumulate, and  (Kemp & O'Hanley, 2010) would
           the legacy dam layer already        then allow barrier scenarios
           restricts 52% of links.             to be ranked on a sound basis.

  E13      **The 2006--2009 presences fix      Use a single *r* per species
           species distributions.** Absences   and let thermal and habitat
           receive 10% of *r*                  suitability determine
           (data_preparation.jl, lines         occupancy; if the 0.1 factor
           1090--1106), and 39% of sites had   is kept, treat it as a
           no fish. This prevents invasive     sensitivity parameter.
           expansion by construction.          

  E14      **Carrying capacity is driven by    Exclude or cap pools when
           isolated pools.** The 84 pool sites estimating *K*, standardise
           hold 56% of total density, and a 30 densities by season, and
           m² pool contains 9,533 units of     document the density units.
           *Gambusia*. Sampling consisted of a 
           single visit, in different seasons. 

  E15      **Fishless sites are treated as     Set *K* ≈ 0 or treat these as
           habitat that can be colonised.** Of transit-only nodes where the
           the 294 fishless sites, 56 were     field team indicated so.
           flagged by the field team as unable 
           to hold fish (SIN_PECES = DIFÍCIL). 

  E16      **Species that do not reproduce in  Set *r* = 0 for these species
           the basin are modelled as local     and represent recruitment as
           populations.** Eel and both mullets immigration from the estuary;
           have REPRODU_WITHIN_THE_BASIN = no  remove the fish-farm records;
           but receive a local *r*. Four of    treat rainbow trout as a
           the nine eel records are fish-farm  stocked species.
           escapees in the Guadiato, according 
           to the data README.                 

  E17      **Growth and dispersal rates lack   Rebuild the values from
           reliable sources.** The cited       verifiable Iberian sources or
           sources include a study of mange in from trait-based dispersal
           kit foxes (for *Gambusia* growth)   models (Radinger & Wolter,
           and one of channel catfish (for its 2014), keep a single source of
           dispersal). Species coded as        values in parameters.toml, and
           strictly sedentary receive 6--25 km add *Alburnus*.
           per year, and *Alburnus* is missing 
           from both dictionaries and takes    
           default values.                     

  E18      **Available habitat filters are not Add elevation and salinity (or
           used.** Salinity (brackish for AB,  conductivity) envelopes to the
           LR and MC) and elevation ranges are suitability function.
           read but not applied, so brown      
           trout is eligible at 477 sites      
           below 500 m.                        

  E19      **Hydrology is absent.** The 261    Keep dry reaches as seasonal
           dry or effluent-only sites are      barriers or transit nodes, and
           removed, so fish cross these        scale *K* or connectivity with
           reaches at no cost. The CEDEX 2045  CEDEX runoff by hydrological
           files are loaded but not used,      unit.
           although they project runoff        
           declines of 17--55%.                

  E20      **Coded values contaminate the      Remove them, add a physical
           water-temperature calibration.**    consistency rule and refit;
           257 Waterbase observations are      the residual standard
           exactly 0.1667 or 0.3333 °C,        deviation falls from 2.65 to
           including summer days with air at   2.34 °C.
           25--34 °C. They pass quality        
           control, and the outlier flag is    
           never applied.                      

  E21      **The sensitivity of water to air   Calibrate on lagged or
           warming is probably                 distributed-lag air
           underestimated.** The slope is      temperature, which better
           estimated from same-day air         reflects the thermal inertia
           temperature (0.51); with 7- and     of rivers (Caissie, 2006), and
           30-day means it rises to 0.60 and   propagate the range of slopes
           0.64 and fits better, and the model to the projected warming.
           runs 0.5 °C too cold in recent      
           years. The size of the bias         
           (15--30%) is our inference.         

  E22      **Cold bias in the Guadalquivir is  Apply in script 13 the same
           not corrected in the series used by correction as in script 10,
           the fish model.** The bias in the   adjust the monthly effects to
           held-out basin is −1.1 °C (−2.2 °C  the basin, and report the bias
           in February--March), and the        by month.
           seasonal amplitude on the main stem 
           is exaggerated by 1--2 °C. Script   
           13 omits the bias and elevation     
           correction that script 10 applies.  

  E23      **Limited empirical support and     Complete the high-elevation
           extrapolation in elevation.** 64%   extraction, flag extrapolated
           of the Guadalquivir observations    sites and propagate the
           are from 2024 and two sites have a  uncertainty of the elevation
           single observation. The warming     term.
           figures rest on six stations        
           between 0.5 and 537 m, and 17 model 
           sites lie above the highest         
           calibration site in the basin       
           (1,127 m).                          

  E24      **The report misdescribes the       Correct the text; read the
           forcing.** Section 3.4 states that  coefficients from the
           all sites received the same         temperature model's JSON file
           anomaly; in fact it varies by site  and set the correct reference
           (0.82--1.08 °C in one run, r =      elevation.
           −0.79 with elevation). If elevation 
           scaling were switched on, it would  
           use a hard-coded reference          
           elevation of 0 m.                   
  ---------------------------------------------------------------------------

# Minor issues

None of these points changes the conclusions, but it would be convenient
to address them in the same revision.

-   **Rare species.** IO, AB and AM occur at a single site, EL at two,
    GL, TT and MC at four and LR at five; any extinction statistic for
    these species depends on one or a few sites.

-   **Juvenile-only presences.** 25 cases with measured fish but zero
    density are treated as absences, and IO has no juvenile column.

-   **Row-order alignment.** The initial state is aligned by row
    position rather than by site code; it works today only because the
    three files share the same order. Indexing by CODIGO or using
    innerjoin(\...; order=:left) would make it robust.

-   **Rate conversion.** r/365 assumes instantaneous rates; for finite
    annual rates log1p(r)/365 is appropriate. The *Gambusia* conversion
    described in the code comments is not the one applied.

-   **Calibration of *k*.** It is calibrated on the mean year, so the
    worst baseline year loses 5.9% rather than the stated 5%.

-   **Calendar.** The forcing has 7,305 days and the simulation 7,300,
    giving a seasonal drift of up to five days by 2045; the climatology
    drops day 366.

-   **Monthly steps.** The monthly fixed effects of the temperature
    model create steps of up to 1.55 °C between months in the daily
    series; a harmonic formulation is equivalent in AIC and smoother.

-   **Hot tail.** The linear air--water relationship overestimates by
    about 0.7 °C above 25 °C; a logistic form could be considered
    (Mohseni et al., 1998).

-   **Mixed-model details.** The likelihood-ratio test compares REML
    fits with different fixed effects (the conclusion holds, and
    strengthens, with ML; Zuur et al., 2009); the intraclass correlation
    is evaluated at 0 °C air temperature; and two in-sample models
    appear in the cross-validation table.

-   **Numerical details.** The logistic term is clipped to \[−1, 2\]
    without documentation; immigration uses unclipped densities while
    emigration uses clipped ones; and the positivity callback creates
    mass.

-   **Obstacle count and coordinates.** The number of matched obstacles
    (973 in the report) varies between 966 and 986 with tie order. We
    also suspect, without having verified it, that sites are in ED50
    while water-body layers are in ETRS89, which could imply a shift of
    about 200 m.

-   **Auxiliary data.** Ten headers of the environmental matrix are
    mistranslated (for example PERIMETRO_MOJADO contains precipitation),
    although these columns are not yet used; the CEDEX files contain a
    non-existent "SSP285" scenario; site 1.30.20 has no network
    distances and site 1.14.10 differs by 863 m between two coordinate
    sources.

-   **Code and tests.** SPECIES_CODES lists "AA" twice and omits "AAL"
    (it is unused); several tests are vacuous; and no test currently
    checks mass conservation, network topology, the sign of α or the
    temperature level.

-   **Viewer.** The all-metrics CSV mixes the 15 metrics and 20 years in
    a single slider; files aggregated by sub-catchment, water body or
    basin are not rendered; no species-level output is exported; and the
    README mentions 1,037 sites while the model uses 775.

# Statements in the reports that would benefit from revision

Independently of the results that will change once the model is
corrected, the technical report and the brief contain several statements
that already differ from the code or the data. We list them so that they
can be revised in the next version.

  -----------------------------------------------------------------------
  **Statement**          **Where**     **What the code or data show**
  ---------------------- ------------- ----------------------------------
  775 sites grouped into Report §2,    288 water bodies plus one
  774 water bodies       §3.1, Table 2 unassigned site (1.2.2), and 77
                                       sub-catchments; 774 is the number
                                       of network links

  *Alburnus alburnus*    Report, Table Its code is AAL and it is an
  coded "Aa" and classed 1             exotic species; "Aa" is *Anguilla*
  as migratory                         in the interaction file

  Daily absolute water   Abstract,     The level is the 1961--1990
  temperature at each    §3.4          sub-catchment mean
  site, with anomalies                 (TEMP_MEDIA_SC), not an observed
  added to observed site               water temperature (C5)
  temperatures                         

  All sites receive the  §3.4,         The daily anomaly varies by site
  same basin-scale       Limitation 8, and with elevation (E24)
  anomaly                brief         

  Habitat index derived  §3.2          IET is the bank-stability index
  from the trophic-state               (E11)
  index                                

  Basin mean of 16.0 °C; §4, Figure 1  16.0 °C is the mean of
  trout 1.5 σ on the                   TEMP_MEDIA_SC; with water
  warm side of its                     temperature, trout lies close to
  optimum                              its optimum (C5)

  Maximum days above 20  §4, Figure    This is the basin maximum, at a
  °C rising from 131 to  12, brief     lowland site without trout; at
  139--174               Figure 2      trout sites the control gives a
                                       mean of 5.2 days and a maximum of
                                       65 (C7)

  Basin-mean warming of  Abstract, §4, This is not the forcing applied in
  0.74--0.97 °C (range   Table 3       the simulations (C7)
  0.28--1.37)                          

  Invasive richness of   Figure 2      It is *Lepomis* at 168 sites; the
  0.2168 because no                    other nine invasive species are
  invasive population                  absent from the initial state (C3)
  crosses the threshold                

  No warming-driven      §4            Expansion is prevented by
  invasive expansion was               construction (C1, C3, C6, E13)
  detected                             

  Blocking helps by      §5.2, brief   In the original matrix invasive
  restricting                          species have no effect on natives
  competitively                        (C1)
  effective invasive                   
  taxa                                 

  The invasive-favouring §3.5, §5.3,   Its magnitudes are extreme, but it
  matrix is a            brief         is the only matrix with the
  deliberately extreme                 correct direction of effects (C1)
  case                                 

  3-D dendritic network  Abstract,     It is a chain by distance to the
  connecting adjacent    §3.3          Guadalquivir and by elevation, not
  sites                                the river tree (C4)

  Migratory taxa have    §3.3          Only MC and LR; AA and AAL have
  dispersal rates above                half the median
  five times the median                

  Final-year (2045)      Figures and   These are the states on 1 January
  metrics                tables        2045 (E5)

  Trout quasi-extinction §3.6,         This is an artefact of the
  is near saturation     Limitation 10 denominator (E6)

  Carrying capacity      §3.2          This matches what was run, but
  equal to 1 × observed                parameters.toml sets 10 × (E10)
  density                              

  1,037 sites; 12 native Repository    The model uses 775 sites and 10
  and 11 exotic species  README        native, 10 invasive and 4
                                       migratory species

  Observed state used as Code          The survey dates from 2006--2009
  the starting point     ("observed    and the starting state does not
                         2026          resemble it (C3); the survey date
                         snapshot")    should be stated
                         and report    
  -----------------------------------------------------------------------

The three messages of the brief depend on critical issues: the trout
message on C5, C6 and C7; the barrier-context message on C4 and E3; and
the recommendation to pair connectivity with invasive control on C1, C2
and C3. We would kindly suggest not circulating the brief until the
corrected model has been re-run.

# What stands

The statistics of the water-temperature model, the quality of the
distance data and the numerical machinery of the code are sound; the
difficulties lie in how these components are connected. Everything below
was recomputed from the repository files or the published outputs.

-   **Water-temperature model.** The cross-validation NSE (0.785 against
    0.773 for Spain; 0.762 against 0.729 for the Guadalquivir), the
    held-out-basin NSE (0.704, RMSE 2.92 °C), the coefficients (b1 =
    0.510; b3 = −0.163) and the 2036--2055 warming (0.71, 0.80, 0.83 and
    0.94 °C) are all reproduced. The elevation × air-temperature
    interaction is robust and becomes stronger when the test is repeated
    with ML. There are no duplicated rows or unit errors, each GCM is
    compared with its own historical run, and all use a standard
    calendar.

-   **Data.** The network-distance matrix is complete, symmetric and in
    metres. The density columns are correctly labelled (AA is eel and
    AAL is bleak). No species record is lost when moving from 1,037 to
    775 sites, because the excluded sites were dry or received only
    effluent. The distributions of the rare species are consistent with
    their known ranges: *Anaecypris* in the Bembézar, *I. oretanum* in a
    Jándula tributary, *Aphanius* near the river mouth and trout between
    524 and 1,997 m.

-   **Code.** The dispersal matrix conserves mass; the sum-of-squares
    formulas are correct for the balanced designs; and the numerical
    settings (Tsit5, tolerances of 10⁻⁶), the median dispersal rate of
    0.0274 km per day, the mean *K* of 86.6 and the richness threshold
    all match the report.

-   **Outputs.** The slopes and R² values in Table 3 are reproduced
    exactly from runs_index.csv, although they are computed on the wrong
    axis (C7).

-   **Hypotheses worth re-testing.** That brown trout acts as a thermal
    sentinel, that the effect of a barrier depends on its network
    context, and that aggregate indicators can hide species-level
    declines are all ecologically sensible ideas. The current model does
    not demonstrate them, but the corrected model may well confirm them.

**Part II. Proposed corrections**

*A complete solution for every issue, from the easiest to those that
depend on the most complex.*

This second part sets out, for every issue identified in Part I. We
would like to stress one point above all: **most of the issues are easy
to fix.** They are local coding or configuration matters, many of them a
few lines long, and the knowledge needed to correct them is already in
the repository or in the data. Only a handful require biological
decisions by the team, data that are not in the repository, or
substantial computing time, and these are gathered at the end of this
part.

# Easy corrections

These corrections are local and can be made and tested quickly. We give
the complete solution for each; where code is shown, it is the essential
part of the change. Code excerpts in this part are written in Julia, the
language of the model. They are simplified sketches that show the logic
of each fix and indicate where it belongs (function and file); they are
not drop-in replacements for the original code.

## Read the direction of each interaction from its text

**Proposed solution**

-   In load_interaction_matrix (src/data_preparation.jl), parse each
    cell into directed effects: the species named in the text is the one
    affected (target) and the other species of the pair acts on it
    (source). Store the result as α\[target, source\], which is what the
    model equation expects.

-   Treat undirected competition (*"interfere through competition"*) as
    symmetric, and read the Spanish clause *"A depreda B"* as predation
    of A on B.

-   Treat *"No coexist"* as absence of interaction (allopatry); the
    thermal and habitat filters already decide whether two species can
    co-occur.

-   Iterate over all rows of the CSV, so that *Lepomis gibbosus* (row
    without a column) is no longer skipped.

-   Flag the few cells whose direction cannot be read (for example
    *"interfere through competition and predation"*) for the team to
    decide (Section II.6).

-   Export the result as an explicit long-format table (source, target,
    mechanism, α, ambiguous) that the loader can read directly once
    validated.

> const TARGET_RX =
> r\"(?:affects\|displaces\|interfere\[s\]?)\\s+(\[A-Z\]\[a-z\]{1,2})\\s+through\"
>
> \# cell in row A, column B; v = strength of the stated mechanism
>
> if occursin(\"no coexist\", lowercase(txt)) \# allopatry: no effect
>
> elseif (m = match(TARGET_RX, txt)) === nothing
>
> α\[i(A), i(B)\] = α\[i(B), i(A)\] = v \# undirected competition
>
> else
>
> tgt = m.captures\[1\]; src = tgt == B ? A : B
>
> α\[i(tgt), i(src)\] = v \# effect of src on tgt
>
> end

**How to check it**

Unit tests on known pairs: largemouth bass on *S. alburnoides* gives
α\[SA, MS\] = −0.8 and α\[MS, SA\] = 0; the summed effect of invasive on
native species is negative; the test that currently asserts the inverted
value is updated.

## Place the interaction term inside the growth term

**Proposed solution**

-   In \_metacommunity_ode! (src/ode_model.jl) replace N \* (r_eff \*
    logistic + interaction − heat) by N \* (r_eff \* (logistic +
    interaction) − heat). This is the competitive Lotka--Volterra form,
    with competition coefficients c = 1 − α on top of the shared
    carrying capacity.

-   Keep the previous form available as an option, for comparison only.

-   If predation is to be represented explicitly later, add it as a
    separate consumer--resource term (for example a Holling type II
    response) so that predators benefit from their prey.

> \# before: dU\[i,s\] = N_is \* (r_eff \* logistic_term +
> interaction_term - heat_stress)
>
> dU\[i, s\] = N_is \* (r_eff \* clamp(1 - total_i / K_i +
> interaction_term, -1, 2) - heat_stress)

**How to check it**

For one site and two species the per-capita growth rate equals r(1 −
(N₁ + N₂)/K + α₁₂N₂/K) and doubles when r doubles, which was not the
case before.

## Use water temperature at each site, and one series throughout

**Proposed solution**

-   In extract_site_temperatures, read the site-level water temperature
    already produced by guadex_tw
    (guadex_tw/outputs/tables/water_temp_baseline_guadex_sites.csv,
    column tw_baseline_mean, 1986--2005), selecting the column by its
    exact name rather than by the substring "TEMP".

-   With this level, the temperature seen by the model equation equals
    the daily guadex_tw series, so the heat-stress calibration and the
    exposure diagnostics automatically use the same temperatures as the
    model (E1).

-   Keep TEMP_MEDIA_SC only as an optional legacy source, and fail
    loudly if the expected water-temperature file or column is missing.

**How to check it**

Mean baseline temperature at the 36 trout sites becomes 11.7 °C instead
of 14.4 °C; a test fails if the wrong column is read; site temperatures
are compared with the 46 monitoring sites.

## Record the forcing actually applied, and compute exposure where species live

**Proposed solution**

-   In run_climate_scenarios.jl, store for every run the mean anomaly
    actually applied to the 775 sites (for example over 2036--2045) and
    use it as the regressor of the dose--response analyses; keep
    warming_end_degc only as a labelled proxy.

-   Fit the dose--response with a random intercept per climate model and
    report the spread among models.

-   In the export, flag the sites where each species is established in
    the baseline state and summarise exposure (mean, median, maximum
    days above the thermal limit) over those sites only.

-   Re-select the cool and warm endpoints on the realised forcing.

**How to check it**

The run index contains the realised forcing; every exposure figure
states the set of sites on which it is based.

## A convergence criterion that also covers community composition

**Proposed solution**

-   Add to spin_up a criterion based on the 95th percentile, over all
    site × species cells, of the annual change in log density, together
    with a minimum number of years, and record the value reached in the
    run metadata.

**How to check it**

Burn-ins stop only when both total biomass and composition have settled;
the 21-year burn-in of the interaction experiment is repeated with this
criterion.

## Report end-of-year states

**Proposed solution**

-   Set the reporting offsets to 1, ..., 20 instead of 0, ..., 19 in
    run_climate_scenarios.jl, run_climate_experiments.jl and
    sensitivity_core.jl, so that the row labelled 2045 is the state
    after the 2045 forcing and aligns with the warming of that year.

**How to check it**

The first reported row is no longer identical to the initial state.

## Quasi-extinction only where a species was established

**Proposed solution**

-   In quasi_extinction_flags (src/outputs.jl), never flag a site where
    the baseline density of the species is at or below the presence
    threshold.

-   In quasi_extinction_summary and in the site table, divide by the
    number of sites (or species) established in the baseline, not by all
    of them.

**How to check it**

The trout quasi-extinct fraction no longer equals 1 − occupancy, and the
metric can re-enter the factor analysis.

## Reproducible runs and caches

**Proposed solution**

-   Compute a stable digest of every model parameter (interaction
    matrix, growth and dispersal rates, temperatures, capacities,
    thermal parameters, structural options) and add it, with the git
    commit, to burn-in cache keys and run fingerprints.

-   In run_climate_scenarios.jl, reuse an existing run only if its
    stored digest matches; otherwise re-run it.

-   Add a model-options section to parameters.toml, write all resolved
    settings to run_metadata.json, and commit a configuration file that
    reproduces the September runs (reconstructed from their run
    metadata).

**How to check it**

Changing any parameter changes the digest and triggers a new burn-in; a
reference run is reproduced exactly from the committed configuration.

## Correct the label of the habitat index

**Proposed solution**

-   Rename IET as the bank-stability index in the code comments, the
    documentation and the report. If a water-quality index is preferred,
    build it from conductivity, pollution indicators and land use with
    robust scaling, as a separate option.

**How to check it**

Sensitivity of the results to the habitat multiplier is reported.

## Explicit options for biological assumptions

**Proposed solution**

-   E13: make the 10% growth factor at unsurveyed sites an explicit
    parameter. For projecting the establishment and spread of alien
    species under warming we would recommend 1.0, letting temperature
    and habitat decide where each species can grow.

-   E15: optionally set the capacity to (near) zero at sites the field
    team judged unable to hold fish (SIN_PECES = DIFÍCIL).

-   E16: optionally set local growth to zero for species that do not
    reproduce in the basin, and remove the four fish-farm eel records.

-   E18: add a salinity envelope for brackish species. We would not
    recommend an elevation envelope for climate projections, because it
    would prevent upslope range shifts.

**How to check it**

Each option has a unit test, and its value is recorded in the run
metadata. The defaults are left to the team (Section II.6).

## Water-temperature calibration and basin bias

**Proposed solution**

-   E20: in 08_models.py, remove values of exactly 1/6 or 1/3 °C and
    readings below 2 °C on days with air above 12 °C, and apply the
    outlier flag that is currently computed but never used (the flag in
    06a should also be two-sided, as in 02).

-   E22: in 13_project_guadex_sites.py, add to every projected value the
    mean held-out-basin residual for its month. From the saved
    validation this correction is +1.1 °C on average and about +2.2 °C
    in February--March. It is identical for historical and scenario
    series, so projected anomalies are unchanged, while the level and
    the seasonal shape are corrected.

-   Repeat the mixed-model likelihood-ratio test with ML fits
    (reml=False).

**How to check it**

Recomputed validation metrics are reported alongside the current ones;
the basin bias by month is reported.

## Correct the text of the reports

**Proposed solution**

-   Revise the statements listed in Section 7 of Part I.

-   The figure of "774 water bodies" has a simple origin: the run
    metadata field crosswalk_water_bodies counts the *sites* with an
    assigned water body (775 − 1), not the water bodies. Renaming the
    field and adding the true count (288) removes the confusion.

**How to check it**

Each corrected statement is checked against the outputs.

## Minor. Minor issues

**Proposed solution**

-   Align initial states by site code (initial_state_from_density) and
    use innerjoin(...; order=:left).

-   Offer log1p(r)/365 for finite annual rates; handle the calendar
    consistently (365-day years or explicit leap days); document the
    clipping of the logistic term.

-   Replace the monthly steps of the water-temperature model by a
    harmonic seasonal term in the daily products.

-   Add the missing tests (mass conservation, network topology, sign of
    α, temperature level) and remove vacuous ones.

-   Viewer: export one time series per metric and per species, so that
    the viewer can show species-level projections.

**How to check it**

Covered by the new unit tests.

# Moderate corrections

These corrections require a new routine or a rebuilt input, but can be
developed and tested with the data already in the repository.

## Rebuild the river network

**Proposed solution**

-   Preferred: derive the directed graph from the river-line layer
    (SW_Line_4C\_.shp) by snapping sites to their segments and following
    flow direction.

-   Immediate alternative, using only the network-distance matrix: link
    each site to the nearest site lying downstream on the same flow
    path, i.e. with D(i,k) ≈ d(i) − d(k), where d is the network
    distance to the Guadalquivir; then join the remaining sites with a
    minimum spanning tree on network distance. A prototype of this rule
    gives a tree of about 8,640 km with only 4 elevation
    inconsistencies, against 31,215 km and 238 inconsistencies in the
    current graph.

-   Rebuild the dam and obstacle overlays on the new graph, and resolve
    ties deterministically.

> for i in 1:n, k in 1:n \# D, dG in metres
>
> (k == i \|\| dG\[k\] \>= dG\[i\]) && continue
>
> on_path = abs(D\[i,k\] - (dG\[i\] - dG\[k\])) \<= 0.02\*D\[i,k\] + 200
>
> on_path && (parent\[i\] == 0 \|\| D\[i,k\] \< D\[i,parent\[i\]\]) &&
> (parent\[i\] = k)
>
> end \# then an MST joins the roots

**How to check it**

The graph has n − 1 links, every link lies on a flow path, and total
length is close to the minimum spanning tree.

## Obstacles from the inventory

**Proposed solution**

-   Assign each obstacle to the corrected network (ideally by snapping
    it to the river lines).

-   Derive its upstream passability from the inventory index IF. We
    assume a 0--10 scale in which higher values are more passable; this
    interpretation should be confirmed.

-   Multiply passabilities along a link so that several obstacles
    accumulate, exclude demolished structures and those outside the
    basin, and de-duplicate against the legacy dam layer.

-   Barrier scenarios can then be ranked with established connectivity
    indices (Cote et al., 2009; Kemp & O'Hanley, 2010).

**How to check it**

Checks on matched, excluded and duplicated structures, and on the number
of obstacles per link.

## Like-for-like controls and burn-ins

**Proposed solution**

-   E2: add a no-warming control for each GCM, driven by that GCM's own
    1986--2005 daily series as anomalies from its own mean, and report
    every scenario as scenario minus the control of the same GCM. A
    present-day baseline window would additionally avoid the step of
    already-realised warming at 2026.

-   E3: compute a burn-in for each combination of upstream cost and
    passability, or pair each run with a matched control.

**How to check it**

Controls show no trend, and differences between scenario and control are
reported per GCM.

## Global sensitivity analysis

**Proposed solution**

-   Replace the factor ranking by within-design sums of squares plus a
    global design (Latin hypercube sampling with Sobol indices; Saltelli
    et al., 2010) over plausible parameter ranges agreed with the team.

**How to check it**

Sobol indices with bootstrap confidence intervals.

## Inputs that need rebuilding

**Proposed solution**

-   E14: exclude or cap isolated pools when estimating carrying
    capacity, and standardise densities by season.

-   E17: move all growth and dispersal rates to a single table
    (data/species_parameters.csv) that lists every species explicitly,
    including *Alburnus*, with its source and status. The values
    themselves are a team decision (Section II.6).

-   E19: keep dry reaches as seasonal barriers or transit nodes, and
    scale capacity or connectivity with the CEDEX runoff projections by
    hydrological unit.

-   E21: refit the air--water model with lagged or distributed-lag air
    temperature (Caissie, 2006) and propagate the range of slopes to the
    projected warming.

**How to check it**

Before-and-after comparisons of the affected outputs.

# Items that depend on the team's decisions

  ------------------------------------------------------------------------
  **Item**        **What we would prepare**       **What the team
                                                  decides**
  --------------- ------------------------------- ------------------------
  Interaction     The long-format draft generated Validation of each
  table (C1)      from the existing matrix, with  entry, particularly who
                  ambiguous entries flagged       preys on whom where the
                                                  text does not say, with
                                                  references

  Thermal niches  A table per species with        The values adopted,
  (C6)            optimum growth temperature and  especially for the alien
                  critical thermal maximum, with  species, whose response
                  sources, and an asymmetric      to warming is central to
                  thermal performance curve       the project
                  (Schoolfield et al., 1981;      
                  Elliott & Elliott, 2010 for     
                  brown trout)                    

  Growth and      Drafts from Iberian literature  Final values and
  dispersal rates and trait-based dispersal       plausible ranges
  (E17)           models (Radinger & Wolter,      
                  2014)                           

  Calibration     Code for the validation         The acceptance criteria,
  criteria (C3)   statistics (per-species TSS,    fixed before calibrating
                  rank correlation of abundance,  
                  richness) and a proposal of     
                  acceptance thresholds (Allouche 
                  et al., 2006)                   

  Biological      The options described in        The default for each
  options         Section II.3                    option
  (E13--E16)                                      

  Alien species   A mechanism to add species      Which species to include
                  absent from the 2006--2009      (for example species
                  survey and introduction events  recorded in the basin
                  (stocking, translocations) at   after the survey) and
                  given sites and years           where they are
                                                  introduced
  ------------------------------------------------------------------------

# Items requiring data not in the repository or computing resources 

  -----------------------------------------------------------------------
  **Item**            **Reason**                   **What would make it
                                                   possible**
  ------------------- ---------------------------- ----------------------
  Re-running the      The daily forcing files for  Access to those files
  44-run climate      each climate model are       
  ensemble            excluded from the repository 
                      and are not among the shared 
                      outputs                      

  Full calibration    They require thousands of    Running the prepared
  (C3) and stochastic simulations                  code on the team's
  replicates (E8)                                  computing cluster

  High-elevation      The regional climate         The team runs the
  water-temperature   projections are not in the   extraction script, or
  projections (E23)   repository and their         shares its outputs
                      download is probably not     
                      reachable from our           
                      environment                  

  Reproducing the     Some settings were supplied  Confirmation of the
  September runs      as environment variables     reconstructed
  exactly             that were not recorded       configuration

  Integrating the     This is the maintainers'     Review and integration
  changes into the    decision                     by the maintainers, if
  repository                                       they agree
  -----------------------------------------------------------------------
