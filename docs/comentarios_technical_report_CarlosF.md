*Thermal Change, Network Fragmentation and Biotic Structure in a
Mediterranean Fish Metacommunity. Final integrated report for the
Guadalquivir River Basin, 2026--2045 simulation horizon* (corrected
October 2026 release, 40 pp.). Code checked against the public GuadeX
repository (commit b8ee15a, 3 November 2026).

## 1. References

  -----------------------------------------------------------------------------------
  «Altermatt, F. (2013). Dendritic  The correct reference is: Altermatt, F. (2013).
  network structure and regional    Diversity in riverine metacommunities: a network
  dynamics of metapopulations and   perspective. *Aquatic Ecology*, 47, 365--377. DOI
  metacommunities. Critical Reviews 10.1007/s10452-013-9450-3
  in Environmental Science and      
  Technology, 43(16), 1734--1775.»  
  --------------------------------- -------------------------------------------------
  «Perkin, J. S., & Gido, K. B.     The correct reference probably is: Perkin, J.S.
  (2012). Fragmentation and         and Gido, K.B. (2012), Fragmentation alters
  defragmentation of dendritic      stream fish community structure in dendritic
  ecological networks. Freshwater   ecological networks. Ecological Applications, 22:
  Science, 31(1), 36--44.»          2176-2187. <https://doi.org/10.1890/12-0318.1>.

  «Tonkin, J. D., et al. (2018).    The correct reference is Tonkin, J. D., Poff, N.
  Preparing river ecosystems for an L., Bond, N. R., Horne, A., Merritt, D. M.,
  uncertain future. Nature,         Reynolds, L. V., Olden, J. D., Ruhi, A., & Lytle,
  556(7702), 427--433.»             D. A. (2019). Prepare river ecosystems for an
                                    uncertain future. *Nature*, *570*(7761),
                                    301--303.
                                    <https://doi.org/10.1038/d41586-019-01877-1>

  *"Herrera, M., Garcia-Berthou,    **The reference does not exist**
  E., Minano, P. S.,                
  Garcia-Meseguer, A. J., and       
  Fernandez-Delgado, C. (2022).     
  Long-term changes in the fish     
  community of a mediterranean      
  river basin: the role of dams and 
  potamodromous species. Water,     
  14(15):2345."*                    
  -----------------------------------------------------------------------------------

**\**

**R2. Equation 1**

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **§3.2 "Local population dynamics", Eq. 1 and the | **Issue.** Circle 1 shows the term  |
| sentence that follows it, p. 6**                  | +Σ c_sj N_i,j / K_i, and circle 2   |
|                                                   | defines c_sj = 1 − α_sj. Combined   |
| ![](media/image1.png){width="4.803149606299213in" | with −Σ N_i,j / K_i, the bracket    |
| height="0.9408234908136482in"}                    | reduces to 1 − Σ α_sj N_i,j / K_i:  |
|                                                   | the density-dependent term cancels  |
|                                                   | (α_ss = 0), and a negative α        |
|                                                   | (competition or predation, as       |
|                                                   | defined in the same sentence) would |
|                                                   | increase growth. This is the        |
|                                                   | opposite of the intended effect and |
|                                                   | of the implementation, 1.0 −        |
|                                                   | total_biomass_i / k_i +             |
|                                                   | interaction_term with               |
|                                                   | interaction_term = Σ α_sj N_i,j /   |
|                                                   | K_i (src/ode_model.jl, lines        |
|                                                   | 192--202 and 226).                  |
|                                                   |                                     |
|                                                   | **Correction.** Rewrite the bracket |
|                                                   | as 1 − Σ_j N_i,j / K_i + Σ_j α_sj   |
|                                                   | N_i,j / K_i, as in Eq. 7, and state |
|                                                   | once that this is equivalent to 1 − |
|                                                   | Σ_j c_sj N_i,j / K_i with c_sj = 1  |
|                                                   | − α_sj.                             |
+---------------------------------------------------+-------------------------------------+
| **§4.7 "Phase 2 re-run and structural             | **Issue.** Eq. 7, which is correct  |
| sensitivities", paragraph 3 ("Competition         | and matches the code, uses +Σ α_sj  |
| structure"), Eq. 7, p. 14**                       | N_i,j / K_i, whereas Eq. 1 uses +Σ  |
|                                                   | c_sj N_i,j / K_i, although both are |
| *"The primary specification keeps a shared        | presented as the primary            |
| per-site carrying capacity, with csj = 1 − αsj    | specification.                      |
| :"*                                               |                                     |
|                                                   | **Correction.** Make Eq. 1          |
|                                                   | identical to Eq. 7.                 |
+---------------------------------------------------+-------------------------------------+

**R3. Consistency between Methods and code**

+--------------------------------+-------------------------------------+
| **Technical report --- exact   | **Issue and how to correct it**     |
| text and position**            |                                     |
+================================+=====================================+
| **§3.5 "Experimental design",  | **Issue.** In the code, "No         |
| paragraph following Table 1,   | coexist" is treated as allopatry    |
| sentence 1, p. 8**             | and mapped to α = 0, not to −1.0    |
|                                | (src/data_preparation.jl, lines     |
| *"The original interaction     | 874, 1078 and 1104--1108).          |
| matrix maps qualitative        |                                     |
| categories to coefficients     | **Correction.** Replace "−1.0 no    |
| (−1.0 no coexistence, −0.8     | coexistence" with "no coexistence   |
| displacement, −0.5 predation,  | treated as allopatry (α = 0)".      |
| −0.3 competition or            |                                     |
| interference, −0.2 weak        |                                     |
| effect, 0.0 neutral), now      |                                     |
| parsed in the correct          |                                     |
| direction and with undirected  |                                     |
| competition made symmetric     |                                     |
| (C1)."*                        |                                     |
+--------------------------------+-------------------------------------+
| --- Not described in the       | **Issue.** In the corrected form,   |
| report (§3.2 is the            | the code truncates the growth       |
| appropriate place).            | bracket to the interval \[−1, 2\]:  |
|                                | clamp(growth_bracket, −1.0, 2.0)    |
|                                | (src/ode_model.jl, line 228). This  |
|                                | bounds both net growth and          |
|                                | over-capacity mortality, and it     |
|                                | shapes the dynamics at sites where  |
|                                | density exceeds K (for example the  |
|                                | near-zero-capacity DIFICIL sites).  |
|                                |                                     |
|                                | **Correction.** Add one sentence in |
|                                | §3.2 describing the truncation and  |
|                                | its rationale.                      |
+--------------------------------+-------------------------------------+

**\**

**R4. Brown-trout thermal exposure**

+--------------------------------+-------------------------------------+
| **Technical report --- exact   | **Issue and how to correct it**     |
| text and position**            |                                     |
+================================+=====================================+
| **Abstract, paragraph 3,       | **Issue.** The rise from 23 to 45   |
| sentence 2, p. 3**             | days reads as an effect of the      |
|                                | simulated warming, but most of it   |
| *"Exposure above the 20 °C     | is already present in the first     |
| upper limit at the 72          | scenario year: 42.7 days in 2026,   |
| baseline-established trout     | which §4.3 (paragraph 2, sentence   |
| sites had a                    | 4, p. 11) attributes to the offset  |
| baseline-climatology maximum   | between the projection products and |
| of 90 days per year, while the | the 1986--2005 baseline. Within     |
| 2045 scenario maxima were      | 2026--2045 the mean rises only from |
| 80--135 days (mean of the      | 42.7 to 45.4 days.                  |
| per-run maxima 109 days) and   |                                     |
| the mean exposure rose from 23 | **Correction.** Report the change   |
| to 45 days."*                  | within the simulated horizon (42.7  |
|                                | → 45.4 days) and attribute the step |
|                                | from 23 days to the offset between  |
|                                | climate products.                   |
+--------------------------------+-------------------------------------+
| **§4.3 "Brown trout as a       | **Issue.** As in the Abstract, the  |
| thermal--hydraulic sentinel",  | comparison with the baseline        |
| paragraph 2, sentences 1 and   | climatology mixes the product       |
| 3, p. 11**                     | offset with the simulated warming;  |
|                                | the next sentence of the same       |
| *"Thermal exposure increased   | paragraph makes this clear, but the |
| with warming. \[...\] The      | figures are not separated.          |
| per-run mean exposure rose     |                                     |
| from 23.4 days in the baseline | **Correction.** Add the 2026 → 2045 |
| climatology to a scenario mean | change for the mean (and for the    |
| of 45.4 days (range            | median, if available), so that the  |
| 30.2--65.7), and the median    | warming-attributable part is        |
| from 17.0 to 48.1 days (range  | explicit.                           |
| 24.5--81.5)."*                 |                                     |
+--------------------------------+-------------------------------------+

**R5. Interpretation of the passability effect**

+--------------------------------+-------------------------------------+
| **Technical report --- exact   | **Issue and how to correct it**     |
| text and position**            |                                     |
+================================+=====================================+
| **§4.4 "Hydraulic controls and | **Issue.** The report does not      |
| the passability effect",       | explain why blocking raises         |
| paragraph 2, sentence 3, p.    | biomass, and the result follows     |
| 12**                           | largely from the model structure:   |
|                                | (i) dispersal rates depend only on  |
| *"Blocking therefore raised    | distance, elevation and             |
| basin-mean native biomass at   | passability, not on density or on   |
| every upstream cost and both   | conditions at the destination (Eq.  |
| endpoints, with no change of   | 4; src/ode_model.jl,                |
| sign."*                        | precompute_dispersal_matrix), so    |
|                                | individuals keep leaving suitable,  |
|                                | well-stocked sites; (ii) the 56     |
|                                | DIFICIL sites have near-zero        |
|                                | capacity (§3.2, paragraph following |
|                                | Eq. 1), so immigrants there decline |
|                                | at the truncated rate and cannot    |
|                                | establish; (iii) the model is       |
|                                | deterministic (limitation 5), so    |
|                                | local populations never go extinct  |
|                                | and connectivity cannot act through |
|                                | rescue from local extinction.       |
|                                | Without this explanation the result |
|                                | may be read as evidence that        |
|                                | barriers benefit native fish.       |
|                                |                                     |
|                                | **Correction.** Add two or three    |
|                                | sentences in §4.4 stating these     |
|                                | mechanisms and that the result must |
|                                | not be interpreted as a benefit of  |
|                                | barriers.                           |
+--------------------------------+-------------------------------------+
| **§5.3 "Connectivity is        | **Issue.** An appropriate caveat,   |
| context dependent", paragraph  | but the structural features listed  |
| 1, sentence 3, p. 15**         | are not those that most directly    |
|                                | produce the result.                 |
| *"The absence of a reversal is |                                     |
| a finding about this model's   | **Correction.** Add                 |
| corrected structure, not a     | destination-independent dispersal,  |
| general ecological law: the    | the near-zero-capacity sites and    |
| fixed on-path tree, the        | the absence of                      |
| multiplicative obstacle rule   | extinction--recolonisation dynamics |
| and the uniform passability    | to the list.                        |
| assumption all shape it."*     |                                     |
+--------------------------------+-------------------------------------+
| **§5.6, paragraph 1, sentence  | **Issue.** In the model,            |
| 3, p. 16**                     | restricting movement through        |
|                                | passability raises trout biomass    |
| *"Brown trout remains a useful | ("Blocked passability left more     |
| sentinel because warming and   | trout biomass in the basin at every |
| restricted movement act in the | upstream cost and both endpoints",  |
| same direction for a sparse,   | §4.3, paragraph 3, sentence 3, p.   |
| headwater-associated           | 12; +71.9 units at cost 0.05,       |
| population."*                  | Figure 10), i.e. in the opposite    |
|                                | direction to warming. Only upstream |
|                                | movement cost acts weakly in the    |
|                                | same direction (355.4 → 353.3 units |
|                                | from cost 0.05 to 0.50). The Figure |
|                                | 10 caption itself says that warming |
|                                | and passability "act on trout       |
|                                | through different magnitudes".      |
|                                |                                     |
|                                | **Correction.** Specify which       |
|                                | restriction is meant (higher        |
|                                | upstream movement cost) or          |
|                                | rephrase, e.g. "...because it       |
|                                | responds to both warming and        |
|                                | hydraulic structure".               |
+--------------------------------+-------------------------------------+

**R6. Role of climate in the variance decomposition**

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **§4.4, paragraph 1, sentence 2, p. 12**          | **Issue.** Table 4, cited in the    |
|                                                   | same sentence, reports a climate    |
| *"Upstream cost dominated total biomass (SS =     | fraction of 0.142 for native        |
| 0.666), native biomass (SS = 0.729), native       | biomass in the interaction          |
| richness (SS = 0.935) and realised richness loss  | experiment (next row), so "minor    |
| (SS = 0.970); passability dominated invasive      | for every metric (SS ≤ 0.015)"      |
| biomass (SS = 0.741) and trout biomass (SS =      | holds only for the obstacle sweep.  |
| 0.951); and the climate endpoint was minor for    | The Abstract (paragraph 3, sentence |
| every metric (SS ≤ 0.015) (Table 4)."*            | 3) and the Conclusions (paragraph   |
|                                                   | 3, sentence 1) restrict the         |
|                                                   | statement correctly.                |
|                                                   |                                     |
|                                                   | **Correction.** could be "...minor  |
|                                                   | for every metric in the obstacle    |
|                                                   | sweep (SS ≤ 0.015), whereas it      |
|                                                   | accounted for 0.142 of native       |
|                                                   | biomass in the interaction sweep    |
|                                                   | (Table 4)".                         |
+---------------------------------------------------+-------------------------------------+
| **Table 4, row "Native biomass (interaction)", p. | **Issue.** Circle 1: climate        |
| 38**                                              | fraction of 0.142, an order of      |
|                                                   | magnitude above the obstacle-sweep  |
| ![](media/image2.png){width="4.803149606299213in" | values.                             |
| height="2.7729527559055116in"}                    |                                     |
|                                                   | **Correction.** Cite this value     |
|                                                   | where the role of climate is        |
|                                                   | summarised (§4.4).                  |
+---------------------------------------------------+-------------------------------------+

**\**

**R7. Figure--caption consistency**

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **Figure 6, panels and caption, p. 28**           | **Issue.** Circle 1: the x-axis is  |
|                                                   | realised warming over 2036--2045    |
| ![](media/image3.png){width="4.803149606299213in" | (0.18--1.29 °C). Circles 2--4: the  |
| height="5.624413823272091in"}                     | panels report "GCM-fixed" slopes    |
|                                                   | with R² = 0.19, 0.07 and 0.00.      |
|                                                   | Circle 5: the caption states that   |
|                                                   | the regressor is 2026--2045 warming |
|                                                   | (0.21--1.09 °C), "the same          |
|                                                   | regressor as Table 3". Circle 6:    |
|                                                   | the caption reports R² = 0.14, 0.06 |
|                                                   | and 0.01, which match the text      |
|                                                   | (§4.1, paragraph 2) but not the     |
|                                                   | panels.                             |
|                                                   |                                     |
|                                                   | **Correction.** Regenerate the      |
|                                                   | figure with the 2026--2045          |
|                                                   | regressor and pooled OLS fits       |
|                                                   | stated in the caption.              |
+---------------------------------------------------+-------------------------------------+
| **Figure 13, panels (a) and (c) and caption, p.   | **Issue.** Circle 1: panel (a)      |
| 35**                                              | reports −1.41 % per °C and R² =     |
|                                                   | 0.56. Circle 2: the panels use      |
| ![](media/image4.png){width="4.803149606299213in" | 2036--2045 warming. Circle 3: the   |
| height="3.7880511811023623in"}                    | caption describes the x-axis as     |
|                                                   | warming at 2045. Circle 4: the      |
|                                                   | caption reports 1.76 percentage     |
|                                                   | points per °C and R² = 0.63, the    |
|                                                   | value used in the text (2026--2045  |
|                                                   | regressor, §3.6). Circle 5: panel   |
|                                                   | (c) shows maxima up to about 150    |
|                                                   | days. Circle 6: the caption gives   |
|                                                   | 80--135 days. In addition, the      |
|                                                   | in-panel annotations overlap the    |
|                                                   | data in (a) and are cut at the      |
|                                                   | panel edge in (c).                  |
|                                                   |                                     |
|                                                   | **Correction.** Regenerate the      |
|                                                   | figure with the 2026--2045          |
|                                                   | regressor used in the text, and     |
|                                                   | align caption and panels (x-axis    |
|                                                   | description, slope, R² and exposure |
|                                                   | range).                             |
|                                                   |                                     |
|                                                   | Modify the figure to have a clarer  |
|                                                   | view of the data                    |
+---------------------------------------------------+-------------------------------------+
| **§4.3, paragraph 1, sentence 4, p. 11**          | **Issue.** The dose--response       |
|                                                   | relationship is shown in Figure 13  |
| *"Across the 44 runs, trout relative biomass      | (p. 35); Figure 10 (p. 32) shows    |
| declined by 1.76 percentage points of baseline    | the trout response to hydraulic     |
| per °C of realised warming (R² = 0.63), and trout | structure.                          |
| final biomass declined by 10.08 model units per   |                                     |
| °C (R² = 0.63) (Figure 10)."*                     | **Correction.** Change "(Figure     |
|                                                   | 10)" to "(Figure 13)".              |
+---------------------------------------------------+-------------------------------------+

**\**

**R8. Median warming in Figures 1 and 2**

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **Figure 1, label and caption, p. 23**            | **Issue.** Circle 1: "+0.62 °C →    |
|                                                   | 16.4 °C". Circle 2: described as    |
| ![](media/image5.png){width="4.803149606299213in" | the "ensemble-median warming to     |
| height="4.15146653543307in"}                      | 2045", without stating the          |
|                                                   | statistic or period. Figure 2 (next |
|                                                   | row) gives +0.74 °C for the         |
|                                                   | ensemble-median end-of-2045 warming |
|                                                   | and +0.58 °C for the 2026--2045     |
|                                                   | mean over the same runs; the Figure |
|                                                   | 1 value matches neither. (Circles   |
|                                                   | 3--5 are addressed in S5.)          |
|                                                   |                                     |
|                                                   | **Correction.** State which         |
|                                                   | quantity Figure 1 shows, and use    |
|                                                   | the same definition and value in    |
|                                                   | Figures 1 and 2.                    |
+---------------------------------------------------+-------------------------------------+
| **Figure 2, label and caption, p. 24**            | **Issue.** Circle 1: "+0.74 °C →    |
|                                                   | 16.6 °C". Circles 2--3: the caption |
| ![](media/image6.png){width="4.803149606299213in" | distinguishes the end-of-2045       |
| height="5.303477690288714in"}                     | warming (+0.74 °C) from the         |
|                                                   | 2026--2045 mean (+0.58 °C).         |
|                                                   |                                     |
|                                                   | **Correction.** Align with Figure 1 |
|                                                   | as indicated above; Figure 2 could  |
|                                                   | be removed, as it largely repeats   |
|                                                   | Figure 1.                           |
+---------------------------------------------------+-------------------------------------+
| **§4.1, paragraph 4 ("E21 warming uncertainty     | **\[Suggested\]**                   |
| band"), p. 10**                                   |                                     |
|                                                   | **Issue.** A further set of         |
| *"On the GuadeX network for the 2026--2045        | "ensemble-median" values            |
| horizon the central/ low/ high ensemble-median    | (0.58--0.73 °C) that differs from   |
| water warming by pathway is 0.580/0.368/0.741 °C  | the pathway means in Table 3        |
| (SSP1-2.6), 0.677/0.429/0.864 °C (SSP2-4.5),      | (0.54--0.70 °C) and from Figures    |
| 0.647/0.411/0.826 °C (SSP3-7.0) and               | 1--2, without explanation.          |
| 0.730/0.463/0.931 °C (SSP5-8.5)."*                |                                     |
|                                                   | **Correction.** Add one sentence    |
|                                                   | explaining how these E21 values     |
|                                                   | relate to those in Table 3 and      |
|                                                   | Figures 1--2.                       |
+---------------------------------------------------+-------------------------------------+

**R9. Interpretation of the aggregate dose--response slopes**

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **Abstract, paragraph 2, sentence 3, p. 3**       | **Issue.** The 44 runs are not      |
|                                                   | independent (they share 11 GCMs),   |
| *"Over the 44 scenario runs, end-of-2045 total    | and no test of the slopes is        |
| biomass changed by +0.075 model units per °C of   | reported. Even under ordinary least |
| realised warming (R² = 0.14), native richness by  | squares, which ignores that         |
| −0.0053 species per site per °C (R² = 0.06) and   | clustering, the richness and loss   |
| the deterministic realised richness loss by       | slopes are not significant (n = 44: |
| −0.00067 per °C (R² = 0.01)."*                    | p ≈ 0.11 and p ≈ 0.52, derived from |
|                                                   | the reported R²). The biomass slope |
|                                                   | is nominally significant (p ≈       |
|                                                   | 0.01), but its between-GCM SD is    |
|                                                   | twice the slope (Figure 5, next     |
|                                                   | row). §4.1 (paragraph 2, sentence   |
|                                                   | 2, p. 9) nevertheless reads them as |
|                                                   | directional ("Total biomass         |
|                                                   | increased with warming ... native   |
|                                                   | richness declined").                |
|                                                   |                                     |
|                                                   | **Correction.** Test the slopes     |
|                                                   | accounting for GCM (mixed model     |
|                                                   | with GCM as a random effect, or     |
|                                                   | GCM-clustered standard errors),     |
|                                                   | report the result, and interpret    |
|                                                   | the sign only where the test        |
|                                                   | supports it.                        |
+---------------------------------------------------+-------------------------------------+
| **Figure 5, panels (d)--(f), p. 27**              | **Issue.** Circle 1: richness slope |
|                                                   | −0.005 ± 0.028 (GCM SD), i.e. an SD |
| ![](media/image7.png){width="4.803149606299213in" | more than five times the slope.     |
| height="4.319378827646545in"}                     | Circle 2: biomass slope 0.075 ±     |
|                                                   | 0.16, i.e. an SD about twice the    |
|                                                   | slope. These SDs are shown in the   |
|                                                   | figure but not discussed in the     |
|                                                   | text.                               |
|                                                   |                                     |
|                                                   | **Correction.** Mention them in     |
|                                                   | §4.1 when the slopes are            |
|                                                   | interpreted.                        |
+---------------------------------------------------+-------------------------------------+
| **§6, limitation 15, p. 19**                      | **Issue.** This limitation already  |
|                                                   | establishes that the sign of the    |
| *"The magnitude of the weak-warming result is     | richness response is not robust,    |
| robust (basin native richness changes by at most  | but the Abstract does not say so.   |
| 0.0129 species, ≤ 10 of 775 sites, and the        | Values such as −0.006452 imply a    |
| warm-endpoint delta is numerically identical at − | precision that the uncalibrated     |
| 0.006452), but the sign of the ensemble-mean      | model cannot support.               |
| richness response depends on the structure:       |                                     |
| shared-K gives a mixed, slightly positive mean    | **Correction.** Carry this caveat   |
| (+0.0020; 26 positive, 4 zero, 14 negative),      | into the Abstract, and round        |
| whereas the alternative is uniformly negative     | reported values to three            |
| (44/44; mean −0.0071) and 30 of 44 runs flip      | significant figures throughout.     |
| sign."*                                           |                                     |
+---------------------------------------------------+-------------------------------------+

**R10. Scope of the corrections stated in the Abstract**

+--------------------------------+-------------------------------------+
| **Technical report --- exact   | **Issue and how to correct it**     |
| text and position**            |                                     |
+================================+=====================================+
| **Abstract, paragraph 1,       | **Issue.** The parenthesis implies  |
| sentence 6, p. 3**             | that all corrections C1--C7 and     |
|                                | E1--E22 were applied, whereas       |
| *"All results in this release  | Section 6 states that several were  |
| are recomputed from the Phase  | not: physiological thermal niches   |
| 2 post-fix model outputs       | (C6, limitation 1), calibration and |
| (corrections C1--C7, E1--E22,  | per-GCM controls (C3 and E2,        |
| applied to the regenerated     | limitation 2) and the air--water    |
| interaction table); they       | refit ("E20 remains a stated        |
| supersede the September 2026   | limitation", limitation 3).         |
| values, and the parser-fix,    |                                     |
| competition-structure and      | **Correction.** Write, for example: |
| passability-floor              | "(corrections C1--C7 and E1--E22 as |
| sensitivities are reported in  | implemented or deferred; see        |
| Section 4.7."*                 | Section 6)".                        |
+--------------------------------+-------------------------------------+

**R11. Internal inconsistencies**

+--------------------------------+-------------------------------------+
| **Technical report --- exact   | **Issue and how to correct it**     |
| text and position**            |                                     |
+================================+=====================================+
| **§4.7, paragraph 4            | **Issue.** 33 obstacles with        |
| ("Passability floor and        | upstream passability 0.1 give       |
| on-path-only obstacles"), p.   | 10⁻³³, the value given in §2.2      |
| 14**                           | (paragraph 2, sentence 5, p. 4) and |
|                                | limitation 8 (p. 18).               |
| *"The extreme link is          |                                     |
| 1.33.2--1.33.58, a parent edge | **Correction.** Change "10−34" to   |
| with 33 obstacles (combined    | "10−33".                            |
| upstream dam factor 10−34);    |                                     |
| parent-only leaves it          |                                     |
| untouched, while the floor     |                                     |
| raises it to 10−3, and it      |                                     |
| contributes only 0.42--1.12%   |                                     |
| of the dispersal rate added by |                                     |
| the floor."*                   |                                     |
+--------------------------------+-------------------------------------+
| **Title page, p. 1**           | **Issue.** §3.6 (paragraph          |
|                                | introducing the list of estimators, |
| *"Corrected October 2026       | p. 8) states that "All headline     |
| release"*                      | numbers in this report were         |
|                                | recomputed from the Phase 2         |
|                                | post-fix output files on            |
|                                | 2026-11-01".                        |
|                                |                                     |
|                                | **Correction.** Update the release  |
|                                | date on the title page.             |
+--------------------------------+-------------------------------------+

## 2. Suggested additions

+---------------------------------------------------+-------------------------------------+
| **Technical report --- exact text and position**  | **Issue and how to correct it**     |
+===================================================+=====================================+
| **§4.2 heading and Table 2 caption, pp. 10--11**  | **\[Suggested\]**                   |
|                                                   |                                     |
| *"Post-spin-up validation against the 2006--2009  | **Issue.** S1. The model is         |
| survey"*                                          | initialised with the same survey,   |
|                                                   | so high recall (pooled sensitivity  |
|                                                   | 0.975) is expected by construction. |
|                                                   | The text itself says that the       |
|                                                   | result is "not treated as a         |
|                                                   | validation success" (§4.2,          |
|                                                   | paragraph 2, sentence 4, p. 10).    |
|                                                   |                                     |
|                                                   | **Correction.** Rename the section  |
|                                                   | and Table 2 "Post-spin-up           |
|                                                   | consistency check against the       |
|                                                   | 2006--2009 survey".                 |
+---------------------------------------------------+-------------------------------------+
| **§3.3 "Dendritic dispersal and obstacles",       | **\[Suggested\]**                   |
| paragraph 3, sentence 1, p. 7**                   |                                     |
|                                                   | **Issue.** S2. "Blocked" is a ×0.1  |
| *"The four passability states were baseline       | attenuation, not a closure, which   |
| (multiplier 1.0), improved (1.5), reduced (0.5)   | may be misread in the management    |
| and blocked (0.1), with passability capped at     | discussion.                         |
| 1.0."*                                            |                                     |
|                                                   | **Correction.** Add "(strong        |
|                                                   | attenuation, not complete closure)" |
|                                                   | after "blocked (0.1)".              |
+---------------------------------------------------+-------------------------------------+
| **§3.2, paragraph following Eq. 2, sentence 4, p. | **\[Suggested\]**                   |
| 6**                                               |                                     |
|                                                   | **Issue.** S3. Limitation 1 covers  |
| *"Thermal optima were set to the midpoint of each | the general point. For brown trout  |
| empirical temperature range and thermal breadth   | the practical consequence is that   |
| to one-sixth of that range (T\* = 12.0 °C and σ = | heat-stress mortality (Eq. 3)       |
| 2.67 °C for brown trout, whose empirical range is | starts at 20 °C, whereas the        |
| 4--20 °C)."*                                      | literature places the upper limit   |
|                                                   | for feeding and growth near 19--20  |
|                                                   | °C and the upper incipient lethal   |
|                                                   | temperature near 25 °C (Elliott &   |
|                                                   | Elliott, 2010).                     |
|                                                   |                                     |
|                                                   | **Correction.** Add a sentence      |
|                                                   | noting that the trout upper limit   |
|                                                   | used in Eq. 3 is a sublethal        |
|                                                   | threshold, so heat-stress mortality |
|                                                   | is probably overestimated for this  |
|                                                   | species.                            |
+---------------------------------------------------+-------------------------------------+
| **Table 1 "Corrected simulation design", row      | **\[Suggested\]**                   |
| "Upstream costs", p. 8**                          |                                     |
|                                                   | **Issue.** S4. Sum-of-squares       |
| *"0.01, 0.05, 0.10 and 0.50 in the obstacle       | fractions depend on the levels      |
| experiment; 0.01 and 0.50 in the interaction      | explored: upstream cost spans a     |
| experiment"*                                      | 50-fold range, against two climate  |
|                                                   | endpoints over 20 years.            |
|                                                   |                                     |
|                                                   | **Correction.** Add a sentence in   |
|                                                   | §4.4 noting that the factor ranking |
|                                                   | is conditional on the ranges        |
|                                                   | explored.                           |
+---------------------------------------------------+-------------------------------------+
| **Figure 1, legend and caption, p. 23**           | **\[Suggested\]**                   |
|                                                   |                                     |
| ![](media/image5.png){width="4.803149606299213in" | **Issue.** S5. The caption states   |
| height="4.15146653543307in"}                      | that "the dominant cyprinids peak   |
|                                                   | at 19--20 °C" (circle 5), but the   |
|                                                   | curves peaking in that range        |
|                                                   | include *Aphanius baeticus*         |
|                                                   | (Aphaniidae; 19.0 °C group, circle  |
|                                                   | 3) and *Cobitis paludica*           |
|                                                   | (Cobitidae; 20.0 °C, circle 4),     |
|                                                   | which are not cyprinids.            |
|                                                   |                                     |
|                                                   | **Correction.** Replace "dominant   |
|                                                   | cyprinids" with "dominant native    |
|                                                   | species" in the caption.            |
+---------------------------------------------------+-------------------------------------+
| **§3.1 "Study domain and assemblage", paragraph   | **\[Suggested\]**                   |
| 1, sentence 4, p. 5**                             |                                     |
|                                                   | **Issue.** S6. *Liza ramada* is no  |
| *"Migratory or diadromous taxa are Anguilla       | longer the accepted name.           |
| anguilla, Liza ramada and Mugil cephalus."*       |                                     |
|                                                   | **Correction.** Replace with        |
|                                                   | *Chelon ramada* throughout.         |
+---------------------------------------------------+-------------------------------------+
| **§2.2 "The Guadalquivir Basin as a model         | **\[Suggested\]**                   |
| system", paragraph 1, sentence 6, p. 4**          |                                     |
|                                                   | **Issue.** S7. With zero local      |
| *"Migratory taxa remain in the dynamical system   | growth (§3.2, paragraph following   |
| and contribute to total biomass, but are excluded | Eq. 1) their basin totals cannot    |
| from native-richness reporting because their      | change (ΔN = 0.0 for AA, LR and MC  |
| estuarine or marine life histories are not        | in Table 2, p. 11); at most they    |
| represented by the modelled freshwater            | are redistributed among sites.      |
| corridor."*                                       |                                     |
|                                                   | **Correction.** Add that their      |
|                                                   | basin totals are held constant by   |
|                                                   | construction.                       |
+---------------------------------------------------+-------------------------------------+
| **§3.1, paragraph 1, sentence 3, p. 5**           | **\[Suggested\]**                   |
|                                                   |                                     |
| *"Exotic/invasive species are Gambusia holbrooki, | **Issue.** S8. Limitation 13 refers |
| Micropterus salmoides, Lepomis gibbosus, Cyprinus | to the "current species pool", but  |
| carpio, Carassius gibelio, Ameiurus melas,        | the pool is that of the 2006--2009  |
| Oncorhynchus mykiss, Esox lucius, Gobio lozanoi,  | survey, so exotics established in   |
| Tinca tinca and Alburnus alburnus."*              | the basin but not recorded in the   |
|                                                   | survey are excluded (e.g. *Silurus  |
|                                                   | glanis* in Iznájar Reservoir).      |
|                                                   |                                     |
|                                                   | **Correction.** State in §3.1 or    |
|                                                   | limitation 13 that the pool         |
|                                                   | reflects the 2006--2009 survey.     |
+---------------------------------------------------+-------------------------------------+
| **§4.6 "Other species and aggregate indicators",  | **\[Suggested\]**                   |
| paragraph 1, sentence 1, p. 13**                  |                                     |
|                                                   | **Issue.** S9. Table 2 (p. 11)      |
| *"Warm-adapted native cyprinids and most invasive | shows large changes in the          |
| and migratory taxa changed comparatively little   | no-warming control over 2026--2045  |
| over the study horizon in aggregate."*            | for most taxa, e.g. L. sclateri     |
|                                                   | −45.9 %, P. willkommii −66.9 %, M.  |
|                                                   | salmoides −87.9 % and G. holbrooki  |
|                                                   | +241.0 %. The sentence presumably   |
|                                                   | refers to scenario-minus-control    |
|                                                   | differences, but does not say so.   |
|                                                   |                                     |
|                                                   | **Correction.** Specify that the    |
|                                                   | statement refers to                 |
|                                                   | scenario-minus-control differences, |
|                                                   | and add below Table 2 that the      |
|                                                   | large control changes reflect the   |
|                                                   | uncalibrated transient.             |
+---------------------------------------------------+-------------------------------------+
| **§7 "Conclusions", paragraph 1, sentence 2, p.   | **\[Suggested\]**                   |
| 20**                                              |                                     |
|                                                   | **Issue.** S10. No per-species      |
| *"This weak aggregate response is not evidence of | scenario-minus-control responses    |
| thermal safety: it reflects cancellation between  | are reported except for trout, so   |
| warm-adapted dominant natives and a cold-water    | the cancellation is not shown (see  |
| specialist moving further above its optimum."*    | also §4.1, paragraph 3, sentence 3: |
|                                                   | "Aggregate metrics consequently     |
|                                                   | conceal opposing species            |
|                                                   | responses."). Trout occupies about  |
|                                                   | 9 % of sites and loses 1.1--1.5     |
|                                                   | percentage points to warming, which |
|                                                   | is unlikely to offset gains         |
|                                                   | elsewhere at the basin scale.       |
|                                                   |                                     |
|                                                   | **Correction.** Report the          |
|                                                   | per-species scenario-minus-control  |
|                                                   | deltas, or present the cancellation |
|                                                   | as a hypothesis ("may reflect").    |
+---------------------------------------------------+-------------------------------------+
| --- Not reported in §3.5 or §4.4, which the       | **\[Suggested\]**                   |
| authors' response to E7 cites as the evidence for |                                     |
| a Latin-hypercube/Sobol pilot.                    | **Issue.** S11. The response to E7  |
|                                                   | states that a global sensitivity    |
|                                                   | pilot was run                       |
|                                                   | (results/global_sensitivity_pilot), |
|                                                   | but the report does not mention it. |
|                                                   |                                     |
|                                                   | **Correction.** Summarise the pilot |
|                                                   | briefly, or state in Section 6 that |
|                                                   | a global sensitivity analysis       |
|                                                   | remains future work.                |
+---------------------------------------------------+-------------------------------------+
| **§3.4 "Temperature forcing and initial           | **\[Suggested\]**                   |
| conditions", paragraph 1, sentence 3, p. 7**      |                                     |
|                                                   | **Issue.** S12. The authors'        |
| *"The corrected basin-mean baseline site          | response to E23 states that the     |
| temperature across 775 sites is 15.82 °C (2026    | extrapolation of the                |
| annual means; range 9.79--18.06 °C)."*            | water-temperature model above the   |
|                                                   | highest calibration site (1,127 m;  |
|                                                   | 17 model sites) "remains a stated   |
|                                                   | limitation", but Section 6 does not |
|                                                   | mention it.                         |
|                                                   |                                     |
|                                                   | **Correction.** Add this limitation |
|                                                   | to Section 6.                       |
+---------------------------------------------------+-------------------------------------+

## References

Altermatt, F. (2013). Diversity in riverine metacommunities: a network
perspective. *Aquatic Ecology, 47*(3), 365--377.

Elliott, J. M., & Elliott, J. A. (2010). Temperature requirements of
Atlantic salmon *Salmo salar*, brown trout *Salmo trutta* and Arctic
charr *Salvelinus alpinus*: predicting the effects of climate change.
*Journal of Fish Biology, 77*, 1793--1817.
https://doi.org/10.1111/j.1095-8649.2010.02762.x

Perkin, J. S., & Gido, K. B. (2012). Fragmentation alters stream fish
community structure in dendritic ecological networks. *Ecological
Applications, 22*, 2176--2187. https://doi.org/10.1890/12-0318.1

Ramos-Merchante, A., & Prenda, J. (2018). The ecological and
conservation status of the Guadalquivir River Basin (S Spain) through
the application of a fish-based multimetric index. *Ecological
Indicators, 84*, 45--59. https://doi.org/10.1016/j.ecolind.2017.08.034

Tonkin, J. D., Poff, N. L., Bond, N. R., Horne, A., Merritt, D. M.,
Reynolds, L. V., Olden, J. D., Ruhi, A., & Lytle, D. A. (2019). Prepare
river ecosystems for an uncertain future. *Nature, 570*(7761), 301--303.
https://doi.org/10.1038/d41586-019-01877-1
