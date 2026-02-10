# Assessing extinction dynamics in changing metacommunities: a case study of the Guadalquivir river basin

## Introduction
M
editerranean river basins, home to a significant number of endemic fish species facing conservation challenges, are among the most vulnerable types of dendritic metacommunities (Smith & Darwall, 2006; Reyjol et al., 2007). The structure and colonization-extinction dynamics of fish metacommunities are naturally constrained by the connectivity within their habitats, shaped by both local and regional factors (Brown & Swan 2010). These indigenous species are increasingly threatened by factors such as reservoir-induced connectivity reduction, increased pollution, habitat loss, and a surge in invasive species, which pose significant threats to freshwater ecosystems globally (Valerio et al. 2022).
Mediterranean rivers' vulnerability to climate change (Jarić et al. 2018) is exacerbated by their species-poor communities and periods of natural drought-flood (Gasith & Resh, 1999). Increasingly intense summer heat is compounding the fragmentation of fish populations, which naturally isolate into pools during the drought (Bernardo et al., 2003, Maghallaes et al., 2007). Although native species have adapted to these fluctuations, escalating climate change poses an unprecedented threat. The region has already experienced a 1.5ᵒC increase in summer temperatures, outpacing global rates by 20% (Ali et al., 2022). Predictions suggest that droughts will worsen, potentially leading to complete dry-outs of headwater streams. In a scenario where temperatures rise by 3ᵒC, there is projected to be a 12% reduction in expected mean rainfall (XIth IPCC Report, Ali et al., 2022). These conditions increase the risk of local extinction, with uncertain impacts on the metacommunity.
Despite these mounting threats, the primary drivers of extinction risks in these communities remain uncertain, particularly regarding the combined effects of multiple problems. To prevent catastrophic biodiversity loss in Mediterranean fish communities, it is essential to understand how extinction risk will vary across the basin under these compounded threats. It is essential not only to recognize the issues that have the most significant impact but also to consider how the interaction of multiple threats could magnify extinction risks in certain areas.
Thus, it is critical to develop robust, data-driven simulation tools that can identify and monitor at-risk regions. Such tools will inform targeted conservation efforts, including exotic species management, creation of refuge zones for native and endemic species' range expansion, and designation of priority areas for heightened protection and monitoring. It will be an invaluable resource to address some of the pressing objectives of the recent Spanish National Strategy for River Restoration 2022-2030, including climate change effects and selection of hydrological reserves (MITECO 2022).

## Aims

Our research aims to address these gaps through the following objectives:
1) What are the native fish extinction risk levels in the Guadalquivir River basin under the current threat scenario?
Using a large empirical metacommunity dataset (Fig. 1), we will model the probability of species extinction in the Guadalquivir basin. Our approach will incorporate empirical data, taking into account the unique disturbances within each sub-basin, such as the prevalence of exotic species, dam presence, and land use patterns. Findings will be compared with a 'null model' scenario, which instead uses 'slices' that follow longitude bands across the whole basin.
2) How will extinction risk change with increased threat levels?
We intend to study the influence of exotic species and climate change-driven desiccation on the range dynamics and extinction risk of native and endemic species. Species-specific biological parameters, i.e., migration, life history traits, will be used to parametrize the dynamics of the simulations.
3) What tools can we provide to predict extinction risks and assess resilience under different scenarios in Mediterranean basins?
We aim to develop a spatially-explicit modeling tool to enable practitioners and stakeholders to predict extinction risks and assess the resilience of various freshwater fish populations under varying drought and anthropogenic impact scenarios in the Guadalquivir. This tool will have a user-friendly interface, with the possibility of being adapted for use in other Mediterranean basins.

## Methods

To achieve these aims, we will employ the following methodologies:
Extinction probability curves
Using an existing database of fish communities, habitat characteristics, and river connectivity in the Guadalquivir river basin, we will model extinction probability curves influenced by different disturbance regimes.
This will be a two-step approach: 1. Random spatial longitude-based slices along the basin will be used to examine
the extinction patterns (Figure 2). 2. Each sub-basin will be used to set the limits for the curves' extent, to reflect the connected
reality of each watershed.

By comparing these two patterns, we aim to understand if there are significant differences in extinction patterns and relate them to current disturbances in each sub-basin. When comparing the spatial slices vs. sub-basin approach, patterns of greater intra-basin variability vs. general variability might emerge.We will calculate the extinction probability curves by:
1. Developing a map detailing the ratio of native to exotic species in each sampling point shown in Figure 1.
2. Generating profiles of native species proportions across the Guadalquivir basin using longitudinal spatial slices. A preliminary approach is shown in Figures 1 and 2.
3. Calculating the extinction probability curves using both the "spatial slices" and sub-basin methods to define curve limits.
4. Comparing the extinction probabilities and identify “cold spots” – regions with a notably high extinction probability. Preliminary analyses suggest local native extinction peaks at around 10% for all random longitudinal slices (Figure 2).
Figure 1: Sampled sites in the Guadalquivir basin (red). Coloured stripes in the background represent slices created along longitude-based bands used to select random sampling points for the null model of extinction probability. Selected sampled sites in each band are used for the density plot in Figure 2.

Figure 2: Density plot (density of sampling points) as a function of the fraction of native species (x-axis) in each longitudinal slice used for the null model, shown in Figure 1. In this hypothetical scenario, native populations collapse at around the 10% (0.1) mark in all cases and there is a high proportion of sites potentially at risk, with values between 30-10% (0.3 and 0.1).
Global change scenarios
Using existing data, we will conduct dynamic simulations to model the extinction probability of multiple species assemblages under different scenarios. These scenarios will encompass a range of factors, including:
• Dendritic fragmentation: We will assess the impact of fragmentation due to connectivity changes and the creation of reservoirs.
• Exotic species dynamics: Our model will simulate responses of native and endemic species to exotic species invasions, considering factors such as competition and predation.
• Climate change effects: We will account for increased drought and desiccation, reduction in channel width, and habitat loss, all of which are projected effects of climate change.
• Anthropogenic land-use changes: Our simulations will consider the gradient of anthropic uses as derived from previous results from our research group and the specific sensitivity of each species to these disturbances.

## Platform for knowledge transfer

The progress of our research will be tracked via a git repository with a final interface that visualizes all steps of the results. This software will form the core product for reporting and will facilitate efficient knowledge transfer during the project's timeline.