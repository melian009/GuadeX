# GuadalX: Assessing Extinction Dynamics in Changing Metacommunities of the Guadalquivir River Basin

## Project Summary

This project addresses the pressing challenge of **catastrophic biodiversity loss** in Mediterranean freshwater systems by studying the fish metacommunity dynamics of the **Guadalquivir River basin**. Mediterranean river basins are classified as among the **most vulnerable types of dendritic metacommunities**, home to a significant number of endemic fish species. We aim to develop **robust, data-driven simulation tools** to understand how the **combined effects of multiple threats** --including climate change, exotic species, and fragmentation-- drive extinction risk in these fragile ecosystems.

The final product will be a spatially-explicit modeling tool with a user-friendly interface, designed to inform targeted conservation efforts such as exotic species management and the selection of **hydrological reserves**.

## Research Gaps and Objectives (Aims)

Despite mounting threats (e.g., reservoir-induced connectivity reduction, increased pollution, and extreme climate change effects), the **primary drivers of extinction risks in these communities remain uncertain**. It is particularly critical to understand **how the interaction of multiple threats could magnify extinction risks** in certain areas.

Our research aims to address these critical gaps:

| Objective | Description |
| :--- | :--- |
| **Aim 1: Current Threat Scenario** | What are the native fish extinction risk levels in the Guadalquivir River basin under the **current threat scenario**? This models extinction probability incorporating unique disturbances like **dam presence, exotic species prevalence, and land use patterns** within each sub-basin. |
| **Aim 2: Increased Threat Levels** | How will extinction risk change with **increased threat levels**? This objective specifically studies the influence of **exotic species and climate change-driven desiccation** on the range dynamics and extinction risk of native and endemic species. |
| **Aim 3: Predictive Tool Development** | What tools can we provide to **predict extinction risks and assess resilience** under different scenarios in Mediterranean basins? This involves creating a spatially-explicit modeling tool for practitioners and stakeholders. |

- Questions:
  * What determines coexistence of native fish species within a "Pristine" Mediterranean (environmentally constrained) dendritic network?
  * How is the community dynamics and -ultimately - composition affected by the introduction of dams (large obstacles to connectivity) and exotic species?
  * Does the native community change when dams and/or exotics are introduced?
  * What is the probability of extinction  of native species after exotics have been introduced, and what (local) factors can increase that probability?
  * Asymmetric basin (current species assemblages are very different in right vs. left hand margin), does it play a role in the dynamics of coexistence?

## Methodology

The project uses a blend of comparative analysis and advanced dynamic simulations, leveraging a large existing dataset of the Guadalquivir river basin.

### 1. Extinction Probability Curve Analysis

We will model extinction probability curves influenced by different disturbance regimes using a two-step comparative approach:
1.  **Null Model:** Use **random spatial longitude-based slices** across the basin to examine general extinction patterns.
2.  **Connected Reality Model:** Use **each sub-basin** to set the limits for the curves' extent, reflecting the connected reality of each watershed.

By comparing these two patterns, we aim to understand if differences in extinction patterns relate to **current disturbances** (dam presence, exotic species prevalence, land use) in each sub-basin. This includes calculating the ratio of native to exotic species at each sampling point and identifying **"cold spots"**—regions with a notably high extinction probability.

### 2. Dynamic Simulations under Global Change Scenarios

We will conduct dynamic simulations using **species-specific biological parameters** (e.g., migration, life history traits) to model the extinction probability of multiple species assemblages. The model is based on niche modeling to predict species assemblages in **3D dendritic networks**, accounting for abiotic, biotic, migration, and fragmentation factors.

Scenarios will encompass:
*   **Dendritic fragmentation:** Assessing the impact of connectivity changes and the creation of **reservoirs**.
*   **Exotic species dynamics:** Simulating native and endemic species' responses to exotic species invasions, considering **competition and predation**.
*   **Climate change effects:** Accounting for increased **drought and desiccation**, reduction in channel width, and habitat loss.
*   **Anthropogenic land-use changes:** Considering the gradient of anthropic uses and the specific sensitivity of each species to these disturbances.

**Migration Modeling (3D Network):**
The model incorporates the physical reality of rivers in high-relief landscapes by extending the 2D network model to **3D**, where dispersal is a function of both hydrological distance and **elevation distance**. A critical parameter is the **upstream migration cost ($c$)**. Downstream migration is assumed to have no cost.

**Model-Data Comparison:**
We will use **Approximate Bayesian Computation (ABC)** or likelihood methods to infer the factors that best explain the empirical data patterns.

## Data Sources and Files

The research relies on a comprehensive dataset from 1,037 sampled sites in the Guadalquivir basin.

| File Name | Description and Key Variables | Data Context |
| :--- | :--- | :--- |
| `ConnectivityUTM.csv` | **Spatial and Network Structure**. Includes `ALTITUD` (Elevation), `CODIGO_S` (Subcatchment code), `EMBALSES_SC` (Number of dams per subcatchment), and distances to upstream/downstream dams. | Used to define the 3D dendritic network structure for fragmentation and migration modeling. |
| `Matrix_Ambiental_Data.csv` | **Environmental and Landscape Data**. Includes over 140 variables covering site characteristics, water quality (e.g., Conductivity), physical habitat (rapids/pools), riparian conditions, land use percentages, and historical climate averages. | Provides the abiotic and anthropogenic factors used in the niche models and simulation scenarios. |
| `FishDensity_and_Juveniles_Matrix.csv` | **Abundance and Demographics**. Contains **Standardized Density (`[Species]_DEN`)** and **Juvenile Presence Codes** for all 24 species (12 native, 11 exotic, 1 uncertain). | Essential for calculating current extinction risks (Aim 1) and assessing recruitment patterns and population resilience. |
| `caracteristicas_peces_Guadalquivir_03-04-2018.csv` | **Species Traits and Life History**. Includes max size, diet, spawning environment, elevation/temperature tolerance, and movement range. | Provides species-specific biological parameters required to parametrize dynamic simulations (Aim 2). |
| `Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv` | **Biotic Interactions**. Asymmetric matrix detailing relationships (coexistence, competition, predation) between native and exotic species. | Crucial for modeling the effects of **exotic species dynamics** (predation and competition) in simulations. |
| `FishSizeMatrix.csv` | **Individual-Level Data**. Fish body length in millimeters (`LONGITUD_mm`) for individual specimens. | Provides input for demographic analysis. |
| `muestreados_FINAL.shp` | **Spatial Points**. Shapefile of the 1,037 sampled sites for GIS analysis. | Used for visualizing the Guadalquivir basin and defining spatial slices. |

## Platform for Knowledge Transfer and Reproducibility

The final deliverable is a **spatially-explicit modeling tool** with a user-friendly interface to facilitate efficient knowledge transfer to practitioners and stakeholders.

*   **Software Language:** Codes will be developed using the **Julia computing language**.
*   **Version Control:** The progress of the research will be tracked via a **git repository**.
*   **Reproducible Research:** A **Jupyter notebook** will be created to make all research steps easily reproducible.

## Project Team and Timeline

The project is planned to run for **1 year**, from Fall 2023 until Fall 2024.

### Core Team:
*   Dr. Carlos Fernandez Delgado (Universidad de Córdoba, Spain)
*   Dr. Lucía Galvez Bravo (Liverpool John Moores University, UK)
*   Dr. Carlos J. Melián (EAWAG, ETH-Domain, Switzerland)
*   Dr. Ramón de Miguel Rubio (Tragsatec, Spain)
*   Dr. Ali R. Vahdati (University of Zurich, Switzerland)

### Key Activities (Timeline):
| Activity | Duration (Months) |
| :--- | :--- |
| Random and sub-basin patterns (Aim 1) | M1 – M4 |
| Global change scenarios (Aim 2) & Simulations | M5 – M9 |
| Numerical analysis & Robustness evaluation | M5 – M11 |
| Interface knowledge transfer (Aim 3) | M9 – M12 |
| Writing manuscript and reporting | M11 – M12 |
| Computer Scientist (Coding, maintenance, simulations) | Full Year |

## Original sources

* [Working draft](https://de.sharelatex.com/project/5ac1f0c58dd6a14ec01055e3)
* [data](https://drive.switch.ch/index.php/s/rNd3V73S2ca6MkT)
* [Other data](https://www.dropbox.com/sh/fndht5q3bxoyv05/AADjAs5uQrO5V4SnjvJC3pEZa?dl=0)
