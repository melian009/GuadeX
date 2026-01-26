# Modelling Approach for Fish Metacommunities in 3D Dendritic Networks

This document outlines a proposed modeling approach for simulating fish metacommunities in the Guadalquivir River basin. The goal is to create a site-specific Ordinary Differential Equation (ODE) model that can be "plugged into" a larger 3D dendritic network, capturing both local population dynamics and spatial dispersal processes.

By focusing on a site-specific ODE that can be networked, we are effectively building a **spatially-explicit modeling tool** designed to handle complex environmental and biological interactions.

## 1. Single Site Model: Parameters and Construction

To build this in Julia using `DifferentialEquations.jl`, the single-site model must be defined as a system of ODEs where the rate of change for a species population ($dN_s/dt$) is the sum of local dynamics and spatial flux.

### Structure of the model and its parameters

#### Local Population Dynamics

- Intrinsic Growth: Each species should have a growth rate (denoted as $\lambda$ or local birth rate).
- Biotic Interactions: incorporates an **asymmetric interaction matrix**. This matrix defines the relationships (competition, predation, or coexistence) between native and exotic species. For example, the presence of an exotic predator would negatively impact the $dN/dt$ of a native prey species. We have one matrix for the whole network, but each site will reference it for local interactions. Not all species may be present at each site.
- Environmental Filters: Population size is modified by abiotic factors. Key variables include:
  - Temperature: Species have specific "thermal windows" or tolerances.
  - Habitat Suitability: This includes river width, mesohabitat type (rapids, riffles, or pools), and substrate composition.
  - Anthropogenic Stressors: Factors like land use (urban/agricultural) and water quality (conductivity) act as "disturbance regimes" that influence extinction probability.

#### Spatial Flux (The "Plug-in" Mechanism)

To ensure sites can be "plugged into" a network, the model requires an **immigration** and **emigration** component based on the 3D network topology:

- Emigration: The rate at which individuals leave a site, often a function of current density and species-specific movement ranges.
- Immigration (The 3D Dispersal Rate): This is the sum of arrivals from neighboring nodes. The dispersal rate ($m_{ij}^k$) for species $k$ from site $j$ to site $i$ using a specific 3D extension formula:
  - Hydrological Distance: Dispersal is inversely related to the watercourse distance between nodes.
  - Upstream Migration Cost ($c$): If the immigration node is higher than the emigration node ($e_{im} > e_{em}$), a cost is applied: $x = (1 + c(e_{im} - e_{em}))^{-1}$.
  - Downstream Movement: Generally assumed to have **no cost** ($x = 1$).

#### Model Parameters for the AI Prompt

The following will be provided to the ODE model as parameters:

1. $N_{i,s}$: Population size of species $s$ at site $i$.
2. $T_i$: Temperature at site $i$.
3. $h_{i}$: Habitat suitability index at site $i$. This is called conductivity, a primary water quality parameter. High conductivity or chemical changes in the water are anthropogenic stressors that can influence the "specific sensitivity of each species" and their subsequent extinction probability. It should also include  mesohabitat types (rapids, riffles, or pools), substrate composition, and river width. In Mediterranean systems, the reduction in channel width and the "complete dry-outs of headwater streams" are critical factors in assessing if a site remains suitable for native fish.
4. $A_{s}$: Abiotic suitability can be calculated as a function of temperature, elevation, and habitat data (envi).
5. $B_{s,j}$: Interaction coefficient between species $s$ and species $j$ (from the interaction matrix).
6. $m$: Intensity of dispersal rate.
7. $c$: Upstream migration cost.
8. $d_{ij}$: Hydrological distance between site $i$ and $j$.
9. $e_i$: Elevation of site $i$.

#### Suggested Mathematical Form: The Spatially-Explicit Lotka-Volterra ODE

For modeling population sizes ($N$) of multiple species across a network of sites, the most robust form is a **Generalised Lotka-Volterra (GLV) model** extended with **environmental filtering** and **3D spatial flux**.

The rate of change for species $s$ at site $i$ should be structured as follows:

$$\frac{dN_{i,s}}{dt} = \text{Local Dynamics} + \text{Net Migration (Spatial Flux)}$$

##### A. Local Dynamics Component

The local component determines how the population grows and interacts within its specific site.

* **Mathematical Form:** $N_{i,s} \cdot [r_s(E_i) - \sum_{j} \alpha_{sj} N_{i,j}]$
* Parameters:
  * **$r_s(E_i)$ (Environmental Filtering):** The intrinsic growth rate must be a function of local environmental variables ($E_i$). Temperature, elevation, and conductivity are key variables. Use a **Gaussian-shaped "thermal window"** or tolerance function, as the species move when temperatures fall outside their "optimal aerobic window".
  * **$\alpha_{sj}$ (Biotic Interactions):** This uses the **asymmetric interaction matrix** provided in the sources. For example, if species $j$ is an exotic predator of native species $s$, $\alpha_{sj}$ will be a large positive value reducing the growth of $s$.

#### B. Spatial Flux Component (The "Plug-in" Mechanism)

This is the most critical part for your network model, as it connects site $i$ to its neighbors.
*   **Mathematical Form:** $\sum_{j \in \text{neighbors}} [m_{ij}^s N_{j,s} - m_{ji}^s N_{i,s}]$
*   **The 3D Dispersal Formula (Required for the AI):**
    A specific formula for the dispersal rate ($m_{ij}^s$) must be used to ensure the "3D" aspect is captured:
    $$m_{ij}^s = m \left( \frac{1}{d_{ij}} \right) \cdot x$$
    Where:
    *   **$d_{ij}$** is the hydrological distance.
    *   **$x = 1$** if moving downstream or to the same elevation ($e_{emigration} \ge e_{immigration}$).
    *   **$x = \frac{1}{1 + c(e_{im} - e_{em})}$** if moving upstream ($e_{im} > e_{em}$), where **$c$** is the **upstream migration cost**.

### Implementation details

1.  **State Vector ($u$):** Define the state as a matrix or flattened vector where $u[i, s]$ represents the population of species $s$ at site $i$.
2.  **Parameter Struct ($p$):** Create a structured object containing:
    *   The **Adjacency/Distance Matrix** (from `ConnectivityUTM.csv`).
    *   The **Elevation Vector** ($z_i$) for each site.
    *   The **Interaction Matrix** ($\alpha$).
    *   The **Upstream Cost Constant ($c$)**.
3.  **Boundary Conditions (Dams):** treat dams as **connectivity multipliers**. If a dam exists between site $i$ and $j$, the dispersal rate $m_{ij}$ should be multiplied by a "passability" factor (0.0 to 1.0) based on the "Demb" variables in the data.
4.  **Solver Recommendation:** the system will likely be **stiff** due to the combination of fast migration and slower growth rates, so it should use a solver like `Rodas5()` or `TRBDF2()` in `DifferentialEquations.jl`.
"