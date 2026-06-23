# Mathematical Formulation of the Fish Metacommunity Model

This document provides a detailed description of the spatially-explicit modeling framework used to simulate fish metacommunity dynamics in the Guadalquivir River basin. The model is implemented as a system of Ordinary Differential Equations (ODEs) that integrate local population dynamics with logistic growth constraints, environmental filtering, species-specific dispersal, and governance interventions.

## 1. General Model Structure

The rate of change for the population size $N_{i,s}$ of species $s$ at site $i$ is defined by the sum of local dynamics and net spatial flux:

$$\frac{dN_{i,s}}{dt} = \underbrace{f_{local}(N_{i,s}, \mathbf{N}_i, E_i)}_{\text{Local Dynamics}} + \underbrace{\sum_{j \in \text{neighbors}} (m_{ji} N_{j,s} - m_{ij} N_{i,s})}_{\text{Spatial Flux}}$$

where:
- $N_{i,s}$ is the population density of species $s$ at site $i$.
- $\mathbf{N}_i$ is the vector of all species populations at site $i$.
- $E_i$ represents the local environmental conditions.
- $m_{ij}$ is the dispersal rate from site $i$ to site $j$.

## 2. Local Dynamics with Logistic Growth

The local component combines logistic growth (constrained by carrying capacity) with biotic interactions, modified by environmental suitability filters.

### 2.1. Population Growth with Carrying Capacity

The local dynamics are expressed as:

$$\frac{dN_{i,s}}{dt} = N_{i,s} \left( r_{i,s}^{eff} \cdot \underbrace{\left(1 - \frac{\sum_j N_{i,j}}{K_i}\right)}_{\text{Logistic Term}} + \underbrace{\frac{\sum_j \alpha_{sj} N_{i,j}}{K_i}}_{\text{Interaction Term}} \right)$$

where:
- $r_{i,s}^{eff}$ is the effective intrinsic growth rate.
- $K_i$ is the site-specific carrying capacity for total biomass.
- $\alpha_{sj}$ is the interaction coefficient representing the effect of species $j$ on species $s$ (asymmetric interaction matrix).

The logistic term $(1 - \sum_j N_{i,j}/K_i)$ ensures that population growth slows as total biomass approaches the carrying capacity, implementing density-dependent regulation at each site.

### 2.2. Environmental Filtering

The effective growth rate $r_{i,s}^{eff}$ is determined by the base growth rate $r_s$ and a site-specific suitability factor $W_{i,s}$:

$$r_{i,s}^{eff} = r_{s} \cdot W_{i,s}$$

The suitability factor $W_{i,s}$ incorporates thermal tolerance and habitat quality:

$$W_{i,s} = \text{ThermalFilter}(T_i, \text{opt}_s, \sigma_s) \cdot h_i$$

#### Thermal Filter
We use a Gaussian-shaped "thermal window" to model species-specific temperature tolerances:

$$\text{ThermalFilter}(T_i, \text{opt}_s, \sigma_s) = \exp\left( -\frac{(T_i - \text{opt}_s)^2}{2\sigma_s^2} \right)$$

where:
- $T_i$ is the temperature at site $i$.
- $\text{opt}_s$ is the optimal temperature for species $s$.
- $\sigma_s$ is the thermal tolerance (standard deviation) of species $s$, derived from the temperature range as $\sigma_s \approx \text{range}/6$.

The parameter $\sigma_s$ (thermal breadth) controls the shape of the thermal niche:
- At $T_i = \text{opt}_s$, the filter equals 1.0 (maximum growth), regardless of $\sigma_s$.
- At $T_i = \text{opt}_s \pm \sigma_s$, the filter drops to $e^{-1/2} \approx 0.607$.
- At $T_i = \text{opt}_s \pm 2\sigma_s$, the filter drops to $e^{-2} \approx 0.135$.
- A smaller $\sigma_s$ means a narrower thermal niche (specialist species: growth decays rapidly as temperature deviates from optimum).
- A larger $\sigma_s$ means a broader thermal niche (generalist species: tolerant of wider temperature ranges).
- Two species with the same optimum but different $\sigma_s$ at a given temperature will experience different growth reductions. For example, with $T_i - \text{opt}_s = +3^\circ\text{C}$: a species with $\sigma_s = 3^\circ\text{C}$ grows at 61% of maximum, while a species with $\sigma_s = 6^\circ\text{C}$ grows at 88%.

#### Habitat Suitability
- $h_i$ is the habitat suitability index at site $i$, derived from the Trophic State Index (IET). Higher IET values indicate lower suitability, and $h_i$ is normalized to range from 0.1 to 1.0.

## 3. Species-Specific Dispersal (3D Dendritic Dispersal)

The spatial component models the movement of individuals through the river network, accounting for hydrological distance, elevation changes, and physical barriers. The model implements species-specific dispersal using a two-stage approach.

### 3.1. Base Dispersal Matrix

The base dispersal rate $m_{ij}$ from site $i$ to site $j$ is defined as:

$$m_{ij} = D_{\text{median}} \cdot \frac{x_{ij} \cdot p_{ij}}{d_{ij}}$$

where:
- $D_{\text{median}} \approx 0.0274 \ \text{km/day}$ is the median daily dispersal speed (derived from literature annual dispersal rates: 10 km/year ÷ 365 days/year).
- $d_{ij}$ is the hydrological distance between site $i$ and $j$ (in km).
- $x_{ij}$ is the 3D connectivity factor (elevation cost).
- $p_{ij}$ is the dam passability factor (ranging from 0.0 to 1.0).

Species-specific dispersal rates are obtained by scaling the base matrix per species: $m_{ij}^s = \text{dispersal\_scaling}_s \cdot m_{ij}$, where $\text{dispersal\_scaling}_s = D_s / D_{\text{median}}$.

### 3.2. Elevation Cost (Upstream vs. Downstream)

The cost of movement depends on the relative elevation of the sites:

- **Downstream or Equal Elevation ($e_j \le e_i$):**
  $$x_{ij} = 1.0$$
- **Upstream ($e_j > e_i$):**
  $$x_{ij} = \frac{1}{1 + c(e_j - e_i)}$$

where:
- $e_i$ and $e_j$ are the elevations of the source and destination sites, respectively.
- $c$ is the upstream migration cost constant (default: 0.05).

### 3.3. Species-Specific Dispersal Scaling

The base dispersal matrix uses a median species dispersal rate. Species-specific scaling factors convert this to species-specific rates:

$$\text{dispersal\_rate}_{s,ij} = \text{dispersal\_scaling}_s \cdot m_{ij}$$

The `dispersal_scaling` vector contains relative dispersal rates normalized so that the median species has a scale factor of 1.0. Literature-derived annual dispersal rates (km/year) are converted to these relative scaling factors. For example:
- Migratory species (e.g., *Pseudochondrostoma willkommii*, *Luciobarbus sclateri*): scale > 5
- Highly sedentary species (e.g., *Anaecypris hispanica*, *Aphanius baeticus*): scale < 0.1
- Invasive generalists (e.g., *Gambusia holbrooki*): scale > 4

### 3.4. Dam Passability

Dams act as semi-permeable barriers that reduce effective dispersal rates. The passability factor $p_{ij}$ is set to 0.1 for sites with downstream dams and 1.0 otherwise. This can be modified through governance interventions (see Section 6).

## 4. Data Preparation

The model relies on a comprehensive data preparation pipeline, implemented in `src/data_preparation.jl`, to integrate heterogeneous datasets into the ODE framework.

### 4.1. Species Characteristics (`load_species_characteristics`)

Loads data from `data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv`:
- Thermal optima are calculated as the midpoint of temperature ranges (`TEMPERATURE_C`)
- Thermal breadth ($\sigma_s$) is approximated as $\sigma_s \approx \text{range}/6$
- Max size (`MAX_SIZE_mm`) is available for reference
- Species codes in the characteristics file use title case (e.g., `Ah`, `Sa`), while the density matrix uses uppercase (e.g., `AH`, `SA`). The lookup uses case-insensitive matching to correctly assign per-species thermal parameters.

### 4.2. Site Data (`load_site_data`)

Merges connectivity data (`data/ConnectivityUTM.csv`) and environmental data (`data/ABIOTIC/Matriz_Ambiental_Data.csv`) on site codes (`CODIGO`). Missing values in dam distance columns (`Demb arr.(m)`, `Demb ab.(m)`) are handled, and data is converted to numeric formats.

### 4.3. Interaction Matrix (`load_interaction_matrix`)

Parses qualitative interaction descriptions from `data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv` into a quantitative matrix using `parse_interaction_string`. Values are assigned based on interaction type:
- "No coexist": -1.0
- "Displaces": -0.8
- "Predation": -0.5
- "Competition" or "Interfere": -0.3
- "Affects" (general): -0.2
- "Coexist" or "Neutral": 0.0

### 4.4. Spatial Connectivity (`build_distance_matrix`)

Constructs a sparse distance matrix from `data/Matrix_distances_1037puntos_BRUTO_FINAL.csv` by:
1. Grouping sites by subcatchment (`CODIGO_S`)
2. Sorting sites within each subcatchment by distance to the main river (`Dist.Guadalq.(m)`)
3. Connecting only adjacent sites within the same subcatchment
4. Connecting outlet sites (closest to main river) across subcatchments based on elevation

### 4.5. Environmental Filtering

- **Site temperatures (`extract_site_temperatures`):** Uses `TEMP_MEDIA_SC` (mean temperature of subcatchment) if available, otherwise estimates from elevation using a lapse rate of 6.5°C per 1000m.
- **Habitat suitability (`extract_habitat_suitability`):** Derived from the Trophic State Index (IET), where higher IET values indicate lower suitability.

### 4.6. Intrinsic Growth Rates (`build_intrinsic_growth_rates`)

Base growth rates are derived from literature values (annual rates converted to daily: $r_{daily} = r_{annual}/365$). Species present at a site use full growth rates; absent species use 10% of the rate to allow potential colonization.

### 4.7. Carrying Capacity (`build_carrying_capacity`)

Site-specific carrying capacity $K_i$ is derived from observed total fish density, scaled by a factor of 10 to represent maximum sustainable biomass. A floor value (50.0 or 10th percentile of non-zero values) ensures minimum capacity for sites with low observed density.

### 4.8. Dispersal Matrix (`precompute_dispersal_matrix`)

A sparse dispersal matrix $M$ is pre-computed incorporating:
- Hydrological distance
- Elevation-dependent migration costs (upstream cost constant $c$)
- Dam passability factors (0.1 for dammed sites, 1.0 otherwise)

## 5. Implementation Details

- **State Vector:** The system state $u$ is a flattened vector representing a matrix of size $N_{sites} \times N_{species}$, where $u[i, s]$ is the population of species $s$ at site $i$.
- **Sparse Matrix Optimization:** Dispersal rates are pre-calculated into a sparse matrix $M$, where $M_{ji}$ represents the rate from $i$ to $j$.
- **Numerical Integration:** The system is solved using `Tsit5()` (explicit Runge-Kutta method) with adaptive time stepping and relative/absolute tolerances of $10^{-6}$. A positivity callback projects negative populations to zero after each step.
- **Time Units:** 1 time unit = 1 day. Default simulation is 1 year (365 days).

## 6. Governance: Modeling Management Interventions

The model includes a governance framework to simulate the effects of management interventions on fish populations. Two primary levers are implemented: **exploitation pressure** and **dam passability**.

### 6.1. Exploitation (Harvesting/Mortality Pressure)

Exploitation factors reduce intrinsic growth rates to simulate fishing pressure, bycatch mortality, or other anthropogenic sources of mortality.

**Formulation:**
$$r_{i,s}^{eff,managed} = r_{i,s}^{eff} \cdot e_i$$

where $e_i$ is the exploitation factor for site $i$ (ranging from 0.0 to 1.0):
- $e_i = 1.0$: No exploitation (baseline)
- $e_i = 0.5$: 50% exploitation (growth rate reduced by half)
- $e_i = 0.0$: Complete closure (no reproduction)

**Implementation:**
- Exploitation is applied at the subcatchment level via `build_exploitation_vector()`
- A dictionary maps subcatchment IDs to exploitation factors
- Factors are applied in-place to the intrinsic growth rate matrix via `apply_exploitation_to_growth_rates!()`

**Example scenario:**
```julia
exploitation_scenario = Dict{Float64, Float64}(
    3.1 => 0.7,   # 30% exploitation in subcatchment 3.1
    4.2 => 0.5,   # 50% exploitation in subcatchment 4.2
)
```

### 6.2. Dam Passability (Connectivity Restoration)

Dam passability multipliers modify the effective passability of dams to simulate fish passage improvements (e.g., fish ladders, fishways) or degradation (dam repairs, blockages).

**Formulation:**
$$p_{ij}^{managed} = \min(1.0, p_{ij} \cdot \pi_j)$$

where $\pi_j$ is the passability multiplier for origin site $j$:
- $\pi_j = 1.0$: Baseline passability (no management)
- $\pi_j > 1.0$: Improved passability (e.g., 1.5 = 50% improvement via fish ladder)
- $\pi_j < 1.0$: Reduced passability (e.g., 0.5 = 50% reduction, blocked dam)

**Implementation:**
- Passability is applied at the subcatchment level via `build_dam_passability_vector()`
- A dictionary maps subcatchment IDs to passability multipliers
- The modified dam matrix is used to recompute the dispersal matrix via `precompute_dispersal_matrix()`

**Example scenario:**
```julia
passability_scenario = Dict{Float64, Float64}(
    2.0 => 1.5,   # Improve passability by 50% in subcatchment 2.0
    3.1 => 2.0,   # Double passability (fish ladder) in subcatchment 3.1
)
```

### 6.3. Combined Management Scenarios

Governance scenarios can combine exploitation and passability interventions to test trade-offs between harvest pressure and habitat connectivity:

```julia
# Example: Moderate fishing in tributaries with improved connectivity
combined_scenario = (
    exploitation = Dict(1.1 => 0.8, 1.2 => 0.8, 2.2 => 0.7),
    passability = Dict(3.0 => 1.5, 6.0 => 1.5)
)
```

### 6.4. Expected Outcomes

- **Increased exploitation:** Generally reduces population sizes, with stronger effects on species with low intrinsic growth rates.
- **Improved dam passability:** Enhances connectivity, benefiting migratory species (high dispersal scaling) more than sedentary species. Can facilitate range expansion and rescue effects for fragmented populations.
- **Combined effects:** Management interventions may have non-additive effects due to density-dependent feedbacks and species interactions.

## 7. Model Assumptions

1. **Constant Interaction Matrix:** Biotic interaction coefficients ($\alpha_{sj}$) are assumed constant across the basin, though realized interactions depend on local densities.
2. **Passive/Active Dispersal Balance:** Dispersal is modeled as a density-dependent flux proportional to the population at the source site.
3. **Gaussian Thermal Niche:** Species performance decays symmetrically as temperatures deviate from the optimum.
4. **Dam Passability:** Dams act as semi-permeable filters that reduce effective dispersal rates between connected nodes.
5. **Logistic Growth:** Population growth is constrained by site-specific carrying capacities derived from observed fish densities.
6. **Species-Specific Dispersal:** Dispersal rates vary among species based on literature-derived relative scaling factors.
7. **Subcatchment-Level Governance:** Management interventions are applied uniformly within subcatchments, allowing for spatial targeting of conservation efforts.
