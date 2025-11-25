# Summary of the data

## File headers

### FishSizeMatrix.csv

```csv
CODIGO,ESPECIE,ESPECIE,LONGITUD_mm
1.1.10,Squalius alburnoides,Nativo,52
1.1.10,Squalius alburnoides,Nativo,54
1.1.10,Squalius alburnoides,Nativo,60
1.1.10,Squalius alburnoides,Nativo,61
1.1.10,Squalius alburnoides,Nativo,63
1.1.10,Squalius alburnoides,Nativo,63
1.1.10,Squalius alburnoides,Nativo,67
1.1.10,Squalius alburnoides,Nativo,67
```

* Description: Individual fish records with measurements. Each row represents a single fish specimen.
  - `CODIGO`: Site identifier
  - `ESPECIE`: Species name (scientific)
  - `ESPECIE_1`: Native/exotic status ("Nativo" or "Exótico")
  - `LONGITUD_mm`: Fish body length in millimeters

### ConnectivityUTM.csv

```csv
CODIGO,ALTITUD,CODIGO_S,EMBALSES_SC,Demb arr.(m),Demb ab.(m),Dist.Guadalq.(m),UTMX,UTMY
1.1.2,797,1.1,2,4448,No existe,5122,515058.43,4205080.41
1.1.3,689,1.2,0,No existe,No existe,2040,516227.38,4211770.44
1.1.4,716,1.1,2,8725,No existe,2100,512490.41,4206303.4
1.1.8,693,1.3,1,No existe,2199,2150,522862.3,4223459.53
1.1.9,673,1.4,1,No existe,1712,1859,524755.24,4231066.56
1.1.10,680,1.1,2,8497,No existe,1077,512435.4,4207324.4
1.1.11,671,1.4,1,No existe,1796,2057,523762.23,4231766.55
1.2.2,350,2.2,0,No existe,No existe,2821,468770.36,4198706.02
1.3.4,354,3,8,3725,6317,62023,475735.2,4220288.13
```

* Description: Site characteristics and spatial information for each sampling location.
  - `CODIGO`: Site identifier (matches FishSizeMatrix) - each site belongs to exactly one subcatchment
  - `ALTITUD`: Elevation in meters
  - `CODIGO_S`: Subcatchment code - identifies which tributary drainage area the site belongs to within the larger Guadalquivir basin. A subcatchment is a smaller watershed unit that drains into a larger river system, representing distinct hydrological units with their own connectivity patterns. **Note:** Multiple sites can belong to the same subcatchment (e.g., sites 1.1.2, 1.1.4, and 1.1.10 all belong to subcatchment 1.1), but each site belongs to only one subcatchment.
  - `EMBALSES_SC`: Number of dams/reservoirs within the subcatchment
  - `Demb arr.(m)`: Distance from the site to nearest upstream dam in meters ("No existe" = no upstream dam present)
  - `Demb ab.(m)`: Distance from the site to nearest downstream dam in meters ("No existe" = no downstream dam present)
  - `Dist.Guadalq.(m)`: Distance from the site to the main Guadalquivir river in meters
  - `UTMX`, `UTMY`: UTM coordinates (spatial location) of the site

### Shape file

* `muestreados_FINAL.shp`: Shapefile containing the sampling points with spatial attributes and connectivity information for GIS analysis.

## Data Content Summary

### Species Composition
Based on `ListNativeExotic.md`, the dataset includes:
- **Native species (Autóctono):** 16 species including *Squalius alburnoides*, *Luciobarbus sclateri*, *Salmo trutta*, etc.
- **Exotic species (Exótico):** 11 species including *Oncorhynchus mykiss*, *Gobio lozanoi*, *Gambusia holbrooki*, etc.
- **Uncertain status:** 1 species (*Anguilla anguilla*)

### Spatial Coverage
The analysis in `dendritic.jl` shows:
- **Longitude range:** 180,000 - 580,000 UTM coordinates
- **Latitude range:** Variable (MINY to MAXY)
- **Spatial divisions:** Data analyzed in 10 longitudinal sectors for regional comparison

### What the Data Shows

#### 1. Native Species Diversity Patterns
The analysis reveals longitudinal gradients in native species diversity:
- Native species richness varies across the basin
- Fraction of native species ranges from 0 to 0.65 across sites
- Different patterns emerge when examining density distributions by longitude

#### 2. Dendritic Network Structure
The connectivity data captures:
- **Dam impacts:** Distance to nearest upstream/downstream dams
- **Network position:** Distance to main Guadalquivir river
- **Elevation effects:** Altitude range from 350m to 797m
- **Subcatchment connectivity:** Number of barriers per subcatchment - each subcatchment represents a distinct tributary system with varying levels of fragmentation (0-8 dams per subcatchment)

#### 3. Fish Community Structure
Individual fish records provide:
- **Size distributions:** Body length measurements per species
- **Species occurrence:** Presence/absence at each site
- **Native-exotic composition:** Status classification for each species

#### 4. Disturbance Gradients
The dataset captures multiple anthropogenic impacts:
- **Fragmentation:** Dam presence and distances
- **Species introductions:** Native vs. exotic species occurrence
- **Spatial heterogeneity:** Varying conditions across the dendritic network

This comprehensive dataset enables analysis of how dendritic network structure, dam fragmentation, and exotic species introductions affect native fish metacommunity dynamics in Mediterranean river systems.

## Note on Code Files

The codebase contains files from multiple projects:

### Current GuadeX Fish Project
- **Primary data files:** `FishSizeMatrix.csv`, `ConnectivityUTM.csv`, `ListNativeExotic.md`, `muestreados_FINAL.shp`
- **Analysis code:** `dendritic.jl` - analyzes native fish species patterns across longitudinal gradients

### Other Project (Plant Communities)
- **Referenced in:** `data_analysis.jl` 
- **Missing data files:** `coccurrence.csv`, `multitrait.csv` - these files are from a different study on plant species co-occurrence and traits
- **Function library:** `functions.jl` - contains general analysis functions for Shannon diversity calculations

The `coccurrence.csv` and `multitrait.csv` files are not present because they belong to a separate plant ecology project that was analyzed using similar metacommunity frameworks.