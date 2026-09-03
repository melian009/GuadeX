**the archive is predominantly for the Guadalquivir River Basin District**, but it is **not exclusively limited to it**.

## Geographic scope

The key identifier is the five-character `EUDHCod` or `EUMASCod` prefix:

- `ES050` = **Guadalquivir River Basin District**
- `ES150` = **Ceuta**
- `ES160` = **Melilla**

The archive contains:

| Layer | Guadalquivir (`ES050`) | Other areas |
|---|---:|---:|
| River/basin polygons | 467 of 468 | 1 record with a blank code |
| Groundwater bodies | 90 of 90 | None |
| Lakes/reservoirs | 97 of 97 | None |
| Transitional waters | 13 of 13 | None |
| Coastal waters | 3 of 9 | 3 Ceuta, 3 Melilla |
| Rivers/streams | 359 of 360 | 1 Melilla feature |

The three `ES050` coastal features are relevant to the Guadalquivir:

- **Pluma del Guadalquivir**
- **Parque Nacional de Doñana**
- **Doñana-Matalascañas**

The unrelated records are:

- Ceuta: **Bahía Norte**, **Puerto de Ceuta**, **Bahía Sur**
- Melilla: **Aguadú–Horcas Coloradas**, **Puerto de Melilla**, **Horcas Coloradas–Cabo Trapana**
- One river feature in Melilla: **Río de Oro**

So, for a Guadalquivir project, you would normally filter the data to `ES050`, while retaining the three `ES050` coastal records because they belong to the Guadalquivir district’s coastal area.

## What the data show

This is primarily a **GIS inventory and classification of water bodies**, not a water-quality-monitoring dataset. It shows the location, geometry, identity, type, and administrative/planning classification of water bodies.

### 1. Rivers and streams

The `SW_Line_4C_` layer contains **360 linear features**, classified as rivers. Examples include:

- Guadalquivir reaches
- Genil
- Guadiamar
- Guadajoz
- Viar
- Jandulilla
- Numerous arroyos and tributaries

It includes information such as:

- Water-body code and name
- River category and typology
- Whether the reach is natural, modified, or canalised
- Surface-water type codes such as `R-T08`, `R-T09`, etc.
- Permanence:
  - `PER` — permanent
  - `TMP` — temporary
  - `INT` — intermittent
- Latitude and longitude
- Length and area fields
- Related water bodies or zones
- Registration and deactivation dates
- Association with other hydrological features

### 2. Main river and basin water bodies

`Cuencas_masas_agua_4c` contains **468 polygon features**. These appear to represent the main delineated surface-water bodies, including rivers and streams represented as polygonal water-body areas.

Important fields include:

- `EUMASCod` — European water-body identifier
- `MAS_Nombre` — Spanish name
- `COD_SISEXP` — planning or management-system code
- `COD_SZONA` — sub-zone code
- `COD_ZONA` — zone code
- `PERI_CUENC` — basin perimeter
- `AREA_CUENC` — basin area
- `MASA_CICLO` — planning-cycle-related field
- Geometry-derived length and area

This layer is useful for mapping and grouping the basin’s surface-water bodies by hydrological zone.

### 3. Groundwater bodies

`GWB_4C` contains **90 groundwater bodies**, all coded `ES050`.

Examples include:

- Sierra de Cazorla
- Sierra de Cañete–Corbones
- Marismas de Doñana
- Manto Eólico Litoral de Doñana
- Aluvial del Guadalquivir – Curso Medio
- Cabra–Gaena
- Lebrija

The data show:

- Groundwater-body boundaries
- Groundwater codes
- Names in Spanish and English
- Associated surface-water bodies
- Geological and hydrogeological relationships
- Whether the body is associated with another water body
- Latitude, longitude, and area
- Date of inclusion in the dataset

For example, the `Associated` field indicates that **75 of the 90 groundwater bodies are associated with another feature**, while 15 are not.

### 4. Lakes and reservoirs

`SWB_Lago` contains **97 lake/reservoir features**, all within `ES050`.

Examples include:

- Embalse de Guadalén
- Embalse de Giribaile
- Embalse de Quiebrajano
- Embalse de Guadalmena
- Embalse de Canales
- Embalse del Negratín
- Embalse de La Bolera
- Embalses Doña Aldonza y Pedro Marín

The layer distinguishes between:

- Natural lakes
- Artificial water bodies
- Heavily modified water bodies
- Reservoirs formed from originally riverine water bodies
- Artificial reservoirs in areas that were not originally water bodies

In this layer:

- 61 features are classified as **very modified**
- 32 as **natural**
- 4 as **artificial**
- About 60 are described as reservoirs formed from water bodies that were originally rivers

It also includes:

- Lake/reservoir type
- Permanence
- Mean depth
- Surface area
- Associated groundwater or surface-water bodies
- Coordinates
- Relevant dates

### 5. Transitional waters

`SWB_Transicion` contains **13 transitional-water features**, all in `ES050`.

These cover the lower Guadalquivir and estuarine/marsh areas, including:

- Marismas de Bonanza
- Desembocadura Guadalquivir–Bonanza
- Guadiamar and Brazo del Oeste
- Brazo del Este
- Encauzamiento del Guadaira
- Dársena Alfonso XII
- Corta San Jerónimo–Presa de Alcalá del Río

All 13 are classified as:

- **Surface water**
- **Transitional water**
- **Very modified**
- **Permanent**
- Associated with another water body

This layer is particularly relevant for studying the Guadalquivir estuary, marshes, artificial channels, and tidal or transitional zones.

### 6. Coastal waters

`SWB_Costera` contains nine coastal-water features, but only three belong to `ES050`:

- Pluma del Guadalquivir
- Parque Nacional de Doñana
- Doñana-Matalascañas

These indicate how the dataset includes the coastal and marine-influenced area connected with the Guadalquivir district. It also records whether each coastal body is natural or heavily modified, its area, coordinates, water-body type, and related features.

## What it does **not** appear to show

Based on the available fields, the archive does **not** appear to contain:

- Measured chemical water-quality values
- Nutrient concentrations
- Pollutant concentrations
- Flow time series
- Rainfall or discharge measurements
- Groundwater levels
- Biological monitoring results
- Ecological-status scores
- Historical daily or monthly observations

Those may exist in separate monitoring or hydrological datasets. This archive mainly provides the **spatial framework**: the official water-body boundaries and their classifications, codes, types, relationships, and basic geometric attributes.