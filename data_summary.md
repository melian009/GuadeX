# Summary of the data

## File headers

### data/FishSizeMatrix.csv

```csv
CODIGO,ESPECIE,ESPECIE_1,LONGITUD_mm
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
* See also: `FishSizeMatrix_README.csv` for detailed variable definitions, species codes, and sampling notes.

### data/ConnectivityUTM.csv

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
* See also: `ConnectivityUTM_README.csv` for variable definitions and notes on network distances.

### data/BIOTIC/FishDensity_and_Juveniles_Matrix.csv

```csv
CODIGO,LS_DEN,SA_DEN,SP_DEN,PW_DEN,CP_DEN,IL_DEN,IO_DEN,AA_DEN,AH_DEN,GL_DEN,GH_DEN,LG_DEN,AAL_DEN,CG_DEN,CC_DEN,MS_DEN,ST_DEN,OM_DEN,EL_DEN,AM_DEN,TT_DEN,LR_DEN,MC_DEN,AB_DEN,Ls Juveniles,Sa Juveniles,Sp Juveniles,Pw Juveniles,Cp Juveniles,Il Juveniles,Aa Juveniles,Ah Juveniles,Gl Juveniles,Gh Juveniles,Lg Juveniles,Aal Juveniles,Cg Juveniles,Cc Juveniles,Ms Juveniles,St Juveniles,Om Juveniles,El Juveniles,Am Juveniles,Tt Juveniles,Lr Juveniles,Mc Juveniles,Ab Juveniles,
1.1.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,18.22414,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
1.1.3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72.631764,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
1.1.4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.43396,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
```

* Description: Standardized density data (abundance/surface area sampled) for each species at each sampling site, plus juvenile presence data.
  - `CODIGO`: Site identifier (matches other datasets)
  - `[Species]_DEN`: Standardized density for each species (abundance per unit area). Sites with only juveniles show 0 because juveniles weren't counted for density calculations
  - `[Species] Juveniles`: Juvenile presence codes (0 = no juveniles, 1 = only juveniles present, 2 = both juveniles and adults present)
  - **Species abbreviations include:**
    - **Native species:** Ls (*Luciobarbus sclateri*), Sa (*Squalius alburnoides*), Sp (*Squalius pyrenaicus*), Pw (*Pseudochondrostoma willkommii*), Cp (*Cobitis paludica*), Il (*Iberochondrostoma lemmingii*), Io (*Iberochondrostoma oretanum*), Ah (*Anaecypris hispanica*), St (*Salmo trutta*), Lr (*Liza ramada*), Mc (*Mugil cephalus*), Ab (*Aphanius baeticus*)
    - **Exotic species:** Gl (*Gobio lozanoi*), Gh (*Gambusia holbrooki*), Lg (*Lepomis gibbosus*), Aal (*Alburnus alburnus*), Cg (*Carassius gibelio*), Cc (*Cyprinus carpio*), Ms (*Micropterus salmoides*), Om (*Oncorhynchus mykiss*), El (*Esox lucius*), Am (*Ameiurus melas*), Tt (*Tinca tinca*)
    - **Uncertain status:** Aa (*Anguilla anguilla*)
* See also: `FishDensity_and_Juveniles_Matrix_README.csv` for species codes, juvenile coding, and sampling notes.

### data/BIOTIC/FishDensity_and_Juveniles_Matrix_README.csv

* Description: Metadata file explaining the density matrix structure, species codes, and data collection protocols.
  - Contains detailed species information with scientific names, common names (Spanish), and native/exotic status
  - Documents that 11 sites had only juveniles for *Luciobarbus sclateri*
  - Explains juvenile coding system and density calculation methods
  - **Note:** Eels from Guadiato basin (4 points in basin 14) are from fish farms and excluded from native/exotic calculations

### data/ABIOTIC/caracteristicas_peces_Guadalquivir_03-04-2018.csv

head of file:

```csv
SP.;MAX. SIZE (mm);FEED;SPAWNING ENVIRONMENT;ELEVATION (m);TEMPERATURE (�C); SALINITY;MOVEMENT RANGE (km);REPRODU. WITHIN THE BASIN;HOME RANGE
Ah;100;invertebrates, occasionally vegetation;shallow stony waters, occasionally plants;0 to 900;8 to 30;freshwater;short spring upstream journey to small streams;yes;Only Bembezar basin
Sa;140;invertebrates, occasionally vegetation;shallow stony waters, occasionally plants;0 to 900;8 to 30;freshwater;short spring upstream journey to small streams;yes;Most Guadalquivir
Sp;160;invertebrates, occasionally vegetation;shallow very oxigenated stony waters;0 to 1300;8 to 25;freshwater;short spring upstream journey to small streams;yes;Most Guadalquivir
Il;150;zooplancton, algae and detritus;low flow, gravel-sand and vegetation;0 to 900;8 to 30;freshwater;very short upstream journey;yes;Most Guadalquivir
Io;150;zooplancton, algae and detritus;low flow gravel bed;0 to 900;9 to 30;freshwater;very short upstream journey;yes;Only Jandula basin
Cp;140;benthic invertebrates, algae and detritus;low flow sandy-gravel bed;0 to 900;10 to 30;freshwater;stric sedentay;yes;Most Guadalquivir
St;400;invertebrates, occasionally small fish;shallow very oxigenated gravel bed;500 to 2500;4 to 20;freshwater;short winter upstream journey to small streams;yes;Only high mountain reaches
Ab;40;zooplancton, algae and detritus;quiet shallow loamy bed with vegetation;0 to 100;8 to 30;brackish;stric sedentay;yes;only high salinity waters close to the mouth
Aa;1000;omnivore, ocasionally small fish;ocean;0 to 900;8 to 30;brackish/freshwater;ocean;no;most Guadalquivir (originally), mouth to first dam (current)
```

* Description: Species ecological characteristics and life history traits.
  - Species codes match the density matrix
  - Includes maximum size, feeding habits, spawning requirements, elevation/temperature tolerances
  - Documents movement patterns and reproductive behavior within the basin
  - Shows current vs. historical distributions (e.g., eels now limited by dams)

### data/ABIOTIC/Matrix_Ambiental_Data.csv and data/ABIOTIC/Matrix_Ambiental_README.csv

head of the data file:

```csv
CODIGO;ANCHURA;FECHA;ESTACION;CODIGO_S;UTM X;UTM Y;EN_POZAS;OLOR_RESID;CANALIZADO;SIN_PECES;PRES/ABS_PECES;LONG.TR (m);SUP.TR(m�); % MAs;% MAe;% MAf;CONDUCTIV.;CVR (%);CLARIDAD;NIVEL;ACEITES;ESCORR.;ACC.GAN.;ACC.HUM.;TUB. Y P.;DEF.RIP.;EXTRACC.;GRAV.;COL.;BN;ZNPA;CSP;CR;AF;ZRR;ZUIAI;IET;PRES_RAPIDOS;N� R�PIDOS;SUP. R�PIDOS (m�);% R�PIDOS;R_INCR. (%);%R_BN;R_FONDO R.;R_ME�OS;R_CANTOS;R_GRAVAS;R_ARENAS;R_LIMOS;R_ARCILLAS;R_DETRITUS;R_CIENO;PRES_TABLAS;N� TABLAS;SUP. TABLAS (m�);% TABLAS;T_PROF.(m);T_FONDO R.;T_ME�OS;T_CANTOS;T_GRAVAS;T_ARENAS;T_LIMOS;T_ARCILLAS;T_DETRITUS;T_CIENO;PRES_POZAS;N� POZAS;SUP. POZAS (m�);% POZAS;P_PROF.(m);P_FONDO R.;P_ME�OS;P_CANTOS;P_GRAVAS;P_ARENAS;P_LIMOS;P_ARCILLAS;P_DETRITUS;P_CIENO;% estandarizado refugios en r�pidos;% estandarizado refugios en tablas;% estandarizado refugios en pozas;*PRES_R;D_PRES_R;E_PRES_R;G_PRES_R;H_PRES_R;*PRES_T;D_PRES_T;E_PRES_T;G_PRES_T;H_PRES_T;*PRES_P;D_PRES_P;E_PRES_P;G_PRES_P;H_PRES_P;Am valle_T(m);Am VR_T(m);% L VR_T;Dist.nac._S(m);Dist.Guadalq.(m);SHREVE arr.;N� emb. arr.;Mc arr.;Demb arr.(m);N� emb. Ab.;Mc ab.;Demb ab.(m);ALTITUD;OBSTRUCCIONES TRANSVERSALES;Pm;ESC;ETR;�REA_DP(Km�);LONGITUD_DP(Km);Pad;Pmad;% ZU;% ZICT;% ZMVC;% ZV;% ZAS;% ZAR;% B;% EVAH;% EASV;% ZH;% SA;�REA_SC (km�);LONGITUD_SC (Km);PENDIENTE_SC;PM_SC;Lap_SC (Km);DENSIDAD_SC (Km-1);Lr�o_SC(Km);EMBALSES_SC;ML_SC ;AT_SC;% ZU;% ZICT;% ZMVC;% ZV;% ZAS;% ZAR;% B;% EVAH;% EASV;% ZH;% SA;N� capt. Superf.;Volumen de agua (m3);N� pozos;N� pert. transversales;T�m_SC;Pm_SC;ESC_SC;ETR_SC
1.1.2;7,8;10.10.2006;OTO06;1,1;515058,43;4205080,41;NO;NO;NO;;1;120,00;1207,19;5;0;0;300;6,6;1;5;1;0;0;1;0;0;0;0;0;1;0;0;0;0;0;0;7,17;1;3;310,10;25,69;21,67;0,00;26,67;28,33;13,33;10,00;1,67;10,00;10,00;8,33;0,00;1;2;798,29;66,13;0,88;15,00;25,00;22,50;7,50;2,50;12,50;15,00;7,50;0,00;0;1;98,80;8,18;1,10;30,00;30,00;5,00;5,00;0,00;15,00;15,00;5,00;0,00;11,83;15,35;28,32;1;0;1;1;0;1;1;1;1;1;1;0;1;1;1;27,40;11,00;95;13860,00;5122,00;1;2;1414,00;4448;0;0,00;No existe;797;0;1081,00;570,00;518,00;58,00;15,78;7,22;6,66;0,00;0,00;0,00;0,00;0,00;0,00;44,70;55,21;0,00;0,00;0,09;95,86;21,12;6,03;5,96;0,00;0,00;9,57;2;1414,00;0,00;0,00;0,00;0,00;0,00;0,12;0,13;65,53;34,17;0,00;0,00;0,05;1;0;0;0;11,18;1092,46;379,84;554,57
1.1.3;5,9;12.10.2006;OTO06;1,2;516227,38;4211770,44;NO;NO;NO;;1;120,00;681,52;0;0;0;260;46,6;1;5;1;0;0;1;0;0;0;0;0;1;0;0;0;0;0;0;10,33;0;1;362,52;53,19;80,00;0,00;0,00;20,00;15,00;30,00;25,00;0,00;10,00;15,00;0,00;0;1;253,20;37,15;0,40;0,00;10,00;10,00;0,00;0,00;60,00;20,00;40,00;0,00;0;1;65,80;9,65;1,30;10,00;10,00;0,00;0,00;0,00;40,00;40,00;40,00;0,00;24,29;16,27;0,00;1;1;1;1;1;1;1;0;1;0;1;1;1;1;0;50,80;15,00;90;7335,00;2040,00;1;0;0,00;No existe;0;0,00;No existe;689;0;961,00;435,00;583,00;21,54;10,27;11,20;11,72;0,00;0,00;0,00;0,00;0,00;0,00;74,50;25,35;0,00;0,00;0,00;33,20;12,42;9,68;10,27;2,00;0,06;9,36;0;0,00;0,00;0,00;0,00;0,00;0,00;0,00;0,00;76,56;23,35;0,00;0,00;0,09;1;1000000;0;0;12,44;954,03;329,67;546,94
```

* Description: Comprehensive environmental data matrix covering physical, chemical, biological, and landscape variables for each sampling site. The README file provides detailed variable descriptions including:

**Key Variable Categories:**
  - **Site characteristics:** River width, sampling date, coordinates, fish presence/absence
  - **Water quality:** Conductivity, water clarity, oil presence, residential odors
  - **Physical habitat:** Mesohabitat types (rapids, riffles, pools), substrate composition, depth measurements
  - **Riparian conditions:** Vegetation cover, channel modifications, human disturbances
  - **Connectivity:** Distance to dams upstream/downstream, barriers, diversions
  - **Landscape context:** Land use percentages (urban, agricultural, forest), drainage area characteristics
  - **Climate variables:** Temperature, precipitation, evapotranspiration (1940-2005 averages)
  - **Refugia availability:** Fish shelter types across different mesohabitats
* See also: `Matrix_Ambiental_README.csv` for detailed variable explanations and metadata.

### data/ABIOTIC/Matriz_Ambiental_README.csv

* Description: Comprehensive metadata file explaining all 140+ environmental variables in the main environmental matrix.
  - **Mesohabitat classification:** Rapids (rápidos), riffles (tablas), and pools (pozas) with detailed substrate and flow characteristics
  - **Connectivity metrics:** Detailed dam impact measures including upstream/downstream distances and barriers
  - **Landscape variables:** Complete land use classification system for both drainage areas and subcatchments
  - **Refugia quantification:** Standardized measures of fish shelter availability by habitat type
  - **Temporal data:** Historical climate series and seasonal sampling information
  - Contains Spanish field descriptions with detailed ecological relevance explanations

### data/BIOTIC/Interacciones_peces_Guadalquivir_03-04-2018_ENG.csv

```csv
;Ah;Sa;Sp;Il;Io;Cp;St;Ab;Aa;Lr;Mc;Ls;Pw;Gl;Ms;Gh;El;Aal;Cc;Cg;Tt;Am;Om
Ah;;;;;;;;;;;;;;;;;;;;;;;
Sa;coexist, neutral.;;;;;;;;;;;;;;;;;;;;;;
Sp;coexist, neutral.;Coexist, affects  Sp through competition and sperm parasitism;;;;;;;;;;;;;;;;;;;;;
Il;coexist, neutral.;coexist, neutral.;coexist, neutral.;;;;;;;;;;;;;;;;;;;;
Io;No coexist.;coexist, neutral.;coexist, neutral.;Coexist, interfere through competition.;;;;;;;;;;;;;;;;;;;
Cp;coexist, neutral.;coexist, neutral.;coexist, neutral.;coexist, neutral.;coexist, neutral.;;;;;;;;;;;;;;;;;
St;No coexist.;Coexist, affects  Sa through predation.;Coexist, affects  Sp through predation.;Coexist, affects  Il through predation.;No coexist.;coexist, neutral.;;;;;;;;;;;;;;;;;
Ab;No coexist.;No coexist.;No coexist.;No coexist.;No coexist.;coexist, neutral.;No coexist.;;;;;;;;;;;;;;;;
Aa;Coexist, affects  Ah through predation.;Coexist, affects  Sa through predation.;Coexist, affects  Sp through predation.;Coexist, affects  Il through predation.;Coexist, affects  Io through predation.;Coexist, affects  Cp through predation.;No coexist.;Coexist, affects  Ab through predation.;;;;;;;;;;;;;;;
```

* Description: Fish species interaction matrix for the Guadalquivir basin. Each cell describes the ecological relationship between two species (e.g., coexistence, competition, predation, or no coexistence). The matrix is asymmetric and includes both native and exotic species. Useful for analyzing community structure, trophic interactions, and potential impacts of species introductions.
* See also: `Interacciones_peces_Guadalquivir_03-04-2018_ENG_README.csv` for species code mapping (abbreviations to scientific names).

### Shape files

* `muestreados_FINAL.shp`: Shapefile containing the sampling points with spatial attributes and connectivity information for GIS analysis.

## Data Content Summary

### Species Composition

Based on the density matrix and species characteristics data, the dataset includes:
- **Native species (Autóctono):** 12 species including *Squalius alburnoides*, *Luciobarbus sclateri*, *Salmo trutta*, *Pseudochondrostoma willkommii*, *Cobitis paludica*, *Iberochondrostoma lemmingii*, *I. oretanum*, *Anaecypris hispanica*, *Liza ramada*, *Mugil cephalus*, *Aphanius baeticus*

- **Exotic species (Exótico):** 11 species including *Oncorhynchus mykiss*, *Gobio lozanoi*, *Gambusia holbrooki*, *Lepomis gibbosus*, *Alburnus alburnus*, *Carassius gibelio*, *Cyprinus carpio*, *Micropterus salmoides*, *Esox lucius*, *Ameiurus melas*, *Tinca tinca*

- **Uncertain status:** 1 species (*Anguilla anguilla*)

### Spatial Coverage
The analysis in `dendritic.jl` shows:
- **Longitude range:** 180,000 - 580,000 UTM coordinates
- **Latitude range:** Variable (MINY to MAXY)
- **Spatial divisions:** Data analyzed in 10 longitudinal sectors for regional comparison

### What the Data Shows

#### 1. Multi-Scale Fish Community Data
The dataset provides three complementary perspectives on fish communities:
- **Individual-level data:** Body size measurements for demographic analysis
- **Population-level data:** Standardized density estimates for abundance comparisons
- **Community-level data:** Species presence/absence and juvenile recruitment patterns

#### 2. Comprehensive Environmental Characterization
The environmental matrix captures multiple scales of habitat variation:
- **Local habitat:** Mesohabitat structure, substrate composition, water quality
- **Reach-scale factors:** Riparian vegetation, channel modifications, flow characteristics
- **Landscape context:** Land use patterns, climate variables, watershed characteristics

#### 3. Dendritic Network Connectivity
Multiple connectivity measures quantify network fragmentation:
- **Dam impacts:** Distance to nearest upstream/downstream barriers
- **Network position:** Distance to main Guadalquivir river and stream order
- **Barrier density:** Number of obstructions per subcatchment and drainage area
- **Historical vs. current distributions:** Species range contractions due to fragmentation

#### 4. Species Ecological Requirements
Life history data enables mechanistic understanding of community patterns:
- **Size constraints:** Maximum body size ranges from 40mm (*Aphanius*) to 1000mm (*Anguilla*)
- **Habitat specialization:** From high mountain streams (*Salmo trutta*) to brackish estuaries (*Aphanius baeticus*)
- **Movement patterns:** From sedentary species to those requiring seasonal migrations
- **Thermal tolerances:** Temperature ranges reflecting altitudinal and seasonal preferences

#### 5. Anthropogenic Impact Assessment
Multiple disturbance gradients captured across sites:
- **Flow regulation:** Dam presence, water extractions, channel modifications
- **Land use intensification:** Agricultural conversion, urbanization, industrial development
- **Water quality degradation:** Pollution indicators, habitat simplification
- **Invasive species:** Exotic species establishment patterns relative to disturbance

This comprehensive dataset enables analysis of how multiple anthropogenic stressors interact across dendritic networks to influence native fish metacommunity dynamics in Mediterranean river systems, with particular focus on the role of dams in fragmenting populations and facilitating exotic species invasions.

## Note on Code Files

The codebase contains files from multiple projects:

### Current GuadeX Fish Project
- **Primary data files:** `FishSizeMatrix.csv`, `ConnectivityUTM.csv`, `FishDensity_and_Juveniles_Matrix.csv`, plus comprehensive environmental and species trait data
- **Analysis code:** `dendritic.jl` - analyzes native fish species patterns across longitudinal gradients

### Other Project (Plant Communities)
- **Referenced in:** `data_analysis.jl`
- **Missing data files:** `coccurrence.csv`, `multitrait.csv` - these files are from a different study on plant species co-occurrence and traits
- **Function library:** `functions.jl` - contains general analysis functions for Shannon diversity calculations

The `coccurrence.csv` and `multitrait.csv` files are not present because they belong to a separate plant ecology project that was analyzed using similar metacommunity frameworks.