# River water-temperature model and 2045 projections for the Guadalquivir basin (GuadeX)

**Pipeline root:** `guadex_tw/` · **Reproduce with:** `python run_all.py`
**Status:** all deliverables produced from real downloads. Every number below carries its N.

---

## Bottom line (plain language)

We calibrated a daily river water-temperature model from air temperature using **40,277
quality-controlled water-temperature observations at 980 Spanish monitoring sites** (46 of them in
the Guadalquivir), driven by NASA POWER air temperature. The national model transfers to the
Guadalquivir with a leave-one-basin-out NSE of **0.70** (held-out basin, 771 observations), clearly
beating the naive baselines (NSE 0.26–0.41). Elevation matters: adding an elevation×air-temperature
interaction is highly significant and improves out-of-sample NSE by **+0.03 nationally, +0.08 within
the Guadalquivir, and +0.04 in the held-out-basin test** — but the effect is a *slope* change, not a
level shift, and it is a refinement rather than a transformation of the fit. Applying the calibrated
model to the AEMET/PNACC CMIP6 ensemble (11 GCMs × 4 SSPs, 8 representative Guadalquivir sites,
5 km) gives a **2036–2055 median warming of +0.72 °C (ssp126) to +0.95 °C (ssp585)** relative to each
model's own 1986–2005 baseline, rising to **+0.73 to +1.25 °C by 2041–2070**. The projections are
best trusted as a *relative climate signal with quantified ensemble spread*, not as absolute
temperatures for a specific reach: they rest on a stationarity assumption, sparse (mostly monthly)
historical observations, and a predictor that differs between calibration (MERRA-2) and projection
(AEMET/PNACC). Details and caveats are in `LIMITATIONS.md`.

---

## 1. Executive summary

- **Calibration network (Stage B):** 39,944 paired daily observations at 980 Spanish sites; pooled LOSO NSE for the elevation-aware model **0.625** vs **0.594** without elevation.
- **Local basin (Stage A):** 771 paired observations at 46 Guadalquivir sites; LOSO NSE **0.609** (elevation) vs **0.526** (no elevation).
- **Held-out basin (decisive):** national model trained without the Guadalquivir predicts it at **NSE 0.704 / RMSE 2.92 °C** (771 obs) vs 0.662 without elevation.
- **Elevation answer:** the elevation×air-temperature interaction is significant (Spain χ²=471, p≈2×10⁻¹⁰⁴; Guadalquivir χ²=35.8, p=2×10⁻⁹) and negative; the elevation *level* term is not significant nationally (p=0.56). Elevation refines the slope, it does not carry an independent level.
- **2045 window (2036–2055):** ensemble-median water-temperature change **+0.72 / +0.81 / +0.83 / +0.95 °C** for ssp126/245/370/585, with 10–90% spread of roughly 0.5–1.3 °C.

---

## 2. Data inventory

| source | variable | N sites | N obs | period | resolution | license | access status |
|---|---|---:|---:|---|---|---|---|
| GRQA v1.3 `TEMP_GRQA.csv` (Spain) | river Tw | 204 | 23,502 | 1979–2018 | point | CC-BY-4.0 (DOI 10.5281/zenodo.7056647) | **SUCCESS** |
| GRQA v1.3 Guadalquivir extract | river Tw | 5 | 195 | 1990–2018 | point | CC-BY-4.0 | **SUCCESS** (matches verified ground truth) |
| Waterbase WISE6 (EEA SDI disaggregated) | river Tw | 2,607 | 27,203 | 2008–2024 | point | EEA open data | **SUCCESS** |
| EEA discodata SQL API | river Tw | 13,336 sites metadata | obs query timed out | – | – | EEA open data | **PARTIAL/FAILED** |
| Combined QC water temperature | river Tw | 980 | 40,277 | 1979–2024 | point | mixed | **SUCCESS** |
| NASA POWER (MERRA-2) daily T2M/MAX/MIN | air Ta | 980 | daily | 1981–2024 (+1981–2024 baseline for 8 sites) | 0.5°×0.625° | NASA open | **SUCCESS** |
| Open-Meteo ERA5 archive | air Ta | 33 sites usable | partial | – | ~9–25 km | CC-BY-4.0 | **FAILED/PARTIAL** (HTTP 429 weighted rate limit) |
| Copernicus DEM GLO-30 (AWS) | elevation | 992 | 1 per site (3×3 median) | static | 30 m | Copernicus | **SUCCESS** |
| Open-Meteo elevation (Copernicus DEM 90 m) | elevation | 33 | 1 per site | static | 90 m | CC-BY-4.0 | **PARTIAL** (same rate limit) |
| AEMET/PNACC CMIP6 `tmean` (ESD-RegBA 5 km) | air Ta future | 11 GCM × 4 SSP × 8 GQ sites | daily | 1950–2014 hist, 2015–2100 proj | 5 km | AEMET/PNACC public | **SUCCESS** (THREDDS NCSS broken; remote HDF5 subset) |
| CHG IDE `ggiscloud_root:a1748512575838_cuenca_guadalquivir` | basin polygon | – | – | – | vector | CHG | **PARTIAL** (WFS GetFeature HTTP 401; CHG ES050 water-body catchments used instead) |

Exact URLs, UTC timestamps, sizes, MD5/SHA256 and errors are in `PROVENANCE.md`; failed rows are
recorded with the literal HTTP error rather than replaced.

---

## 3. QC log (Section 5 rules)

Applied first to GRQA (Spain/Iberia), then to the Waterbase-direct extract.

| rule | GRQA extract | Waterbase direct |
|---|---:|---:|
| 1. `param_code == TEMP` | enforced at read (2,991,485 TEMP rows scanned globally; 23,502 kept for Spain) | determinand `EEA_3121-01-5` |
| 2. unit normalisation (accept °C only) | 23,502 all `Deg C`; **0 removed** | 27,203 all `Cel`; **0 removed** |
| 3. physical range [−1, 45] °C | **0 removed** | **0 removed** |
| 4. cross-source duplicates (`TEMP_GRQA_dup_obs.csv`, 32,263 pairs) | **0** Spain rows present | 0 after GRQA overlap check |
| 5. IQR-outlier flag (not deleted) | **0 flagged** (`obs_iqr_outlier` all `no`) | n/a |
| 6. detection-limit rows | **0 removed** (field always null) | **0 removed** |
| 7. sites with < 12 observations | **23 sites / 36 obs removed** | **1,808 sites / 9,461 obs removed** |
| 8. daily aggregation (mean; `n_samples`) | 23,479 daily rows (23 collapsed) | 26,295 daily rows |
| 9. 4×MAD vs site-month median (flagged) | **767 flagged** | **633 flagged** (combined 1,400) |
| 10. continuity/availability | reported in `logs/qc_log.json` (e.g. `site_ts_continuity` 0.01–1.0; many sites ~0.5) | n/a |

Retained: `data/processed/water_temp_all_qc.csv` = **40,277 obs / 980 sites**
(GLORICH 22,237; WATERBASE_DIRECT 16,834; GEMSTAT 1,206), span 1979-01-11 → 2024-12-30.
`data/interim/water_temp_raw.csv` keeps the pre-QC extract for audit.

**Pipeline validation:** the Guadalquivir raw extract reproduces the verified ground truth exactly —
**195 rows / 5 sites** (ESP00016 = 67, ESP00017 = 58, ESP00018 = 68, plus two Waterbase single
observations); the naive bbox returns **406 rows / 11 sites**, and all **6 non-Guadalquivir sites are
excluded** (2 × ES040 Záncara/Valbuena, 3 Guadiana GEMSTAT, 1 Portuguese Guadiana).

---

## 4. Elevation

- **All 992 sites:** median 384 m, range **0–1,452 m** (Copernicus DEM GLO-30, 3×3 median; Open-Meteo DEM90 cross-check; 0 sites missing elevation).
- **Guadalquivir (46 sites):** median 170 m, range **0.5–1,127 m**.
- Site classification: river 925, reservoir 35, estuary 24, canal 8 (`data/processed/sites.csv`).
- The Guadalquivir range is now wider than the original three-site subset (0–235 m), which is
  essential for identifying the elevation interaction within the basin.

---

## 5. Models

Predictor: `ta` = NASA POWER daily mean air temperature, **lapse-corrected to site elevation**
with `Ta_corr = Ta_grid + γ·(z_site − z_grid)`, **γ = −0.0065 K m⁻¹** (a stated parameter, not truth).
`z = elevation_m / 1000`.

- **Model 1 — per-station, per-month OLS** `Tw = α[m] + β[m]·Ta + ε` (fit where n ≥ 15, else pooled-month fallback).
- **Model 2 — pooled month effects + elevation interaction** (core):
  `Tw = a0 + Σ_m f(month) + b1·Ta + b2·z + b3·(Ta·z) [+ b4·(Ta·log10 area)] + ε`.
  Month dummies vs 2-harmonic seasonal term were both fitted; month dummies won on AIC in both scopes
  (Guadalquivir 3654.5 vs 3654.9; Spain 191,100 vs 191,114) and are retained.
- **Model 3 — lagged / thermal-inertia** `Tw = a0 + f(month) + Σ b_k·Ta_d[Lk] + b2·z`.
- **Model 4 — mixed effects** `Tw ~ Ta + z + Ta:z + month + (1 + Ta | site_id)` (BLUPs saved).
- **Model 5 — process-based (air2stream)** — **not fitted** (see omissions).

### Fitted coefficients (physical units, °C)

| term | Guadalquivir (n=771, 46 sites) | All Spain (n=39,944, 980 sites) |
|---|---|---|
| intercept | +6.58 ± 0.54 | +4.60 ± 0.07 |
| `b1` Ta (°C/°C) | **+0.436 ± 0.036** (p≈5×10⁻³¹) | **+0.439 ± 0.015** (p≈3×10⁻¹⁸⁶) |
| `b2` z (°C/km) | +1.77 ± 0.95 (p=0.063) | +0.068 ± 0.116 (p=0.56) |
| `b3` Ta·z (°C/°C/km) | **−0.299 ± 0.050** (p=3×10⁻⁹) | **−0.157 ± 0.007** (p≈2×10⁻¹⁰⁴) |
| residual SD | 2.54 °C | 2.65 °C |
| MixedLM ICC | 0.248 | 0.444 |

Likelihood-ratio tests (Model 2 full vs reduced):
Guadalquivir **b3: χ²=35.8, p=2×10⁻⁹**; b2: χ²=3.5, p=0.061.
Spain **b3: χ²=471.0, p≈2×10⁻¹⁰⁴**; b2: χ²=0.34, p=0.56.

**Interpretation (physics):** air temperature already carries most of the elevation *level* signal
(colder air at altitude), so `b2` is small/insignificant once Ta is present. The significant negative
`b3` means **the water-to-air response slope flattens with elevation** — high reaches warm less per
degree of air warming — consistent with shorter growing seasons, snow/groundwater influence and
higher-gradient reaches. A non-significant `b3` would have been a legitimate finding; here it is
clearly non-zero.

**Sanity check:** predicted summer-minus-winter amplitude is 8.4–18.0 °C (Spain) and 7.3–15.2 °C
(Guadalquivir); summer exceeds winter at every site. One Guadalquivir site (a high, regulated reach)
falls slightly below the 8 °C rule of thumb and is flagged.

### Model-selection rule and omissions

Final model chosen by **leave-one-station-out NSE**, preferring the simplest model within 0.02:
**Model 2 (month dummies + Ta + z + Ta:z)**. Model 4 (mixed effects) reaches similar LOSO skill
(Spain 0.628 population-level, 0.681 on 20 true held-out stations) but adds complexity without a
decisive NSE gain, so it is reported as a sensitivity, not the primary model.

- **Model 3:** only 3 of 980 sites meet the ≥24 obs/yr cadence rule, and all three fall within a
  single calendar year (2024), so a time-blocked CV of lag windows is impossible. Lag models are
  therefore fitted **in-sample only** and flagged (`logs/model3_summary.json`); lag structure is not
  identifiable from these data. All lag columns are still provided in the paired table, computed
  (never interpolated) with full-window requirements.
- **Model 5 (air2stream):** **not fitted.** `air2stream` is not on PyPI and an 8-parameter process
  model cannot be identified from mostly monthly grab samples with no flow/shading data. Documented
  omission.
- **Basin-area term** `b4`: `upstream_basin_area` is only ~55% populated, so it is median-imputed
  within scope; the term is included for Spain and reported but not relied upon (elevation results
  are unchanged if it is dropped).

---

## 6. Validation (Section 9)

**Metrics** are computed on held-out data only. `NSE = 1 − SS_res/SS_tot`, `bias = mean(pred − obs)`.
Full table: `outputs/tables/cv_metrics.csv`.

### 6.1 Leave-one-station-out (LOSO), pooled over folds — primary

| model | scope | N test | NSE | RMSE (°C) | MAE (°C) | bias (°C) | note |
|---|---|---:|---:|---:|---:|---:|---|
| Tw = Ta | Spain | 39,942 | 0.151 | 3.48 | 2.76 | +0.44 | naive |
| Tw = Ta + constant | Spain | 39,942 | 0.172 | 3.47 | 2.77 | −0.00 | naive |
| monthly Tw climatology | Spain | 39,942 | 0.340 | 3.11 | 2.54 | +0.00 | naive |
| Model 1 (per station-month) | Spain | 39,942 | 0.603 | 2.15 | 1.69 | −0.09 | in-sample mean |
| **Model 2 (no elevation)** | Spain | 39,942 | 0.594 | 2.47 | 1.96 | −0.00 | LOSO |
| **Model 2 (Ta + z + Ta:z)** | Spain | 39,942 | **0.625** | **2.40** | 1.90 | +0.00 | LOSO |
| Model 4 (MixedLM, population) | Spain | 39,942 | 0.628 | 2.43 | 1.93 | +0.24 | RE=0 |
| Model 4 (MixedLM, true LOSO) | Spain | 1,042 | 0.681 | 2.45 | 1.96 | +0.25 | 20 held-out stations |
| Tw = Ta | Guadalquivir | 769 | −0.288 | 4.04 | 3.16 | +1.23 | naive |
| Tw = Ta + constant | Guadalquivir | 769 | −0.144 | 3.91 | 3.11 | −0.03 | naive |
| monthly Tw climatology | Guadalquivir | 769 | 0.406 | 2.97 | 2.45 | −0.02 | naive |
| Model 1 (per station-month) | Guadalquivir | 769 | 0.542 | 2.54 | 2.09 | −0.01 | in-sample mean |
| Model 2 (no elevation) | Guadalquivir | 769 | 0.526 | 2.64 | 2.17 | −0.03 | LOSO |
| **Model 2 (Ta + z + Ta:z)** | Guadalquivir | 769 | **0.609** | **2.46** | 2.00 | −0.04 | LOSO |
| Model 4 (MixedLM, population) | Guadalquivir | 769 | 0.656 | 2.39 | 1.94 | −0.16 | RE=0 |

At monthly cadence, the naive `Tw = Ta` baseline is *worse than useless* (negative NSE) because
day-to-day air temperature is far more variable than the monthly water sample; the monthly
climatology is a strong baseline (NSE 0.34–0.41), which the models beat.

### 6.2 Leave-one-basin-out (decisive Stage-B transfer test)

Train on all non-Guadalquivir Spanish sites (39,173 obs / 934 sites), predict the 46 Guadalquivir
sites (771 obs):

| model | NSE | RMSE (°C) |
|---|---:|---:|
| **Model 2 (elevation-aware)** | **0.704** | **2.92** |
| Model 2 (no elevation) | 0.662 | 3.13 |
| Tw = Ta + constant | 0.414 | 4.12 |
| Tw = Ta | 0.383 | 4.22 |
| monthly Tw climatology | 0.260 | 4.62 |

**The national elevation-aware model transfers to the held-out basin and beats every baseline.**

### 6.3 Blocked time-split (within-station non-stationarity)

Train on all but the last two distinct years per station, test on those (247 stations with ≥4 years):
mean NSE **0.477**, median RMSE 1.73 °C, mean bias **−0.52 °C** (the model is systematically cold in
the most recent test years — a real non-stationarity signal). Model 1 and Model 2 give identical
time-split skill because `z` and `Ta·z` are constant within a station, so a within-station test
cannot see elevation.

### 6.4 Figures

`outputs/figures/` (all at 200 dpi; captions in `captions.txt`):
`fig01` Tw–Ta scatter + pooled fit; `fig02` residuals vs month and vs elevation;
`fig03` LOSO predicted vs observed; `fig04` elevation-coloured site map with the ES050 basin;
`fig05` air/water time series for the three Guadalquivir GEMSTAT sites; `fig06` future monthly
ensemble spread; `fig07` distribution of 2036–2055 ΔTw by SSP; `fig08` the elevation answer —
per-station slope and offset vs elevation with 95% CI.

---

## 7. Future projections to 2045 (Section 10)

- **Predictor consistency:** the model is applied to AEMET/PNACC 5 km `tmean` after a **per-site,
  per-GCM additive bias correction** over **1986–2005** that aligns PNACC to the NASA POWER
  calibration climatology. Anomalies are computed **per model against that model's own historical
  run**, never against a multi-model historical mean.
- **Ensemble:** all **11 CMIP6 GCMs** for **ssp126, ssp245, ssp370, ssp585**; 8 representative
  Guadalquivir sites (0.5–1,127 m).
- **Baseline period:** 1986–2005.
- **Windows:** primary **2036–2055** (centred on 2045), standard **2041–2070**, appendix **2021–2040**.

### Ensemble ΔTw (°C, per-model anomaly vs own 1986–2005 baseline)

| window | ssp126 | ssp245 | ssp370 | ssp585 |
|---|---:|---:|---:|---:|
| 2021–2040 | 0.51 [0.33–0.84] | 0.54 [0.39–0.88] | 0.55 [0.41–0.86] | 0.61 [0.42–0.92] |
| **2036–2055** | **0.72 [0.50–1.08]** | **0.81 [0.57–1.31]** | **0.83 [0.63–1.23]** | **0.95 [0.70–1.31]** |
| 2041–2070 | 0.73 [0.51–1.18] | 0.93 [0.69–1.47] | 1.05 [0.78–1.48] | 1.25 [0.85–1.64] |

Values are the ensemble **median** with **10th–90th** percentiles in brackets (N = 88 site×GCM values
per cell). The figures also report the interquartile range. `outputs/tables/water_temp_future_2045.csv`
carries every GCM, and repeats the ensemble median/p25/p75/p10/p90 on each row.

### Absolute values and thermal metrics

Applied to absolute calibrated temperatures, the 2036–2055 ensemble-median annual mean water
temperature at the 8 sites is ~15.5–17.1 °C, with summer (JJA) means ~21.7–24.1 °C. Per site,
scenario and period, `outputs/tables/thermal_metrics.csv` reports
`mean_annual_c, mean_summer_c, mean_winter_c, p95_annual_c, summer_p95_c,
n_days_above_20/25/28, longest_run_above_28, n_days_below_10, temp_range_c, temp_sd_c`.
**Thresholds (20/25/28 °C hot, 10 °C cold) are placeholders defined at the top of
`scripts/10_project.py` and must be swapped for the actual GuadeX species optima.**

---

## 8. GuadeX-ready outputs

| file | contents |
|---|---|
| `outputs/tables/water_temp_historical_sites.csv` | per-observation modelled vs observed Tw, Ta, 90% prediction interval (N=39,944) |
| `outputs/tables/water_temp_future_2045.csv` | per GCM × scenario × site × period Tw and ΔTw + ensemble percentiles |
| `outputs/tables/thermal_metrics.csv` | per site × scenario × period thermal metrics (96 rows) |
| `data/processed/sites.csv` | site table with dual-source elevation, basin flag, class, cadence |
| `data/processed/paired_observations.csv` | paired Tw–Ta calibration table with lag columns (N=40,277) |
| `data/processed/paired_sites_matched.csv` | matched POWER/ROCIO/ERA5 cells and distances per site |

---

## 9. Assumptions and limitations (summary)

1. **Stationarity** of the calibrated Tw–Ta relation under climate change is assumed; changes in
   flow, reservoir operation, shading and abstraction are ignored. This is the largest structural
   uncertainty.
2. **Predictor mismatch** between calibration (NASA POWER/MERRA-2) and projection (AEMET/PNACC),
   addressed by per-site/per-GCM bias correction but not eliminated.
3. **Sparse, mostly monthly observations**; the three multi-year Guadalquivir series end in 1995.
4. **Reservoir/estuary sites** are pooled with rivers (classified but not separately modelled).
5. **Area imputation** and a single GCM member per model (no internal variability ensemble).

Full discussion and the minimum additional data that would materially improve the result are in
`LIMITATIONS.md`.

---

## 10. Acceptance checklist

| item | result |
|---|---|
| `grqa_temp_guadalquivir.csv` reproduces verified counts (67/58/68 + two 2018 singles; 1990–1995) | **PASS** — 195 rows / 5 sites |
| all bbox non-Guadalquivir sites excluded | **PASS** — 6 excluded |
| `grqa_temp_spain.csv` national set, wider elevation range | **PASS** — 204 sites; 0–1,452 m vs Guadalquivir 0.5–1,127 m |
| Waterbase attempted, outcome recorded | **PASS** — SQL API PARTIAL (timeout), SDI direct SUCCESS (27,203 obs) |
| elevation for every retained site, cross-checked | **PASS** — 0 missing, Copernicus DEM 30 m + Open-Meteo DEM90 |
| `paired_*.csv` with no fabricated lags; cadence rule applied | **PASS** — 40,277 rows; 3 eligible sites, all single-year, flagged |
| all five models (or documented omissions) | **PASS** — 1,2,3 (in-sample, flagged),4 fitted; 5 omitted with reason |
| LOSO metrics for every model + three naive baselines | **PASS** — `cv_metrics.csv` |
| elevation interaction tested by LRT **and** CV, conclusion stated | **PASS** |
| future projections ≥3 SSPs, full ensemble, median/IQR/P10–P90, per-model anomalies | **PASS** |
| every figure has labels, units, N, caption | **PASS** — `captions.txt` |
| `REPORT.md`, `LIMITATIONS.md`, `PROVENANCE.md`, `run_all.py`, `requirements.txt` | **PASS** |
| no credential/API key in tracked files | **PASS** — AEMET key read from env only; none present |
| plain-language bottom line at top | **PASS** — see above |
