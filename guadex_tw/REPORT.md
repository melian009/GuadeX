# River water-temperature model and 2045 projections for the Guadalquivir basin (GuadeX)

**Pipeline root:** `guadex_tw/` · **Reproduce with:** `python run_all.py`
**Status:** all deliverables produced from real downloads. Every number below carries its N.

----------------------------------------

## Missing Water temperature data

## Online Guadalquivir data 
** Historical data https://www.chguadalquivir.es/saih/DatosHistoricos_Avan.aspx
** BUT not temperature data 

## Link to the switch 
** https://drive.switch.ch/public.php/dav/files/rNd3V73S2ca6MkT/?accept=zip
** BUT not temperature data


## Bottom line (plain language) -- Calibration WTM from ATData

We calibrated a daily river water-temperature model from air temperature using **40,277
quality-controlled water-temperature observations at 980 Spanish monitoring sites** (46 of them in
the Guadalquivir), driven by NASA POWER air temperature. The national model transfers to the
Guadalquivir with a leave-one-basin-out **pooled NSE of 0.704** (held-out basin, 771 observations),
clearly beating the naive baselines (pooled NSE 0.26–0.41). Elevation matters: the
elevation×air-temperature interaction is significant under site-clustered inference and a mixed
model (national cluster-robust p = 7×10⁻¹²; MixedLM LRT p = 3×10⁻⁹) and improves out-of-sample
pooled NSE by **+0.011 nationally, +0.034 within the Guadalquivir, and +0.042 in the held-out-basin
test**. It is a *slope* change, not a level shift, and it is a refinement rather than a
transformation of the fit; the effect magnitude is smaller under a mixed model (−0.10 national,
−0.27 basin) than under OLS (−0.16, −0.30), and only the slope should be interpreted. Applying the
calibrated model to the AEMET/PNACC CMIP6 ensemble (11 GCMs × 4 SSPs, Guadalquivir sites, 5 km)
gives a **2036–2055 site-averaged warming of +0.71 °C (ssp126) to +0.94 °C (ssp585)** relative to
each model's own 1986–2005 baseline (across-GCM 10–90% spread 0.51–1.29 °C), rising to **+0.71 to
+1.23 °C by 2041–2070**. The projections are best trusted as a *relative climate signal with
quantified ensemble spread*, not as absolute temperatures for a specific reach: they rest on a
stationarity assumption, sparse (mostly monthly, 2024-dominated) historical observations, and a
predictor that differs between calibration (MERRA-2) and projection (AEMET/PNACC). Details and caveats
are in `LIMITATIONS.md`.

---

## 1. Executive summary

- **Calibration network (Stage B):** 39,944 paired daily observations at 980 Spanish sites; pooled
  LOSO NSE for the elevation-aware model **0.785** vs **0.773** without elevation (mean-of-station
  NSE 0.624 vs 0.591).
- **Local basin (Stage A):** 771 paired observations at 46 Guadalquivir sites; pooled LOSO NSE
  **0.762** (elevation) vs **0.729** (no elevation).
- **Held-out basin (decisive):** national model trained without the Guadalquivir predicts it at
  **NSE 0.704 / RMSE 2.92 °C** (771 obs) vs 0.662 without elevation.
- **Elevation answer:** the elevation×air-temperature interaction is significant under
  cluster-robust inference (Spain Wald p = 7×10⁻¹²; Guadalquivir p = 5×10⁻⁶) and under a MixedLM
  likelihood-ratio test (Spain χ²=34.9, p=3×10⁻⁹; Guadalquivir χ²=12.0, p=5×10⁻⁴) and is negative.
  The elevation *level* term is not stable across specifications; only the slope interaction is
  interpreted.
- **2045 window (2036–2055):** site-averaged water-temperature change **+0.71 / +0.80 / +0.83 /
  +0.94 °C** for ssp126/245/370/585, with across-GCM spread of roughly 0.5–1.3 °C.

---

## 2. Data inventory

| source | variable | N sites | N obs | period | resolution | license | access status |
|---|---|---:|---:|---|---|---|---|
| GRQA v1.3 `TEMP_GRQA.csv` (Spain) | river Tw | 204 | 23,502 | 1979–2018 | point | CC-BY-4.0 (DOI 10.5281/zenodo.7056647) | **SUCCESS** |
| GRQA v1.3 Guadalquivir extract | river Tw | 5 | 195 | 1990–2018 | point | CC-BY-4.0 | **SUCCESS** (matches verified ground truth) |
| Waterbase WISE6 (EEA SDI disaggregated) | river Tw | 2,607 | 27,203 | 2008–2024 | point | EEA open data | **SUCCESS** |
| EEA discodata SQL API | river Tw | 13,336 sites metadata | obs query timed out | – | – | EEA open data | **PARTIAL/FAILED** |
| Combined QC water temperature | river Tw | 980 | 40,277 | 1979–2024 | point | mixed | **SUCCESS** |
| NASA POWER (MERRA-2) daily T2M/MAX/MIN | air Ta | 980 | daily | 1981–2024 (+1981–2024 baseline for Guadalquivir sites) | 0.5°×0.625° | NASA open | **SUCCESS** |
| Open-Meteo ERA5 archive | air Ta | 33 sites usable | partial | – | ~9–25 km | CC-BY-4.0 | **FAILED/PARTIAL** (HTTP 429 weighted rate limit) |
| Copernicus DEM GLO-30 (AWS) | elevation | 992 | 1 per site (3×3 median) | static | 30 m | Copernicus | **SUCCESS** |
| Open-Meteo elevation (Copernicus DEM 90 m) | elevation | **9 usable** | 1 per site | static | 90 m | CC-BY-4.0 | **PARTIAL** (same rate limit) |
| AEMET/PNACC CMIP6 `tmean` (ESD-RegBA 5 km) | air Ta future | 11 GCM × 4 SSP × Guadalquivir sites | daily | 1950–2014 hist, 2015–2100 proj | 5 km | AEMET/PNACC public | **SUCCESS** (THREDDS NCSS broken; remote HDF5 subset) |
| CHG IDE `ggiscloud_root:a1748512575838_cuenca_guadalquivir` | basin polygon | – | – | – | vector | CHG | **PARTIAL** (WFS GetFeature HTTP 401; CHG ES050 water-body catchments used instead) |

Exact URLs, UTC timestamps, sizes, MD5/SHA256 and errors are in `PROVENANCE.md`; failed rows are
recorded with the literal HTTP error rather than replaced.

Note on elevation: Copernicus DEM GLO-30 (3×3 median) is the elevation source for all 992 sites,
but an **independent** second source (Open-Meteo DEM90) exists for only 9 sites, so the
dual-source cross-check is limited to those; a within-source check (3×3 median vs centre cell)
flags 0 sites above 50 m. The pipeline no longer describes the DEM as independently cross-checked
for every site, and the acceptance suite records the actual coverage.

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
(GLORICH 22,237; WATERBASE_DIRECT 16,834; GEMSTAT 1,206), span 1979-01-11 → 2024-12-30. Of these,
39,944 have a matched NASA POWER air temperature and enter the models (333 dropped for missing
POWER). `data/interim/water_temp_raw.csv` keeps the pre-QC extract for audit.

**Pipeline validation:** the Guadalquivir raw extract reproduces the verified ground truth exactly —
**195 rows / 5 sites** (ESP00016 = 67, ESP00017 = 58, ESP00018 = 68, plus two Waterbase single
observations); the naive bbox returns **406 rows / 11 sites**, and all **6 non-Guadalquivir sites are
excluded** (2 × ES040 Záncara/Valbuena, 3 Guadiana GEMSTAT, 1 Portuguese Guadiana).

---

## 4. Elevation

- **All 992 sites:** median 384 m, range **0–1,452 m** (Copernicus DEM GLO-30, 3×3 median; independent
  Open-Meteo DEM90 exists for 9 sites; 0 sites missing elevation).
- **Guadalquivir (46 sites):** median 170 m, range **0.5–1,127 m**.
- Site classification: river 925, reservoir 35, estuary 24, canal 8 (`data/processed/sites.csv`).
- The Guadalquivir range is wider than the original three-site subset (0–235 m), which is
  essential for identifying the elevation interaction within the basin. It rests on a small number
  of high headwater stations, so within-basin elevation inference is fragile.

---

## 5. Models

Predictor: `ta` = NASA POWER daily mean air temperature, **lapse-corrected to site elevation**
with `Ta_corr = Ta_grid + γ·(z_site − z_grid)`, **γ = −0.0065 K m⁻¹** (a stated parameter, not truth).
`z = elevation_m / 1000`.

- **Model 1 — per-station, per-month OLS** `Tw = α[m] + β[m]·Ta + ε` (fit where n ≥ 15, else pooled-month fallback).
- **Model 2 — pooled month effects + elevation interaction** (core, **area-free**):
  `Tw = a0 + Σ_m f(month) + b1·Ta + b2·z + b3·(Ta·z) + ε`.
  Month dummies vs a 2-harmonic seasonal term were both fitted; month dummies won on AIC in both
  scopes (Guadalquivir 3654.5 vs 3654.9; Spain 191,123 vs 191,137) and are retained.
- **Model 3 — lagged / thermal-inertia** `Tw = a0 + f(month) + Σ b_k·Ta_d[Lk] + b2·z`.
- **Model 4 — mixed effects** `Tw ~ Ta + z + Ta:z + month + (1 + Ta | site_id)` (BLUPs saved).
- **Model 5 — process-based (air2stream)** — **not fitted** (see omissions).
- **Basin-area sensitivity (not primary):** `Ta·log10(area)` is fitted on Spain only, because
  `upstream_basin_area` is ~55% populated nationally and 0% populated for the Guadalquivir. It is
  significant there (coefficient +0.0073 ± 0.0015, p ≈ 7×10⁻⁷) but is **not** part of the
  transferable model or the projections.

### Fitted coefficients (physical units, °C)

The basin-area term is excluded from the primary model so that a single formula transfers to and
projects the Guadalquivir. Standard errors are **cluster-robust by site** (observations are strongly
clustered: MixedLM ICC 0.444 nationally, 0.248 in the basin).

| term | Guadalquivir (n=771, 46 sites) | All Spain (n=39,944, 980 sites) |
|---|---|---|
| intercept | +6.58 ± 0.54 | +4.60 ± 0.07 |
| `b1` Ta (°C/°C) | **+0.436 ± 0.039** (cluster p≈8×10⁻²⁹) | **+0.510 ± 0.015** (cluster p≈3×10⁻²³⁹) |
| `b2` z (°C/km) | +1.77 ± 1.49 (p=0.23) | +0.066 ± 0.450 (p=0.88) |
| `b3` Ta·z (°C/°C/km) | **−0.299 ± 0.066** (p=5×10⁻⁶) | **−0.163 ± 0.024** (p=7×10⁻¹²) |
| residual SD | 2.54 °C | 2.65 °C |
| MixedLM ICC | 0.248 | 0.444 |
| MixedLM `b3` | −0.274 ± 0.063 | −0.103 ± 0.016 |

Significance of `b3` by three routes:

| test | Guadalquivir | Spain |
|---|---|---|
| naive OLS LRT | χ²=35.8, p=2×10⁻⁹ | χ²=519, p=6×10⁻¹¹⁵ |
| **cluster-robust Wald** | **t=−4.55, p=5×10⁻⁶** | **t=−6.85, p=7×10⁻¹²** |
| **MixedLM LRT** | **χ²=12.0, p=5×10⁻⁴** | **χ²=34.9, p=3×10⁻⁹** |

The naive OLS p-values are reported only for audit; the cluster-robust Wald and MixedLM LRT are the
defensible evidence. The cluster-robust national interval (±0.024) is more than three times wider
than the naive OLS interval (±0.007), and the MixedLM point estimate is smaller than OLS.

**Interpretation (physics):** air temperature already carries most of the elevation *level* signal
(colder air at altitude), so `b2` is small/insignificant once Ta is present; the MixedLM `z`
coefficient is unstable, confirming the level is not identifiable here. The significant negative
`b3` means **the water-to-air response slope flattens with elevation** — high reaches warm less per
degree of air warming — consistent with shorter growing seasons, snow/groundwater influence and
higher-gradient reaches.

**Sanity check:** predicted summer-minus-winter amplitude is 8.4–18.1 °C (Spain) and 7.3–15.2 °C
(Guadalquivir); **1** Guadalquivir site (a high, regulated reach) falls below the 8 °C rule of thumb
and is counted in `model_coefficients.json` (`n_sites_below_8c`).

### Model-selection rule and omissions

Final model chosen by **leave-one-station-out NSE**, preferring the simplest model within 0.02:
**Model 2 (month dummies + Ta + z + Ta:z)**. Model 4 (mixed effects) reaches similar skill
(Spain pooled NSE 0.780 population-level, 0.821 on 20 sampled held-out stations) but adds complexity
without a decisive gain, so it is reported as a sensitivity, not the primary model. The 20-station
true-LOSO subsample is deliberately small and seed-dependent; its pooled NSE is the number to trust.

- **Model 3:** only 3 of 980 sites meet the ≥24 obs/yr cadence rule, and all three fall within a
  single calendar year (2024), so a time-blocked CV of lag windows is impossible. Lag models are
  therefore fitted **in-sample only** and flagged (`logs/model3_summary.json`); lag structure is not
  identifiable from these data. All lag columns are still provided in the paired table, computed
  (never interpolated) with full-window requirements.
- **Model 5 (air2stream):** **not fitted.** `air2stream` is not on PyPI and an 8-parameter process
  model cannot be identified from mostly monthly grab samples with no flow/shading data.
  Documented omission.

---

## 6. Validation (Section 9)

Metrics are computed on held-out data only. **POOLED_HELDOUT** is the primary metric: a single NSE
over the concatenated held-out predictions. **MEAN_STATION_NSE** (the unweighted mean of per-station
NSE) is reported as a secondary, station-equal-weight metric; the two are not interchangeable, and
only the pooled comparison against baselines is meaningful. Full table: `outputs/tables/cv_metrics.csv`.

### 6.1 Leave-one-station-out (LOSO), pooled over folds — primary

| model | scope | N test | NSE | RMSE (°C) | MAE (°C) | bias (°C) |
|---|---|---:|---:|---:|---:|---:|
| Tw = Ta | Spain | 39,944 | 0.583 | 3.69 | 2.76 | +0.44 |
| Tw = Ta + constant | Spain | 39,944 | 0.589 | 3.67 | 2.77 | −0.00 |
| monthly Tw climatology | Spain | 39,944 | 0.657 | 3.35 | 2.54 | +0.00 |
| Model 1 (per station-month) | Spain | 39,944 | 0.817 | 2.45 | 1.69 | −0.09 |
| **Model 2 (no elevation)** | Spain | 39,944 | 0.773 | 2.72 | 1.97 | +0.00 |
| **Model 2 (Ta + z + Ta:z)** | Spain | 39,944 | **0.785** | **2.65** | 1.90 | +0.00 |
| Model 4 (MixedLM, population) | Spain | 39,944 | 0.780 | 2.68 | 1.93 | +0.24 |
| Model 4 (MixedLM, true LOSO) | Spain | 1,042 | 0.821 | 2.52 | 1.96 | +0.25 |
| Tw = Ta | Guadalquivir | 771 | 0.383 | 4.22 | 3.17 | +1.24 |
| Tw = Ta + constant | Guadalquivir | 771 | 0.431 | 4.06 | 3.13 | −0.02 |
| monthly Tw climatology | Guadalquivir | 771 | 0.669 | 3.10 | 2.46 | −0.01 |
| Model 1 (per station-month) | Guadalquivir | 771 | 0.747 | 2.71 | 2.10 | −0.00 |
| Model 2 (no elevation) | Guadalquivir | 771 | 0.729 | 2.80 | 2.18 | −0.02 |
| **Model 2 (Ta + z + Ta:z)** | Guadalquivir | 771 | **0.762** | **2.62** | 2.01 | −0.03 |
| Model 4 (MixedLM, population) | Guadalquivir | 771 | 0.775 | 2.55 | 1.95 | −0.15 |

Secondary `MEAN_STATION_NSE` (not comparable to the baselines above): Spain Model 2 (elev) 0.624 vs
0.591 without elevation; Guadalquivir 0.609 vs 0.526. Model 1 is in-sample and therefore optimistic.

At monthly cadence the pooled `Tw = Ta` baseline is mediocre (0.38–0.58), and the monthly
climatology is a strong baseline (0.66–0.67), which the models beat. The *mean-of-station* metric
previously used in this report made `Tw = Ta` look worse than useless (negative NSE) in the
Guadalquivir; that was an artifact of equal-weighting 12-observation stations, not a property of the
baseline.

### 6.2 Leave-one-basin-out (decisive Stage-B transfer test)

Train on all non-Guadalquivir Spanish sites (39,173 obs / 934 sites), predict the 46 Guadalquivir
sites (771 obs), pooled metric:

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
`fig01` Tw–Ta scatter (single series, lines labelled inline — no per-site legend); `fig02` residuals
vs month and vs elevation; `fig03` LOSO predicted vs observed (Stage B only, N=39,944);
`fig04` elevation-coloured site map with the ES050 basin; `fig05` air/water time series for the three
Guadalquivir GEMSTAT sites; `fig06` future monthly ensemble spread; `fig07` distribution of
2036–2055 ΔTw by SSP; `fig08` the elevation answer — per-station slope and offset vs elevation with
95% CI.

The historical table `water_temp_historical_sites.csv` now carries an **observation-specific** 90%
OLS prediction interval (leverage-aware) instead of a single constant ±1.645·SD interval; it is
conditional on the fixed-effects model and does not add the between-site random effect.

---

## 7. Future projections to 2045 (Section 10)

- **Predictor consistency:** the model is applied to AEMET/PNACC 5 km `tmean` after a **per-site,
  per-GCM additive bias correction** over **1986–2005** that aligns PNACC to the NASA POWER
  calibration climatology. Because the model is linear in Ta and the correction is an additive
  constant, it **cancels exactly in the anomalies** (future mean minus each model's own historical
  mean); it corrects absolute temperatures only, and the relative signal is driven solely by each
  GCM's own warming. Anomalies are computed **per model against that model's own historical run**,
  never against a multi-model historical mean.
- **Ensemble:** all **11 CMIP6 GCMs** for **ssp126, ssp245, ssp370, ssp585**.
- **Baseline period:** 1986–2005.
- **Windows:** primary **2036–2055** (centred on 2045), standard **2041–2070**, appendix **2021–2040**.
- **Sites:** the shipped run uses the **6 Guadalquivir sites with complete 11-GCM coverage**, spanning
  **0.5–537 m**. The extraction has been re-stratified to span **0.5–1,127 m** (including the two
  highest basin reaches at 850 m and 1,127 m); the full re-extraction is gated on the AEMET/PNACC
  endpoint and has not completed, so the shipped numbers do not yet include the high-elevation
  reaches. `logs/pnacc_extract_summary.json` records the corrected selection.

### Ensemble ΔTw (°C, per-model anomaly vs own 1986–2005 baseline)

Two spreads are reported separately so that between-site differences are not presented as model
uncertainty. `median` is the median across sites of each site's 11-GCM median; `GCM 10–90` is the
median across sites of each site's across-GCM 10th–90th range; `site 10–90` is the spread of the
site medians.

| window | statistic | ssp126 | ssp245 | ssp370 | ssp585 |
|---|---|---|---|---|---|
| 2021–2040 | median | 0.49 | 0.52 | 0.55 | 0.60 |
| | GCM 10–90 | 0.34–0.82 | 0.38–0.85 | 0.43–0.83 | 0.45–0.89 |
| **2036–2055** | **median** | **0.71** | **0.80** | **0.83** | **0.94** |
| | **GCM 10–90** | **0.51–1.07** | **0.56–1.27** | **0.66–1.20** | **0.71–1.29** |
| | site 10–90 | 0.69–0.74 | 0.77–0.81 | 0.80–0.84 | 0.91–0.96 |
| 2041–2070 | median | 0.71 | 0.92 | 1.05 | 1.23 |
| | GCM 10–90 | 0.54–1.14 | 0.71–1.43 | 0.79–1.42 | 0.91–1.59 |

N = 6 sites × 11 GCMs per scenario/period. Full per-GCM values are in
`outputs/tables/water_temp_future_2045.csv`; the site-averaged summary is in
`outputs/tables/water_temp_future_ensemble_summary.csv`.

### Elevation dependence of the projected change

Because the model is linear in Ta with slope `(b1 + b3·z)`, the elevation dependence of the change is
`d(ΔTw)/dz = b3 · ΔTa`. For 2036–2055 this implies **−0.24 °C per km (ssp126) to −0.33 °C per km
(ssp585)**: a 1,000 m reach is projected to warm about a third of a degree less than sea level for
the same air-temperature signal. This is a model-based deduction, not an independent spatial
simulation, and is recorded per period in `logs/projection_summary.json` and at every site in
`outputs/tables/water_temp_future_elevation_profile.csv`.

### Absolute values and thermal metrics

Applied to absolute calibrated temperatures, the 2036–2055 ensemble-median annual mean water
temperature at the 6 shipped sites is ~15.0–18.2 °C, with summer (JJA) means ~21.7–24.7 °C. Per site,
scenario and period, `outputs/tables/thermal_metrics.csv` reports
`mean_annual_c, mean_summer_c, mean_winter_c, p95_annual_c, summer_p95_c,
n_days_above_20/25/28, longest_run_above_28, n_days_below_10, temp_range_c, temp_sd_c`.
**Thresholds (20/25/28 °C hot, 10 °C cold) are placeholders defined at the top of
`scripts/10_project.py` and must be swapped for the actual GuadeX species optima.**

Absolute values inherit the calibration's predictors and the per-(site, GCM) bias correction; the
*changes* do not, and are the more trustworthy output.

---

## 8. GuadeX-ready outputs

| file | contents |
|---|---|
| `outputs/tables/water_temp_historical_sites.csv` | per-observation modelled vs observed Tw, Ta, observation-specific 90% prediction interval (N=39,944) |
| `outputs/tables/water_temp_future_2045.csv` | per GCM × scenario × site × period Tw and ΔTw + per-site ensemble statistics |
| `outputs/tables/water_temp_future_ensemble_summary.csv` | site-averaged GCM median and spread, separated from between-site spread |
| `outputs/tables/water_temp_future_elevation_profile.csv` | ΔTw vs site elevation, per site/scenario/period |
| `outputs/tables/thermal_metrics.csv` | per site × scenario × period thermal metrics (72 rows) |
| `data/processed/sites.csv` | site table with elevation, cross-check diagnostics, basin flag, class, cadence |
| `data/processed/paired_observations.csv` | paired Tw–Ta calibration table with lag columns (N=40,277) |
| `data/processed/paired_sites_matched.csv` | matched POWER/ROCIO/ERA5 cells and distances per site |

---

## 9. Assumptions and limitations (summary)

1. **Stationarity** of the calibrated Tw–Ta relation under climate change is assumed; changes in
   flow, reservoir operation, shading and abstraction are ignored. This is the largest structural
   uncertainty.
2. **Predictor mismatch** between calibration (NASA POWER/MERRA-2) and projection (AEMET/PNACC),
   addressed for absolute values by per-site/per-GCM bias correction; the correction cancels in the
   reported anomalies, so any POWER-vs-PNACC difference in *warming trend* is an unaddressed, shared bias.
3. **Sparse, mostly monthly observations**, and the Guadalquivir record is 2024-dominated (~2/3 of
   the 771 modelled observations); the three multi-year Guadalquivir series end in 1995.
4. **Reservoir/estuary sites** are pooled with rivers (classified but not separately modelled).
5. **Single GCM member per model** (no internal-variability ensemble).
6. **Site-clustered inference**: cluster-robust and MixedLM tests are used because naive OLS
   p-values are far too small; the elevation *level* term is not stable across specifications.
7. **Projection elevation coverage**: the shipped run covers 0.5–537 m; the re-stratified selection
   spans 0.5–1,127 m but the full re-extraction is pending on the PNACC endpoint.

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
| elevation for every retained site, cross-check coverage reported honestly | **PASS** — 0 missing; independent cross-check on 9 of 992 sites (recorded) |
| `paired_*.csv` with no fabricated lags; cadence rule applied | **PASS** — 40,277 rows; 3 eligible sites, all single-year, flagged |
| all five models (or documented omissions) | **PASS** — 1,2,3 (in-sample, flagged),4 fitted; 5 omitted with reason |
| pooled LOSO NSE for every model + three naive baselines | **PASS** — `cv_metrics.csv` (`POOLED_HELDOUT`) |
| elevation interaction tested cluster-robust **and** by MixedLM LRT **and** CV | **PASS** |
| primary model area-free and identical to the projection model | **PASS** — `tw_obs ~ ta + C(month) + z + ta_z` |
| future projections ≥3 SSPs, full ensemble, per-model anomalies | **PASS** — 4 SSPs, 11 GCMs, 6 complete sites; PNACC selection re-stratified to 0.5–1,127 m |
| ensemble spread separated from between-site spread | **PASS** — `water_temp_future_ensemble_summary.csv` |
| every figure has labels, units, N, caption | **PASS** — `captions.txt`; fig01 has no per-site legend |
| `REPORT.md`, `LIMITATIONS.md`, `PROVENANCE.md`, `run_all.py`, `requirements.txt` | **PASS** |
| no credential literal in tracked files | **PASS** — AEMET key and EEA WebDAV token read from env only; acceptance scans for literals |
| plain-language bottom line at top | **PASS** — see above |
