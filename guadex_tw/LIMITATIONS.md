# LIMITATIONS

This document states plainly what the GuadeX water-temperature analysis can and
cannot support. Numbers are drawn from `logs/qc_log.json`,
`logs/waterbase_direct_qc.json`, `logs/site_summary.json` and
`outputs/tables/cv_metrics.csv`.

## 1. How much data actually drives the model

| stage | sites | paired daily obs | notes |
|---|---:|---:|---|
| Guadalquivir (Stage A) | 46 | 771 | 3 GEMSTAT sites (1990–1995, monthly) + 43 Waterbase sites, but ~2/3 of the 771 modelled observations are from **2024 alone** (2009–2010 contribute ~90, 2018 two, and nothing in 1996–2008) |
| All Spain (Stage B) | 980 | 39,944 | GLORICH 22,237 / Waterbase-direct 16,834 / GEMSTAT 1,206 |

The national model is driven by ~40k observations across 980 sites, but the
distribution is extremely uneven: many Waterbase sites contribute the 12–20
observations needed to clear the ≥12 rule and little else, while a handful of
GLORICH/GEMSTAT sites carry dense multi-year records. The **effective** sample
size for the elevation interaction is closer to the number of *sites spanning
elevation* than to 40,000.

## 2. Temporal and spatial gaps

- **Guadalquivir historical record is thin, old and front-loaded.** The only
  true multi-year Guadalquivir river series in GRQA are the three GEMSTAT sites
  (Menjíbar, Peñaflor, Seville), roughly monthly, 1990–1995. The modern record
  is Waterbase monitoring dominated by a single year: of the 771 modelled
  Guadalquivir observations, ~512 are from 2024, ~193 from 1990–1995, and only
  a handful from 2009–2018. The basin month effects and within-basin elevation
  slope therefore rest largely on a one-year snapshot, which is a real risk
  given the −0.52 °C recent-period bias seen in the blocked time-split.
- **No continuous daily water-temperature series exists for free-flowing
  Guadalquivir reaches.** Every site is a grab sample. This is why no credible
  daily lag model can be fitted (Section 6 of the report): the three sites that
  meet the ≥24 obs/yr cadence rule all fall inside a *single calendar year*
  (2024), so a time-blocked cross-validation of lag windows is impossible.
- **Spatial gap.** Site locations are monitoring stations, not a representative
  sample of the ~1,000 GuadeX reaches. Large lowland canalised reaches and
  regulated main-stem sections are poorly covered.
- **Seasonal gap.** Winter sampling is thin; the cold tail (n_days_below_10) is
  model extrapolation more than observation.

## 3. Sources that failed or were not used

| source | status | reason |
|---|---|---|
| Open-Meteo ERA5 archive | **FAILED/PARTIAL** | weighted rate limit returned HTTP 429 for bulk site requests; cached responses retained as a cross-check only |
| CHG IDE WFS basin polygon | **FAILED** | `GetFeature` returned HTTP 401; WMS capabilities available but no geometry. Used the official CHG ES050 water-body catchment layer already in the repo + `ES050` site-id prefix |
| AEMET OpenData stations | **NOT USED** | requires an API key; the no-auth path (NASA POWER) was sufficient. Columns `nearest_aemet_*` are therefore blank |
| ERA5-Land, E-OBS (CDS) | **NOT USED** | Copernicus CDS account + licence acceptance required; no-auth path sufficient |
| ESA CCI Lakes LSWT (supplementary) | **NOT USED** | optional Section 4.6 evidence; not needed for the river calibration and not fetched |
| CHG SAIH `TEMPERATURA` | **EXCLUDED BY DESIGN** | confirmed unusable as a water-temperature record |

The CHG WFS failure and Open-Meteo throttling are recorded as FAILED rows in
`PROVENANCE.md` with the exact URLs and errors.

## 4. Why the Guadalquivir-only elevation range is insufficient

The original Guadalquivir calibration universe was **three low, low-gradient
sites (0–235 m)**. Within that range the elevation term is nearly collinear with
the intercept and the slope cannot be identified. The national Stage-B model is
what identifies elevation: Spanish sites span 0–1,452 m versus 0.5–1,127 m for
the (now enlarged) Guadalquivir set. Even so, the Guadalquivir sites are
concentrated below ~450 m, so elevation dependence *within the basin* rests on a
few high Genil/Guadiana-Menor headwater stations. The projection sample is now explicitly stratified to include the two highest basin sites (850 m and 1,127 m) and two of the lowest, so the projected set spans 0.5–1,127 m rather than the 0.5–537 m of the earlier unstratified selection.

### 4.1 Statistical inference under site clustering

Observations are strongly clustered within monitoring sites (MixedLM ICC ≈ 0.44 nationally, 0.25 in the basin). Treating the ~40,000 observations as independent makes ordinary-OLS standard errors far too small. Model 2 is therefore reported with **cluster-robust standard errors by site**, and the elevation interaction is additionally tested with a MixedLM likelihood-ratio test. The effect survives (basin cluster p = 5×10⁻⁶; national cluster p = 7×10⁻¹²), but the honest confidence interval on the national slope is roughly **three times wider** than the naive OLS interval, and the MixedLM point estimate is smaller (`ta:z ≈ −0.10` national, `−0.27` basin) than the OLS estimate (`−0.16`, `−0.30`). The elevation *level* term (`z`) is not stable across specifications (OLS p ≈ 0.56–0.88; MixedLM strongly negative), so only the slope interaction should be interpreted.

### 4.2 Validation metrics

Two validation summaries are reported and they are not interchangeable:
**POOLED_HELDOUT** is the true pooled NSE over concatenated held-out predictions (the primary number); **MEAN_STATION_NSE** is the unweighted mean of per-station NSE and gives a station with 12 observations the same weight as one with 400. Both are written to `cv_metrics.csv`. The naive baselines look far stronger under the pooled metric (e.g. national `Tw = Ta` NSE 0.58 pooled vs 0.15 as a station mean), so only the pooled comparison against baselines is meaningful.

## 5. Modelling weaknesses

- **Predictor mismatch.** Calibration uses NASA POWER (MERRA-2, 0.5°) daily T2M
  lapse-corrected to site elevation with γ = −0.0065 K m⁻¹. The projection uses
  AEMET/PNACC 5 km `tmean` tuned to ROCIO-IBEB. We align them with a per-site,
  per-GCM additive bias correction over 1986–2005. Because the model is linear
  in Ta and the correction is an additive constant per (site, GCM), **it cancels
  exactly in the reported anomalies** (future mean minus each model's own
  baseline); it corrects absolute temperatures only. The *relative* signal is
  therefore driven solely by each GCM's own warming, and any difference between
  the POWER and PNACC warming trends is an unaddressed, shared bias. Residual
  definitional differences (24-h mean vs downscaled bias-adjusted mean) remain.
- **Air-temperature product choice was forced.** Open-Meteo ERA5 was the intended
  workhorse; the rate limit made NASA POWER the only consistent no-auth source.
- **Monthly aggregation.** The Stage-A record is monthly, so "summer" for those
  sites is two or three samples.
- **Basin area excluded from the primary model.** `upstream_basin_area` is only
  ~55% populated nationally and 0% populated for the Guadalquivir sites, so it
  cannot enter a model that is transferred to the basin. The area interaction is
  fitted only as a Spain sensitivity (`ta·log10(area)` is significant there,
  p ≈ 7×10⁻⁷) and is *not* used for projections; including it moves the national
  `ta` coefficient from 0.51 to 0.44 and is not transferable.
- **Stationarity assumption.** The projection assumes the calibrated Tw~Ta
  relation is stationary. It ignores changes in flow regime, reservoir
  operation, riparian shading, abstraction and land use — all of which can shift
  the water-to-air offset. This is the single largest source of structural
  uncertainty in the 2045 numbers.
- **Reservoir pooling.** 35 sites are classified `reservoir` and 24 `estuary`;
  they are retained in the pooled fit with a `site_class` field but not modelled
  separately, because per-class data are too sparse. This can bias the intercept.
- **No process model.** Air2stream (Model 5) was not implemented; an 8-parameter
  process model cannot be identified from monthly grab samples. No flow, shading
  or channel-geometry predictors were available.

## 6. What would materially improve the result (minimum additional data)

1. **A daily water-temperature logger network on ~15–25 Guadalquivir reaches
   spanning 0–1,300 m for ≥5 years.** This alone would move the basin model from
   inference-by-transfer to direct calibration and make lag/inertia terms
   identifiable.
2. **Recent (post-2014) co-located air temperature at the water sites** (AEMET
   ROCIO 5 km is the natural choice) so calibration and projection share a
   predictor definition.
3. **Reservoir-vs-river stratification** with discharge or residence-time
   covariates; reservoir thermal behaviour is not captured by Ta alone.
4. **A curated Guadalquivir basin polygon via authenticated CHG WFS**, replacing
   the ES050 water-body-catchment approximation.
5. Optionally, **ESA CCI Lakes LSWT** over the Guadalquivir reservoirs as a
   dense multi-elevation shape check (not a river calibration target).

Until (1) and (2) exist, the 2045 projections should be used as a **relative
climate-change signal** (how much warmer, with what ensemble spread), not as
absolute forecasts of a given reach.
