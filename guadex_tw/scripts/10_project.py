"""Project river water temperature to 2045 with the calibrated model (Sections 10-11).

Pipeline:
  1. Fit the Stage-B (all-Spain) elevation-aware Model 2 on the paired table.
  2. Bias-correct PNACC tmean per (site, GCM) to the NASA POWER calibration
     climatology over 1986-2005, so the projection predictor matches the
     predictor the model was calibrated on.
  3. Apply the model to each GCM/scenario, compute per-model anomalies against
     that model's own historical run, and report the 11-model ensemble.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as CM

# ---- tunable thresholds (placeholders for GuadeX species optima) ----
THRESHOLDS = {"hot_c": [20, 25, 28], "cold_c": 10}
# ---------------------------------------------------------------------

GAMMA = -0.0065
GCMS = ["ACCESS-CM2", "CMCC-CM2-SR5", "CNRM-ESM2-1", "EC-Earth3-Veg", "IITM-ESM",
        "KACE-1-0-G", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NorESM2-MM",
        "UKESM1-0-LL"]
SCENARIOS = ["ssp126", "ssp245", "ssp370", "ssp585"]
BASELINE = (1986, 2005)
PERIODS = {"2036-2055": (2036, 2055), "2041-2070": (2041, 2070), "2021-2040": (2021, 2040)}
POWER_MIN = pd.Timestamp("1981-01-01")


def fit_model() -> tuple:
    paired = pd.read_csv(CM.PROCESSED / "paired_observations.csv", dtype={"site_id": str})
    d = paired.dropna(subset=["tw_obs", "ta_mean_corr", "elevation_m", "month"]).copy()
    d["ta"] = d["ta_mean_corr"]
    d["z"] = d["elevation_m"] / 1000.0
    d["ta_z"] = d["ta"] * d["z"]
    d["month"] = d["month"].astype(int)
    m = smf.ols("tw_obs ~ ta + C(month) + z + ta_z", data=d).fit()
    return m, d


def predict_tw(model, ta: np.ndarray, z: np.ndarray, month: np.ndarray) -> np.ndarray:
    p = float(model.params.get("Intercept", 0.0))
    p = p + float(model.params.get("ta", 1.0)) * np.asarray(ta)
    p = p + float(model.params.get("z", 0.0)) * np.asarray(z)
    p = p + float(model.params.get("ta_z", 0.0)) * np.asarray(ta) * np.asarray(z)
    for k, v in model.params.items():
        if k.startswith("C(month)[T."):
            mm = int(k.split("T.")[1].rstrip("]"))
            p = p + float(v) * (np.asarray(month) == mm)
    return p


def power_full(sid: str, lat: float, lon: float) -> dict:
    """POWER 1981-2024 series for the historical baseline (cached)."""
    import requests
    path = CM.RAW / "nasa_power" / f"{sid}_full.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    params = {"parameters": "T2M,T2M_MAX,T2M_MIN", "community": "AG",
              "longitude": lon, "latitude": lat, "start": "19810101", "end": "20241231",
              "format": "JSON", "time-standard": "UTC"}
    r = requests.get("https://power.larc.nasa.gov/api/temporal/daily/point",
                     params=params, timeout=180, headers={"User-Agent": CM.USER_AGENT})
    r.raise_for_status()
    j = r.json()
    if "properties" not in j:
        raise RuntimeError(f"POWER error: {j.get('messages')}")
    j["_retrieved_utc"] = CM.utcnow_iso()
    path.write_text(json.dumps(j), encoding="utf-8")
    return j


def power_series_df(j: dict, z_site: float) -> pd.DataFrame:
    params = j["properties"]["parameter"]
    coords = j.get("geometry", {}).get("coordinates") or [None, None, None]
    dates = sorted(params["T2M"].keys())
    df = pd.DataFrame({"date": pd.to_datetime(dates, format="%Y%m%d"),
                       "ta": [params["T2M"].get(d, np.nan) for d in dates]})
    df["ta"] = pd.to_numeric(df["ta"], errors="coerce").replace(-999.0, np.nan)
    z_cell = coords[2]
    df["ta"] = df["ta"] + (GAMMA * (z_site - z_cell) if np.isfinite(z_cell) else 0.0)
    return df


def stat_metrics(tw: pd.Series, dates: pd.DatetimeIndex) -> dict:
    tw = pd.Series(np.asarray(tw, float), index=dates)
    s = pd.DataFrame({"tw": tw})
    s["month"] = s.index.month
    summer = s[s["month"].isin([6, 7, 8])]["tw"]
    winter = s[s["month"].isin([12, 1, 2])]["tw"]
    hot = THRESHOLDS["hot_c"]
    out = {
        "mean_annual_c": float(s["tw"].mean()),
        "mean_summer_c": float(summer.mean()),
        "mean_winter_c": float(winter.mean()),
        "p95_annual_c": float(s["tw"].quantile(0.95)),
        "p05_annual_c": float(s["tw"].quantile(0.05)),
        "summer_p95_c": float(summer.quantile(0.95)),
        "temp_range_c": float(s["tw"].max() - s["tw"].min()),
        "temp_sd_c": float(s["tw"].std()),
    }
    for th in hot:
        out[f"n_days_above_{th}"] = int((s["tw"] > th).sum())
    out[f"n_days_below_{THRESHOLDS['cold_c']}"] = int((s["tw"] < THRESHOLDS["cold_c"]).sum())
    # longest run above hottest threshold
    hott = max(hot)
    above = (s["tw"] > hott).to_numpy()
    longest = cur = 0
    for a in above:
        cur = cur + 1 if a else 0
        longest = max(longest, cur)
    out[f"longest_run_above_{hott}"] = int(longest)
    out["tw_amplitude"] = float(summer.mean() - winter.mean())
    return out


def main() -> None:
    model, d = fit_model()
    print("Stage-B Model2 fitted on", len(d), "obs /", d["site_id"].nunique(), "sites")
    (CM.MODELS / "projection_model_coefficients.json").write_text(
        json.dumps({"params": {k: float(v) for k, v in model.params.items()},
                    "bse": {k: float(v) for k, v in model.bse.items()},
                    "n": int(len(d)), "n_sites": int(d["site_id"].nunique())}, indent=2),
        encoding="utf-8")

    # --- sites to project: Guadalquivir ---
    sites = pd.read_csv(CM.PROCESSED / "sites.csv", dtype={"site_id": str})
    gq = sites[sites["in_guadalquivir_basin"]].dropna(subset=["elevation_m"])
    print("projecting at", len(gq), "Guadalquivir sites")

    # POWER 1986-2005 climatology per site
    base_lo, base_hi = BASELINE
    power_clim = {}
    for _, s in gq.iterrows():
        try:
            j = power_full(str(s["site_id"]), float(s["lat"]), float(s["lon"]))
            df = power_series_df(j, float(s["elevation_m"]))
            m = df[(df["date"].dt.year >= base_lo) & (df["date"].dt.year <= base_hi)]
            power_clim[s["site_id"]] = float(m["ta"].mean())
        except Exception as exc:  # noqa: BLE001
            print(f"  POWER full FAILED for {s['site_id']}: {exc!r}")
            CM.record(f"NASA POWER full baseline {s['site_id']}",
                     url="https://power.larc.nasa.gov/api/temporal/daily/point",
                     status="FAILED", error=repr(exc))
    print("POWER baseline climatology sites:", len(power_clim))

    pnacc = pd.read_parquet(CM.INTERIM / "pnacc_tmean_guadalquivir.parquet")
    pnacc["date"] = pd.to_datetime(pnacc["date"])
    gq = gq[gq["site_id"].isin(pnacc["site_id"].unique())]
    print("sites with PNACC extraction:", len(gq))
    pnacc["month"] = pnacc["date"].dt.month
    pnacc["year"] = pnacc["date"].dt.year

    zmap = dict(zip(gq["site_id"], gq["elevation_m"] / 1000.0))
    rows_future, metrics_rows, hist_rows = [], [], []

    for sid in gq["site_id"]:
        sub = pnacc[pnacc["site_id"] == sid]
        z = zmap[sid]
        pc = power_clim.get(sid, np.nan)
        # per-GCM historical bias correction (1986-2005)
        deltas = {}
        for gcm in GCMS:
            h = sub[(sub["gcm"] == gcm) & (sub["experiment"] == "historical")]
            hb = h[(h["year"] >= BASELINE[0]) & (h["year"] <= BASELINE[1])]
            deltas[gcm] = (pc - float(hb["tmean"].mean())) if len(hb) and np.isfinite(pc) else 0.0
        # historical modelled Tw for the site (all models, baseline)
        for gcm in GCMS:
            h = sub[(sub["gcm"] == gcm) & (sub["experiment"] == "historical")].copy()
            if h.empty:
                continue
            h["ta_adj"] = h["tmean"] + deltas[gcm]
            h["tw"] = predict_tw(model, h["ta_adj"].to_numpy(), z, h["month"].to_numpy())
            hb = h[(h["year"] >= BASELINE[0]) & (h["year"] <= BASELINE[1])]
            base_mean = float(hb["tw"].mean())
            hist_rows.append(pd.DataFrame({"site_id": sid, "gcm": gcm, "date": h["date"],
                                           "tw_modelled": h["tw"]}))
            for scen in SCENARIOS:
                f = sub[(sub["gcm"] == gcm) & (sub["experiment"] == scen)].copy()
                if f.empty:
                    continue
                f["ta_adj"] = f["tmean"] + deltas[gcm]
                f["tw"] = predict_tw(model, f["ta_adj"].to_numpy(), z, f["month"].to_numpy())
                for pname, (y0, y1) in PERIODS.items():
                    fp = f[(f["year"] >= y0) & (f["year"] <= y1)]
                    if fp.empty:
                        continue
                    mt = stat_metrics(fp["tw"].to_numpy(), pd.DatetimeIndex(fp["date"]))
                    rows_future.append({
                        "site_id": sid, "scenario": scen, "gcm": gcm,
                        "period_start": y0, "period_end": y1,
                        "tw_mean": mt["mean_annual_c"], "tw_summer_mean": mt["mean_summer_c"],
                        "tw_winter_mean": mt["mean_winter_c"], "tw_p95": mt["p95_annual_c"],
                        "tw_p05": mt["p05_annual_c"], "tw_amplitude": mt["tw_amplitude"],
                        "delta_tw_mean": mt["mean_annual_c"] - base_mean,
                    })
                    # ensemble median daily series across GCMs handled after loop
        # store per-site/scenario daily Tw for ensemble metrics
        # (kept in a separate parquet below)

    future = pd.DataFrame(rows_future)
    # ensemble spread across GCMs per site/scenario/period
    grp = future.groupby(["site_id", "scenario", "period_start", "period_end"])
    ens = grp["tw_mean"].agg(
        ensemble_median="median", ensemble_p25=lambda s: s.quantile(0.25),
        ensemble_p75=lambda s: s.quantile(0.75), ensemble_p10=lambda s: s.quantile(0.10),
        ensemble_p90=lambda s: s.quantile(0.90)).reset_index()
    future = future.merge(ens, on=["site_id", "scenario", "period_start", "period_end"], how="left")
    future.to_csv(CM.TABLES / "water_temp_future_2045.csv", index=False)
    print(f"wrote water_temp_future_2045.csv ({len(future)} rows)")

    # --- ensemble-median daily series -> thermal metrics ---
    daily_rows = []
    for sid in gq["site_id"]:
        sub = pnacc[pnacc["site_id"] == sid]
        z = zmap[sid]
        for scen in SCENARIOS:
            all_days = []
            # baseline delta per gcm
            deltas = {}
            for gcm in GCMS:
                hb = sub[(sub["gcm"] == gcm) & (sub["experiment"] == "historical")]
                hb = hb[(hb["year"] >= BASELINE[0]) & (hb["year"] <= BASELINE[1])]
                pc = power_clim.get(sid, np.nan)
                deltas[gcm] = (pc - float(hb["tmean"].mean())) if len(hb) and np.isfinite(pc) else 0.0
            for gcm in GCMS:
                f = sub[(sub["gcm"] == gcm) & (sub["experiment"] == scen)].copy()
                if f.empty:
                    continue
                f["ta_adj"] = f["tmean"] + deltas[gcm]
                f["tw"] = predict_tw(model, f["ta_adj"].to_numpy(), z, f["month"].to_numpy())
                all_days.append(f[["date", "tw"]].assign(gcm=gcm))
            if not all_days:
                continue
            mat = pd.concat(all_days).pivot_table(index="date", columns="gcm", values="tw")
            daily_rows.append(pd.DataFrame({"site_id": sid, "scenario": scen,
                                            "date": mat.index,
                                            "tw_ensemble_median": mat.median(axis=1).values}))
    daily = pd.concat(daily_rows, ignore_index=True) if daily_rows else pd.DataFrame()
    daily.to_parquet(CM.INTERIM / "tw_future_ensemble_daily.parquet", index=False)

    mrows = []
    for (sid, scen), g in daily.groupby(["site_id", "scenario"]):
        for pname, (y0, y1) in PERIODS.items():
            gp = g[(g["date"].dt.year >= y0) & (g["date"].dt.year <= y1)]
            if gp.empty:
                continue
            mt = stat_metrics(gp["tw_ensemble_median"].to_numpy(), pd.DatetimeIndex(gp["date"]))
            mrows.append({"site_id": sid, "scenario": scen, "period": pname, **mt})
    thermal = pd.DataFrame(mrows)
    thermal.insert(2, "period_start", thermal["period"].map(lambda p: PERIODS[p][0]))
    thermal.insert(3, "period_end", thermal["period"].map(lambda p: PERIODS[p][1]))
    thermal.to_csv(CM.TABLES / "thermal_metrics.csv", index=False)
    print(f"wrote thermal_metrics.csv ({len(thermal)} rows)")

    # historical modelled vs observed (Guadalquivir)
    hist = pd.concat(hist_rows, ignore_index=True) if hist_rows else pd.DataFrame()
    obs = d[d["in_guadalquivir_basin"]][["site_id", "obs_date", "tw_obs", "ta_mean_corr", "source"]]
    obs = obs.rename(columns={"obs_date": "date"})
    if not hist.empty:
        hist.to_csv(CM.INTERIM / "tw_historical_modelled.csv", index=False)

    CM.record("Generated future Tw projections", status="SUCCESS",
             path=str(CM.TABLES / "water_temp_future_2045.csv"), checksum=False,
             notes=f"{len(GCMS)} GCMs x {len(SCENARIOS)} SSPs; periods {list(PERIODS)}; "
                   f"baseline {BASELINE}; anomalies per model vs own historical run")
    summary = {"n_gq_sites": int(len(gq)), "n_gcm": len(GCMS), "scenarios": SCENARIOS,
               "periods": PERIODS, "baseline": BASELINE, "thresholds": THRESHOLDS,
               "power_baseline_sites": len(power_clim)}
    (CM.LOGS / "projection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
