"""Produce GuadeX-ready tables and all figures (Sections 11).

Every figure states N in its caption/title and saves a caption to
outputs/figures/captions.txt.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as CM

BASIN_SHP = CM.ROOT.parent / "data" / "Version_02-01-2026-Masas" / "Cuencas_masas_agua_4c.shp"
DPI = 200
CAPTIONS: list[str] = []


def save(fig, name: str, caption: str) -> None:
    cap = f"[{name}] {caption}"
    fig.tight_layout()
    fig.savefig(CM.FIGURES / name, dpi=DPI)
    plt.close(fig)
    CAPTIONS.append(cap)
    print("  figure:", name)


def build_historical_table(paired: pd.DataFrame) -> pd.DataFrame:
    d = paired.dropna(subset=["tw_obs", "ta_mean_corr", "elevation_m", "month"]).copy()
    d["ta"] = d["ta_mean_corr"]
    d["z"] = d["elevation_m"] / 1000.0
    d["ta_z"] = d["ta"] * d["z"]
    m = smf.ols("tw_obs ~ ta + C(month) + z + ta_z", data=d).fit()
    resid_sd = float(np.std(m.resid))
    d["tw_modelled"] = m.predict(d)
    # Observation-specific 90% prediction interval from the OLS fit (leverage
    # aware). This is conditional on the fixed-effects model and does not add
    # the between-site random effect; it replaces the previous constant
    # +/-1.645*resid_sd interval that was identical for every observation.
    sf = m.get_prediction(d).summary_frame(alpha=0.10)
    d["ci90_low"] = sf["obs_ci_lower"].values
    d["ci90_high"] = sf["obs_ci_upper"].values
    out = pd.DataFrame({
        "site_id": d["site_id"], "date": d["obs_date"], "tw_modelled": d["tw_modelled"],
        "tw_obs": d["tw_obs"], "model_name": "Model2_month_elevation",
        "ta_mean": d["ta_mean_corr"], "ta_source": "NASA_POWER_MERRA2_lapse_corrected",
        "ci90_low": d["ci90_low"], "ci90_high": d["ci90_high"], "source": d["source"],
    })
    out.to_csv(CM.TABLES / "water_temp_historical_sites.csv", index=False)
    width = float((out["ci90_high"] - out["ci90_low"]).mean())
    print(f"  water_temp_historical_sites.csv: {len(out)} rows, residual sd={resid_sd:.2f} C, "
          f"mean 90% PI width={width:.2f} C")
    return d


def fig_scatter(d: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    # single-series scatter: a per-site legend previously produced ~600 entries
    # and covered the whole plot, so line labels are annotated inline instead.
    ax.scatter(d["ta"], d["tw_obs"], s=8, alpha=0.25, color="steelblue")
    xs = np.linspace(d["ta"].min(), d["ta"].max(), 50)
    pooled = smf.ols("tw_obs ~ ta", data=d).fit()
    yfit = pooled.params["Intercept"] + pooled.params["ta"] * xs
    ax.plot(xs, yfit, "k-", lw=2)
    ax.plot(xs, xs, "r--", lw=1.2)
    i_fit = int(0.72 * len(xs))
    ax.annotate(f"pooled OLS (slope={pooled.params['ta']:.2f})",
                (xs[i_fit], yfit[i_fit]), textcoords="offset points", xytext=(4, 6),
                fontsize=8, color="k")
    i_id = int(0.10 * len(xs))
    ax.annotate("Tw = Ta", (xs[i_id], xs[i_id]), textcoords="offset points",
                xytext=(4, 6), fontsize=8, color="r")
    ax.set_xlabel("Air temperature (deg C, NASA POWER, lapse-corrected)")
    ax.set_ylabel("Observed water temperature (deg C)")
    ax.set_title(f"Tw vs Ta, Spanish rivers (N={len(d)} obs, {d['site_id'].nunique()} sites)")
    save(fig, "fig01_tw_vs_ta_scatter.png",
         f"Tw vs Ta for all sites (single series); pooled OLS and 1:1 lines labelled inline. "
         f"N={len(d)} observations, {d['site_id'].nunique()} sites. Aerated/canal sites not filtered.")


def fig_residuals(d: pd.DataFrame) -> None:
    d = d.copy()
    d["resid"] = d["tw_obs"] - d["tw_modelled"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    mm = d.groupby("month")["resid"].agg(["mean", "std", "size"]).reset_index()
    axes[0].errorbar(mm["month"], mm["mean"], yerr=mm["std"], fmt="o-", capsize=3)
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set_xlabel("Month")
    axes[0].set_ylabel("Residual (obs - model) deg C")
    axes[0].set_title(f"Residuals vs month (N={len(d)})")
    axes[1].scatter(d["elevation_m"], d["resid"], s=8, alpha=0.4)
    if d["elevation_m"].notna().sum() > 10:
        z = np.polyfit(d["elevation_m"], d["resid"], 1)
        xs = np.linspace(d["elevation_m"].min(), d["elevation_m"].max(), 50)
        axes[1].plot(xs, np.polyval(z, xs), "r-")
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].set_xlabel("Site elevation (m)")
    axes[1].set_ylabel("Residual (obs - model) deg C")
    axes[1].set_title(f"Residuals vs elevation (N={len(d)})")
    save(fig, "fig02_residuals_month_elevation.png",
         f"Model residual diagnostics. Left: monthly mean residual with 1 SD. "
         f"Right: residual vs site elevation with OLS trend. N={len(d)}.")


def fig_loso(preds: pd.DataFrame) -> None:
    # restrict to Stage-B (all-Spain) folds: the Guadalquivir folds are a subset
    # of Spain, so pooling both scopes would double-count them.
    p = preds[(preds["model_name"] == "Model2_month") & (preds["scope"] == "Spain")]
    p = p.dropna(subset=["tw_obs", "pred"])
    if p.empty:
        return

    def nse(o, pr):
        o, pr = np.asarray(o), np.asarray(pr)
        return 1 - np.sum((pr - o) ** 2) / np.sum((o - o.mean()) ** 2)

    n = nse(p["tw_obs"], p["pred"])
    rmse = float(np.sqrt(np.mean((p["pred"] - p["tw_obs"]) ** 2)))
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(p["tw_obs"], p["pred"], s=10, alpha=0.4)
    lims = [min(p["tw_obs"].min(), p["pred"].min()), max(p["tw_obs"].max(), p["pred"].max())]
    ax.plot(lims, lims, "r--", lw=1)
    ax.set_xlabel("Observed Tw (deg C)")
    ax.set_ylabel("LOSO predicted Tw (deg C)")
    ax.set_title(f"Leave-one-station-out, Stage B (N={len(p)})\nNSE={n:.2f}  RMSE={rmse:.2f} C")
    save(fig, "fig03_loso_pred_vs_obs.png",
         f"Predicted vs observed water temperature for leave-one-station-out validation "
         f"of Model 2 (month + Ta + elevation interaction). N={len(p)} held-out observations.")


def fig_map(sites: pd.DataFrame) -> None:
    g = gpd.read_file(BASIN_SHP)
    es050 = g[g["EUMASCod"].astype(str).str.startswith("ES050")]
    poly = unary_union(es050.geometry.values)
    bp = gpd.GeoSeries([poly], crs=25830).to_crs(4326)
    fig, ax = plt.subplots(figsize=(7, 8))
    bp.plot(ax=ax, facecolor="none", edgecolor="grey", linewidth=0.7)
    s = sites.dropna(subset=["lat", "lon", "elevation_m"]).copy()
    gq = s[s["in_guadalquivir_basin"]]
    oth = s[~s["in_guadalquivir_basin"]]
    sc = ax.scatter(oth["lon"], oth["lat"], c=oth["elevation_m"], cmap="terrain",
                    s=18, vmin=0, vmax=2000, alpha=0.7)
    ax.scatter(gq["lon"], gq["lat"], c=gq["elevation_m"], cmap="terrain", s=45,
               vmin=0, vmax=2000, edgecolor="k", linewidth=0.5)
    cb = plt.colorbar(sc, ax=ax, shrink=0.6)
    cb.set_label("Elevation (m, Copernicus DEM GLO-30)")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_title(f"Tw sites (N={len(s)}); ES050 basin outline highlighted")
    save(fig, f"fig04_site_map.png",
         f"Elevation-coloured map of all {len(s)} retained water-temperature sites; "
         f"the {len(gq)} Guadalquivir (ES050) sites are ringed. Basin outline = dissolved "
         f"CHG ES050 water-body catchments.")


def fig_timeseries(d: pd.DataFrame) -> None:
    gq = d[d["site_id"].isin(["ESP00016", "ESP00017", "ESP00018"])]
    if gq.empty:
        return
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for ax, sid in zip(axes, ["ESP00016", "ESP00017", "ESP00018"]):
        s = gq[gq["site_id"] == sid].sort_values("obs_date")
        if s.empty:
            continue
        ax.plot(s["obs_date"], s["tw_obs"], "o-", ms=3, label=f"Tw obs (n={len(s)})")
        ax.plot(s["obs_date"], s["ta_mean_corr"], "s--", ms=3, alpha=0.6,
                label="Ta (lapse-corrected)")
        ax.set_ylabel("deg C")
        ax.set_title(sid)
        ax.legend(fontsize=7)
    axes[-1].set_xlabel("Date")
    save(fig, "fig05_gq_timeseries.png",
         "Air temperature (lapse-corrected NASA POWER) and observed water temperature for "
         "the three Guadalquivir GEMSTAT sites (N=67, 58, 68).")


def fig_future_monthly(daily: pd.DataFrame) -> None:
    f = daily.copy()
    f = f[(f["date"].dt.year >= 2036) & (f["date"].dt.year <= 2055)]
    if f.empty:
        return
    # monthly ensemble spread per scenario; the median series is across GCMs at
    # each site-date, and the shading is the interquartile range across
    # site-days (sites x years), not a GCM-only spread.
    f["month"] = f["date"].dt.month
    fig, ax = plt.subplots(figsize=(9, 5))
    for scen, g in f.groupby("scenario"):
        mm = g.groupby("month")["tw_ensemble_median"].median()
        q = g.groupby("month")["tw_ensemble_median"].quantile([0.1, 0.25, 0.75, 0.9]).unstack()
        ax.plot(mm.index, mm.values, "-o", label=scen)
        ax.fill_between(q.index, q[0.25], q[0.75], alpha=0.2)
    ax.set_xlabel("Month")
    ax.set_ylabel("Modelled Tw (deg C, ensemble median across GCMs)")
    ax.set_title(f"2036-2055 monthly Tw by scenario (shaded IQR across site-days); "
                 f"N sites={f['site_id'].nunique()}")
    ax.legend()
    save(fig, "fig06_future_monthly_ensemble.png",
         f"Monthly ensemble-median Tw for 2036-2055 by SSP; shading is the interquartile "
         f"range across site-days (sites x years). N={f['site_id'].nunique()} sites x 20 years.")


def fig_future_delta(future: pd.DataFrame) -> None:
    if future.empty:
        return
    f = future[(future["period_start"] == 2036) & (future["period_end"] == 2055)]
    fig, ax = plt.subplots(figsize=(8, 5))
    data = [f[f["scenario"] == s]["delta_tw_mean"].dropna().values for s in
            sorted(f["scenario"].unique())]
    ax.boxplot(data, tick_labels=sorted(f["scenario"].unique()))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("delta Tw (deg C) vs each model's 1986-2005 baseline")
    ax.set_title(f"2045-window (2036-2055) water-temperature change (N={len(f)} model-site values)")
    save(fig, "fig07_future_delta_box.png",
         f"Distribution of 2036-2055 Tw change (per GCM x site) by SSP; N={len(f)} values.")


def fig_elevation_response(d: pd.DataFrame) -> None:
    rows = []
    for sid, g in d.groupby("site_id"):
        if len(g) < 15:
            continue
        try:
            m = smf.ols("tw_obs ~ ta", data=g).fit()
        except Exception:  # noqa: BLE001
            continue
        ci = m.conf_int().loc["ta"]
        rows.append({"site_id": sid, "elevation_m": g["elevation_m"].iloc[0],
                     "slope": m.params["ta"], "slope_lo": ci[0], "slope_hi": ci[1],
                     "offset": (g["tw_obs"] - g["ta"]).mean(), "n": len(g)})
    r = pd.DataFrame(rows)
    if r.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].errorbar(r["elevation_m"], r["slope"], yerr=[r["slope"] - r["slope_lo"],
                     r["slope_hi"] - r["slope"]], fmt="o", ms=4, capsize=2, alpha=0.7)
    if len(r) > 2:
        z = np.polyfit(r["elevation_m"], r["slope"], 1)
        xs = np.linspace(r["elevation_m"].min(), r["elevation_m"].max(), 50)
        axes[0].plot(xs, np.polyval(z, xs), "r-")
    axes[0].set_xlabel("Site elevation (m)")
    axes[0].set_ylabel("dTw/dTa (per-station slope)")
    axes[0].set_title(f"Slope vs elevation (N={len(r)} stations)")
    axes[1].scatter(r["elevation_m"], r["offset"], s=20)
    if len(r) > 2:
        z = np.polyfit(r["elevation_m"], r["offset"], 1)
        xs = np.linspace(r["elevation_m"].min(), r["elevation_m"].max(), 50)
        axes[1].plot(xs, np.polyval(z, xs), "r-")
    axes[1].set_xlabel("Site elevation (m)")
    axes[1].set_ylabel("Mean (Tw - Ta) deg C")
    axes[1].set_title(f"Offset vs elevation (N={len(r)} stations)")
    save(fig, "fig08_elevation_response.png",
         f"Per-station dTw/dTa slope and mean (Tw-Ta) offset against elevation, with 95% CI "
         f"on slopes. N={len(r)} stations with >=15 observations.")


def main() -> None:
    paired = pd.read_csv(CM.PROCESSED / "paired_observations.csv", parse_dates=["obs_date"],
                         dtype={"site_id": str}, low_memory=False)
    sites = pd.read_csv(CM.PROCESSED / "sites.csv", dtype={"site_id": str})
    preds = pd.read_csv(CM.MODELS / "loso_predictions.csv", parse_dates=["obs_date"],
                        dtype={"site_id": str}) if (CM.MODELS / "loso_predictions.csv").exists() else pd.DataFrame()
    print("building historical table ...")
    d = build_historical_table(paired)

    print("figures ...")
    fig_scatter(d)
    fig_residuals(d)
    if not preds.empty:
        fig_loso(preds)
    fig_map(sites)
    fig_timeseries(d)
    fig_elevation_response(d)

    fut = pd.read_csv(CM.TABLES / "water_temp_future_2045.csv")
    fig_future_delta(fut)
    dpath = CM.INTERIM / "tw_future_ensemble_daily.parquet"
    if dpath.exists():
        daily = pd.read_parquet(dpath)
        fig_future_monthly(daily)

    (CM.FIGURES / "captions.txt").write_text("\n".join(CAPTIONS), encoding="utf-8")
    CM.record("Generated figures and historical table", status="SUCCESS", checksum=False,
             path=str(CM.FIGURES), notes=f"{len(CAPTIONS)} figures with captions.txt")
    print(f"wrote {len(CAPTIONS)} figures + captions.txt")


if __name__ == "__main__":
    main()
