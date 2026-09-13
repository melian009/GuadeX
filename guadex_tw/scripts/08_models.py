"""Fit Tw~Ta models (Section 8) and run the validation protocol (Section 9).

Stage A = Guadalquivir sites; Stage B = all Spanish sites.
Outputs: outputs/tables/cv_metrics.csv, models/model_coefficients.json,
models/model_summary.txt, models/loso_predictions_*.csv.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as CM

warnings.filterwarnings("ignore")
GAMMA = -0.0065
CADENCE_RULE = 24


def metrics(obs, pred):
    obs = np.asarray(obs, float)
    pred = np.asarray(pred, float)
    m = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[m], pred[m]
    n = len(obs)
    if n < 2:
        return dict(n=n, nse=np.nan, rmse_c=np.nan, mae_c=np.nan, bias_c=np.nan, r2=np.nan)
    resid = pred - obs
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    return dict(
        n=n,
        nse=1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
        rmse_c=float(np.sqrt(np.mean(resid ** 2))),
        mae_c=float(np.mean(np.abs(resid))),
        bias_c=float(np.mean(resid)),
        r2=float(np.corrcoef(obs, pred)[0, 1] ** 2) if np.std(pred) > 0 else np.nan,
    )


def prep(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["ta"] = d["ta_mean_corr"]
    d["z"] = d["elevation_m"] / 1000.0
    d["ta_z"] = d["ta"] * d["z"]
    area = pd.to_numeric(d.get("upstream_basin_area"), errors="coerce")
    # keep the raw log-area; imputation is done per scope below (upstream_basin_area
    # is only ~55% populated nationally and absent for Guadalquivir sites).
    d["log_area_raw"] = np.log10(area.where(area > 0))
    d["log_area"] = d["log_area_raw"]
    d["ta_logarea"] = np.nan
    d["month"] = d["month"].astype(int)
    d = d.dropna(subset=["tw_obs", "ta", "z", "month"])
    return d


def baselines_loso(df: pd.DataFrame, scope: str) -> tuple[list, list]:
    rows, preds = [], []
    for sid, te in df.groupby("site_id"):
        tr = df[df["site_id"] != sid]
        for name in ("baseline_Tw=Ta", "baseline_Ta+offset", "baseline_monthly_climatology"):
            if name == "baseline_Tw=Ta":
                p = te["ta"]
            elif name == "baseline_Ta+offset":
                c = float((tr["tw_obs"] - tr["ta"]).mean())
                p = te["ta"] + c
            else:
                mm = tr.groupby("month")["tw_obs"].mean()
                p = te["month"].map(mm)
            mt = metrics(te["tw_obs"], p)
            rows.append(dict(model_name=name, scope=scope, test_station=sid,
                             n_train=len(tr), n_test=len(te),
                             period_train=f"{tr['year'].min()}-{tr['year'].max()}",
                             period_test=f"{te['year'].min()}-{te['year'].max()}",
                             notes="LOSO", **mt))
            preds.append(pd.DataFrame({"site_id": sid, "obs_date": te["obs_date"].values,
                                       "tw_obs": te["tw_obs"].values,
                                       "pred": np.asarray(p, float), "model_name": name,
                                       "scope": scope}))
    return rows, preds


def loso_ols(df: pd.DataFrame, formula: str, model_name: str, scope: str) -> tuple[list, list]:
    rows, preds = [], []
    for sid, te in df.groupby("site_id"):
        tr = df[df["site_id"] != sid]
        try:
            m = smf.ols(formula, data=tr).fit()
            p = m.predict(te)
        except Exception as exc:  # noqa: BLE001
            rows.append(dict(model_name=model_name, scope=scope, test_station=sid,
                             n_train=len(tr), n_test=len(te), notes=f"FIT_FAILED {exc!r}"))
            continue
        mt = metrics(te["tw_obs"], p)
        rows.append(dict(model_name=model_name, scope=scope, test_station=sid,
                         n_train=len(tr), n_test=len(te),
                         period_train=f"{tr['year'].min()}-{tr['year'].max()}",
                         period_test=f"{te['year'].min()}-{te['year'].max()}",
                         notes="LOSO", **mt))
        preds.append(pd.DataFrame({"site_id": sid, "obs_date": te["obs_date"].values,
                                   "tw_obs": te["tw_obs"].values,
                                   "pred": np.asarray(p, float), "model_name": model_name,
                                   "scope": scope}))
    return rows, preds


def pooled_mean_row(rows: list, model_name: str, scope: str, note: str = "") -> dict:
    """Aggregate fold metrics into a pooled row (mean across folds, N-weighted)."""
    sub = [r for r in rows if r["model_name"] == model_name and r["scope"] == scope
           and "nse" in r and np.isfinite(r.get("nse", np.nan))]
    if not sub:
        return dict(model_name=model_name, scope=scope, test_station="POOLED_FOLDS",
                    notes=note or "no folds")
    w = np.array([r["n_test"] for r in sub], float)
    out = dict(model_name=model_name, scope=scope, test_station="POOLED_FOLDS_MEAN",
               n_train=int(np.sum([r["n_train"] for r in sub])),
               n_test=int(np.sum(w)),
               period_train="", period_test="", notes=note or "mean over LOSO folds")
    for k in ("nse", "rmse_c", "mae_c", "bias_c", "r2"):
        v = np.array([r[k] for r in sub], float)
        out[k] = float(np.nansum(v * w) / np.nansum(np.where(np.isfinite(v), w, 0))) if k != "nse" else float(np.nanmean(v))
    return out


def per_station_month_model1(df: pd.DataFrame, scope: str) -> list:
    """Model 1: per-station per-month OLS, n>=15, else pooled-month fallback."""
    rows = []
    pooled_month = {
        int(mo): smf.ols("tw_obs ~ ta", data=g).fit()
        for mo, g in df.groupby("month") if len(g) >= 15
    }
    for sid, sg in df.groupby("site_id"):
        fitted = []
        for mo, mg in sg.groupby("month"):
            if len(mg) >= 15:
                m = smf.ols("tw_obs ~ ta", data=mg).fit()
                fitted.append(pd.DataFrame({"tw_obs": mg["tw_obs"], "pred": m.predict(mg),
                                            "n": len(mg), "fallback": False}))
            else:
                pm = pooled_month.get(mo)
                if pm is not None:
                    fitted.append(pd.DataFrame({"tw_obs": mg["tw_obs"],
                                                "pred": pm.predict(mg), "n": len(mg),
                                                "fallback": True}))
        if fitted:
            f = pd.concat(fitted)
            mt = metrics(f["tw_obs"], f["pred"])
            rows.append(dict(model_name="Model1_per_station_month", scope=scope,
                             test_station=sid, n_train=len(sg), n_test=len(f),
                             period_train=f"{sg['year'].min()}-{sg['year'].max()}",
                             period_test=f"{sg['year'].min()}-{sg['year'].max()}",
                             notes=f"in-sample; months<15 use pooled fallback "
                                   f"({int(f['fallback'].sum())} obs)", **mt))
    return rows


def blocked_time_split(df: pd.DataFrame, formula: str, model_name: str, scope: str) -> list:
    """Per station: train on all but the last two distinct years, test on those."""
    rows = []
    for sid, sg in df.groupby("site_id"):
        yrs = sorted(sg["year"].unique())
        if len(yrs) < 4:
            continue
        cut = yrs[-2]
        tr, te = sg[sg["year"] < cut], sg[sg["year"] >= cut]
        if len(te) < 2 or len(tr) < 10:
            continue
        try:
            m = smf.ols(formula, data=tr).fit()
            p = m.predict(te)
        except Exception:  # noqa: BLE001
            continue
        mt = metrics(te["tw_obs"], p)
        rows.append(dict(model_name=model_name, scope=scope, test_station=sid,
                         n_train=len(tr), n_test=len(te),
                         period_train=f"{tr['year'].min()}-{tr['year'].max()}",
                         period_test=f"{te['year'].min()}-{te['year'].max()}",
                         notes=f"blocked time-split; cutoff {cut}", **mt))
    return rows


def main() -> None:
    paired = pd.read_csv(CM.PROCESSED / "paired_observations.csv", parse_dates=["obs_date"],
                         dtype={"site_id": str})
    d = prep(paired)
    print(f"paired rows with Tw+Ta+z: {len(d):,} across {d['site_id'].nunique()} sites")

    area_cov = d["log_area"].notna().mean()
    print(f"log(upstream_basin_area) coverage: {area_cov:.1%}")
    use_area = area_cov > 0.5

    scopes = {
        "Guadalquivir": d[d["in_guadalquivir_basin"]],
        "Spain": d,
    }
    cv_rows, pred_frames = [], []
    coefs = {}
    summary_lines = []

    for scope, sdf in scopes.items():
        sdf = sdf.copy()
        print(f"\n=== Scope {scope}: {len(sdf):,} obs / {sdf['site_id'].nunique()} sites ===")
        summary_lines.append(f"=== {scope}: n={len(sdf)}, sites={sdf['site_id'].nunique()} ===")
        if len(sdf) < 20:
            summary_lines.append("  too few observations; skipped")
            continue

        # --- baselines + LOSO ---
        br, bp = baselines_loso(sdf, scope)
        cv_rows += br
        pred_frames += bp
        for nm in ("baseline_Tw=Ta", "baseline_Ta+offset", "baseline_monthly_climatology"):
            cv_rows.append(pooled_mean_row(br, nm, scope, "LOSO baseline"))

        # --- Model 1 per-station/month (in-sample) ---
        m1 = per_station_month_model1(sdf, scope)
        cv_rows += m1
        m1p = pooled_mean_row(m1, "Model1_per_station_month", scope, "in-sample mean over stations")
        cv_rows.append(m1p)

        # --- Model 2 core ---
        scope_area_cov = float(sdf["log_area_raw"].notna().mean())
        scope_use_area = bool(use_area and scope_area_cov > 0.5)
        # impute log-area within the scope so the interaction can be fit
        med_area = sdf["log_area_raw"].median()
        sdf["log_area"] = (sdf["log_area_raw"].fillna(med_area) if scope_use_area
                           else sdf["log_area_raw"].fillna(0.0))
        sdf["ta_logarea"] = sdf["ta"] * sdf["log_area"]
        print(f"  basin-area coverage in scope: {scope_area_cov:.1%} "
              f"(area term {'on' if scope_use_area else 'off'})")
        term_area = " + ta_logarea" if scope_use_area else ""
        f_month = f"tw_obs ~ ta + C(month) + z + ta_z{term_area}"
        f_noz = f"tw_obs ~ ta + C(month){term_area}"
        f_zonly = f"tw_obs ~ ta + C(month) + z{term_area}"      # reduced for b3 (ta:z)
        f_not_z = f"tw_obs ~ ta + C(month) + ta_z{term_area}"   # reduced for b2 (z)
        f_harm = f"tw_obs ~ ta + sin1 + cos1 + sin2 + cos2 + z + ta_z{term_area}"

        sdf = sdf.copy()
        sdf["sin1"] = np.sin(2 * np.pi * sdf["month"] / 12)
        sdf["cos1"] = np.cos(2 * np.pi * sdf["month"] / 12)
        sdf["sin2"] = np.sin(4 * np.pi * sdf["month"] / 12)
        sdf["cos2"] = np.cos(4 * np.pi * sdf["month"] / 12)

        full = smf.ols(f_month, data=sdf).fit()
        harm = smf.ols(f_harm, data=sdf).fit()
        zonly = smf.ols(f_zonly, data=sdf).fit()
        not_z = smf.ols(f_not_z, data=sdf).fit()
        summary_lines.append(f"\n-- Model2 month-dummies AIC={full.aic:.1f} R2={full.rsquared:.3f}")
        summary_lines.append(str(full.summary()))
        summary_lines.append(f"\n-- Model2 harmonics AIC={harm.aic:.1f} R2={harm.rsquared:.3f}")
        summary_lines.append(str(harm.summary()))

        # LRT for ta_z and z
        try:
            lr_ta_z = full.compare_lr_test(zonly)
            lr_z = full.compare_lr_test(not_z)
        except Exception as exc:  # noqa: BLE001
            lr_ta_z = lr_z = (np.nan, np.nan, str(exc))
        coefs[scope] = {
            "formula": f_month,
            "params": {k: float(v) for k, v in full.params.items()},
            "bse": {k: float(v) for k, v in full.bse.items()},
            "pvalues": {k: float(v) for k, v in full.pvalues.items()},
            "aic_month": float(full.aic), "aic_harmonic": float(harm.aic),
            "r2_month": float(full.rsquared), "r2_harmonic": float(harm.rsquared),
            "lr_test_ta_z": {"stat": float(lr_ta_z[0]), "p": float(lr_ta_z[1])},
            "lr_test_z": {"stat": float(lr_z[0]), "p": float(lr_z[1])},
            "resid_sd": float(np.std(full.resid)),
            "n": int(len(sdf)), "n_sites": int(sdf["site_id"].nunique()),
            "use_area": bool(scope_use_area), "area_coverage": float(scope_area_cov),
        }
        best_formula = f_month if full.aic <= harm.aic else f_harm
        best_name = "Model2_month" if full.aic <= harm.aic else "Model2_harmonic"
        summary_lines.append(f"  best seasonal form: {best_name}")

        for fml, nm in [(best_formula, best_name), (f_noz, "Model2_no_elevation")]:
            r, p = loso_ols(sdf, fml, nm, scope)
            cv_rows += r
            pred_frames += p
            cv_rows.append(pooled_mean_row(r, nm, scope))

        # Model 4 mixed effects (full-data fit; fixed-effects-only LOSO approximation)
        try:
            md = smf.mixedlm("tw_obs ~ ta + z + ta_z + C(month)", sdf,
                             groups=sdf["site_id"], re_formula="~ta")
            mf = md.fit(method="lbfgs", maxiter=300)
            re_var = mf.cov_re.iloc[0, 0]
            resid_var = float(mf.scale)
            icc = float(re_var / (re_var + resid_var)) if (re_var + resid_var) > 0 else np.nan
            blups = mf.random_effects
            summary_lines.append(f"\n-- Model4 MixedLM converged={mf.converged} "
                                 f"ICC={icc:.3f} resid_var={resid_var:.3f}")
            summary_lines.append(str(mf.summary()))
            # Population-level (RE=0) predictions for every station from the full fit
            r, p, r_true = [], [], []
            rng = np.random.default_rng(42)
            true_folds = set(rng.choice(sorted(sdf["site_id"].unique()),
                                        size=min(20, sdf["site_id"].nunique()), replace=False))
            for sid, te in sdf.groupby("site_id"):
                pred = _mixedlm_fixed_predict(mf, te)
                mt = metrics(te["tw_obs"], pred)
                r.append(dict(model_name="Model4_mixedlm_population", scope=scope,
                              test_station=sid, n_train=len(sdf) - len(te), n_test=len(te),
                              period_train=f"{sdf['year'].min()}-{sdf['year'].max()}",
                              period_test=f"{te['year'].min()}-{te['year'].max()}",
                              notes="population-level fixed effects (RE=0); full-data fit", **mt))
                p.append(pd.DataFrame({"site_id": sid, "obs_date": te["obs_date"].values,
                                       "tw_obs": te["tw_obs"].values,
                                       "pred": np.asarray(pred, float),
                                       "model_name": "Model4_mixedlm_population", "scope": scope}))
                if sid in true_folds:
                    tr = sdf[sdf["site_id"] != sid]
                    try:
                        mm = smf.mixedlm("tw_obs ~ ta + z + ta_z + C(month)", tr,
                                         groups=tr["site_id"], re_formula="~ta").fit(
                            method="lbfgs", maxiter=300)
                        pred2 = _mixedlm_fixed_predict(mm, te)
                    except Exception:  # noqa: BLE001
                        continue
                    mtt = metrics(te["tw_obs"], pred2)
                    r_true.append(dict(model_name="Model4_mixedlm", scope=scope,
                                       test_station=sid, n_train=len(tr), n_test=len(te),
                                       period_train=f"{tr['year'].min()}-{tr['year'].max()}",
                                       period_test=f"{te['year'].min()}-{te['year'].max()}",
                                       notes="true LOSO fixed-effects (RE=0); 20 held-out stations",
                                       **mtt))
                    p.append(pd.DataFrame({"site_id": sid, "obs_date": te["obs_date"].values,
                                           "tw_obs": te["tw_obs"].values,
                                           "pred": np.asarray(pred2, float),
                                           "model_name": "Model4_mixedlm", "scope": scope}))
            cv_rows += r
            pred_frames += p
            cv_rows.append(pooled_mean_row(r, "Model4_mixedlm_population", scope,
                                           "population-level RE=0"))
            if r_true:
                cv_rows += r_true
                cv_rows.append(pooled_mean_row(r_true, "Model4_mixedlm", scope,
                                               "true LOSO (20 held-out stations)"))
            if scope == "Spain":
                blup_df = pd.DataFrame({
                    "site_id": list(blups.keys()),
                    "blup_intercept": [float(v.iloc[0]) for v in blups.values()],
                    "blup_ta_slope": [float(v.iloc[1]) if len(v) > 1 else np.nan
                                      for v in blups.values()],
                })
                blup_df.to_csv(CM.MODELS / "mixedlm_blups.csv", index=False)
            coefs[scope]["mixedlm"] = {
                "converged": bool(mf.converged), "icc": icc,
                "resid_var": resid_var, "random_intercept_var": float(re_var),
                "fixed_effects": {k: float(v) for k, v in mf.params.items()},
                "fixed_bse": {k: float(v) for k, v in mf.bse.items()},
            }
        except Exception as exc:  # noqa: BLE001
            summary_lines.append(f"\n-- Model4 MixedLM FAILED: {exc!r}")
            coefs[scope]["mixedlm"] = {"error": repr(exc)}

        # --- blocked time split (Model1 & Model2) ---
        cv_rows += blocked_time_split(sdf, "tw_obs ~ ta + C(month)", "Model1_time_split", scope)
        cv_rows += blocked_time_split(sdf, f_month, "Model2_time_split", scope)

        # sanity check on Model 2 predictions
        sdf["pred_m2"] = full.predict(sdf)
        amp = []
        for sid, sg in sdf.groupby("site_id"):
            summer = sg[sg["month"].isin([7, 8])]["pred_m2"].mean()
            winter = sg[sg["month"].isin([1, 2])]["pred_m2"].mean()
            amp.append(summer - winter)
        coefs[scope]["summer_minus_winter_pred_c"] = {
            "min": float(np.nanmin(amp)), "max": float(np.nanmax(amp)),
            "mean": float(np.nanmean(amp)),
        }
        summary_lines.append(f"  summer-winter predicted amplitude: "
                             f"min={np.nanmin(amp):.2f} max={np.nanmax(amp):.2f}")

    cv = pd.DataFrame(cv_rows)
    cv.to_csv(CM.TABLES / "cv_metrics.csv", index=False)
    all_preds = pd.concat(pred_frames, ignore_index=True) if pred_frames else pd.DataFrame()
    all_preds.to_csv(CM.MODELS / "loso_predictions.csv", index=False)
    (CM.MODELS / "model_coefficients.json").write_text(json.dumps(coefs, indent=2, default=str),
                                                      encoding="utf-8")
    (CM.MODELS / "model_summary.txt").write_text("\n".join(summary_lines), encoding="utf-8")
    print("wrote cv_metrics.csv, model_coefficients.json, model_summary.txt, loso_predictions.csv")

    # cadence rule for lag models
    sites = pd.read_csv(CM.PROCESSED / "sites.csv", dtype={"site_id": str})
    elig = sites[sites["median_obs_per_year"] >= CADENCE_RULE]
    lag_note = {
        "cadence_rule_obs_per_year": CADENCE_RULE,
        "n_sites_eligible": int(len(elig)),
        "eligible_site_ids": list(elig["site_id"]),
    }
    (CM.LOGS / "lag_eligibility.json").write_text(json.dumps(lag_note, indent=2), encoding="utf-8")
    print("lag-eligible sites:", lag_note["n_sites_eligible"])
    CM.record("Generated model outputs", status="SUCCESS", checksum=False,
             notes="cv_metrics.csv, model_coefficients.json, loso_predictions.csv")


def _mixedlm_fixed_predict(model, new: pd.DataFrame) -> np.ndarray:
    """Predict from a MixedLM using fixed effects only (random effect = 0).

    Built manually from the known fixed-effect parameter names so it does not
    depend on patsy design-info internals.
    """
    params = model.params
    p = np.full(len(new), float(params.get("Intercept", 0.0)))
    p = p + float(params.get("ta", 0.0)) * new["ta"].to_numpy()
    p = p + float(params.get("z", 0.0)) * new["z"].to_numpy()
    p = p + float(params.get("ta_z", 0.0)) * new["ta_z"].to_numpy()
    months = new["month"].to_numpy()
    for k, v in params.items():
        if k.startswith("C(month)[T."):
            m = int(k.split("T.")[1].rstrip("]"))
            p = p + float(v) * (months == m)
    return p


if __name__ == "__main__":
    main()
