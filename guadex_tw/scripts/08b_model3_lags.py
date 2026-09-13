"""Model 3: lagged / thermal-inertia model for cadence-eligible sites only.

Only sites with median_obs_per_year >= 24 (Section 7 rule) are eligible; three
sites qualify. Candidate lag sets are chosen by a blocked time-split CV, never by
in-sample R2. Results are appended to outputs/tables/cv_metrics.csv.
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

CANDIDATES = {
    "M3_same_day": "tw_obs ~ ta",
    "M3_same+d3d7": "tw_obs ~ ta + ta_mean_d3 + ta_mean_d7",
    "M3_lag_bank": "tw_obs ~ ta + ta_mean_d3 + ta_mean_d7 + ta_mean_d14 + ta_mean_d30",
    "M3_d30": "tw_obs ~ ta_mean_d30",
}


def metrics(obs, pred):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    m = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[m], pred[m]
    if len(obs) < 2:
        return dict(n=len(obs), nse=np.nan, rmse_c=np.nan, bias_c=np.nan)
    resid = pred - obs
    ss_tot = np.sum((obs - obs.mean()) ** 2)
    return dict(n=len(obs), nse=1 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan,
                rmse_c=float(np.sqrt(np.mean(resid ** 2))), bias_c=float(resid.mean()))


def main() -> None:
    elig = json.loads((CM.LOGS / "lag_eligibility.json").read_text())["eligible_site_ids"]
    d = pd.read_csv(CM.PROCESSED / "paired_observations.csv", parse_dates=["obs_date"],
                    dtype={"site_id": str})
    d = d[d["site_id"].isin(elig)].copy()
    d["ta"] = d["ta_mean_corr"]
    rows = []
    single_year = []
    for sid, sg in d.groupby("site_id"):
        yrs = sorted(sg["year"].unique())
        if len(yrs) < 3:
            # no temporal block available -> in-sample only, loudly caveated
            single_year.append(sid)
            for name, f in CANDIDATES.items():
                try:
                    m = smf.ols(f, data=sg.dropna(subset=["tw_obs"])).fit()
                    p = m.predict(sg)
                except Exception:  # noqa: BLE001
                    continue
                mt = metrics(sg["tw_obs"], p)
                rows.append(dict(model_name=name, scope="Spain", test_station=sid,
                                 n_train=len(sg), n_test=len(sg),
                                 period_train=f"{yrs[0]}-{yrs[-1]}",
                                 period_test=f"{yrs[0]}-{yrs[-1]}",
                                 mae_c=np.nan, r2=np.nan,
                                 notes="IN-SAMPLE ONLY - single year, no time-blocked CV "
                                       f"possible | {f}", **mt))
            continue
        cut = yrs[-2]
        tr, te = sg[sg["year"] < cut], sg[sg["year"] >= cut]
        best = None
        for name, f in CANDIDATES.items():
            try:
                m = smf.ols(f, data=tr.dropna(subset=["tw_obs"])).fit()
                p = m.predict(te)
            except Exception:  # noqa: BLE001
                continue
            mt = metrics(te["tw_obs"], p)
            rows.append(dict(model_name=name, scope="Spain", test_station=sid,
                             n_train=len(tr), n_test=len(te),
                             period_train=f"{tr['year'].min()}-{tr['year'].max()}",
                             period_test=f"{te['year'].min()}-{te['year'].max()}",
                             mae_c=np.nan, r2=np.nan,
                             notes=f"blocked time-split | {f}", **mt))
            if best is None or (np.isfinite(mt["nse"]) and mt["nse"] > best[0]):
                best = (mt["nse"], name)
        if best:
            print(f"  {sid}: best lag model = {best[1]} (time-split NSE {best[0]:.3f})")

    out = pd.DataFrame(rows)
    if not out.empty:
        cv = pd.read_csv(CM.TABLES / "cv_metrics.csv")
        cv = cv[~cv["model_name"].isin(CANDIDATES.keys())]  # idempotent re-runs
        cv = pd.concat([cv, out], ignore_index=True)
        cv.to_csv(CM.TABLES / "cv_metrics.csv", index=False)
    summary = {"eligible_sites": elig, "n_sites_fitted": int(d["site_id"].nunique()),
               "candidate_models": list(CANDIDATES), "rows": len(out),
               "single_year_sites_no_timesplit": single_year,
               "note": "All cadence-eligible sites are confined to a single year (2024); "
                       "time-blocked CV is impossible, so lag models are in-sample only and "
                       "not used as a primary model."}
    (CM.LOGS / "model3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    CM.record("Model 3 lagged fits (eligible sites only)", status="SUCCESS", checksum=False,
              notes=json.dumps(summary))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
