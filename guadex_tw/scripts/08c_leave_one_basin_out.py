"""Leave-one-basin-out transfer test (Section 9.3).

Fit on all Spanish sites outside the Guadalquivir basin, predict every
Guadalquivir site. This is the decisive test of whether the national
elevation-aware model transfers to the study basin. Results appended to
outputs/tables/cv_metrics.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as CM

MODELS = {"LOBO_Model2_elevation": "tw_obs ~ ta + C(month) + z + ta_z",
          "LOBO_Model2_no_elevation": "tw_obs ~ ta + C(month)"}


def metrics(obs, pred):
    obs, pred = np.asarray(obs, float), np.asarray(pred, float)
    m = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[m], pred[m]
    if len(obs) < 2:
        return dict(n=len(obs), nse=np.nan, rmse_c=np.nan, mae_c=np.nan, bias_c=np.nan, r2=np.nan)
    r = pred - obs
    ss_tot = float(np.sum((obs - obs.mean()) ** 2))
    return dict(n=len(obs), nse=1 - float(np.sum(r ** 2)) / ss_tot if ss_tot > 0 else np.nan,
                rmse_c=float(np.sqrt(np.mean(r ** 2))), mae_c=float(np.mean(np.abs(r))),
                bias_c=float(r.mean()),
                r2=float(np.corrcoef(obs, pred)[0, 1] ** 2) if np.std(pred) > 0 else np.nan)


def main() -> None:
    d = pd.read_csv(CM.PROCESSED / "paired_observations.csv", parse_dates=["obs_date"],
                    dtype={"site_id": str}, low_memory=False)
    d = d.dropna(subset=["tw_obs", "ta_mean_corr", "elevation_m", "month"]).copy()
    d["ta"] = d["ta_mean_corr"]
    d["z"] = d["elevation_m"] / 1000.0
    d["ta_z"] = d["ta"] * d["z"]
    d["month"] = d["month"].astype(int)
    train = d[~d["in_guadalquivir_basin"].astype(bool)]
    test = d[d["in_guadalquivir_basin"].astype(bool)]
    print(f"train {len(train)} obs / {train['site_id'].nunique()} sites (non-Guadalquivir)")
    print(f"test  {len(test)} obs / {test['site_id'].nunique()} sites (Guadalquivir)")

    rows = []
    preds = []
    for name, f in MODELS.items():
        m = smf.ols(f, data=train).fit()
        p = m.predict(test)
        mt = metrics(test["tw_obs"], p)
        rows.append(dict(model_name=name, scope="Spain_to_Guadalquivir", test_station="ALL_GQ",
                         n_train=len(train), n_test=len(test),
                         period_train=f"{train['year'].min()}-{train['year'].max()}",
                         period_test=f"{test['year'].min()}-{test['year'].max()}",
                         notes="leave-one-basin-out", **mt))
        for sid, sg in test.assign(pred=np.asarray(p, float)).groupby("site_id"):
            mm = metrics(sg["tw_obs"], sg["pred"])
            rows.append(dict(model_name=name, scope="Spain_to_Guadalquivir", test_station=sid,
                             n_train=len(train), n_test=len(sg),
                             period_train=f"{train['year'].min()}-{train['year'].max()}",
                             period_test=f"{sg['year'].min()}-{sg['year'].max()}",
                             notes="leave-one-basin-out", **mm))
            preds.append(pd.DataFrame({"site_id": sid, "obs_date": sg["obs_date"],
                                       "tw_obs": sg["tw_obs"], "pred": sg["pred"],
                                       "model_name": name, "scope": "Spain_to_Guadalquivir"}))
        print(f"  {name}: NSE={mt['nse']:.3f} RMSE={mt['rmse_c']:.2f} C n={mt['n']}")

    # baselines
    for bname in ("baseline_Tw=Ta", "baseline_Ta+offset", "baseline_monthly_climatology"):
        if bname == "baseline_Tw=Ta":
            p = test["ta"]
        elif bname == "baseline_Ta+offset":
            c = float((train["tw_obs"] - train["ta"]).mean())
            p = test["ta"] + c
        else:
            mm = train.groupby("month")["tw_obs"].mean()
            p = test["month"].map(mm)
        mt = metrics(test["tw_obs"], p)
        rows.append(dict(model_name=bname, scope="Spain_to_Guadalquivir", test_station="ALL_GQ",
                         n_train=len(train), n_test=len(test),
                         period_train=f"{train['year'].min()}-{train['year'].max()}",
                         period_test=f"{test['year'].min()}-{test['year'].max()}",
                         notes="leave-one-basin-out baseline", **mt))
        print(f"  {bname}: NSE={mt['nse']:.3f} RMSE={mt['rmse_c']:.2f} C")

    out = pd.DataFrame(rows)
    cv = pd.read_csv(CM.TABLES / "cv_metrics.csv")
    cv = cv[cv["scope"] != "Spain_to_Guadalquivir"]  # idempotent
    cv = pd.concat([cv, out], ignore_index=True)
    cv.to_csv(CM.TABLES / "cv_metrics.csv", index=False)
    if preds:
        pd.concat(preds, ignore_index=True).to_csv(CM.MODELS / "lobo_predictions.csv", index=False)
    CM.record("Leave-one-basin-out transfer test", status="SUCCESS", checksum=False,
              path=str(CM.TABLES / "cv_metrics.csv"),
              notes=f"{len(rows)} rows; train non-GQ {len(train)}, test GQ {len(test)}")
    print("wrote leave-one-basin-out rows")


if __name__ == "__main__":
    main()
