import json

import pandas as pd

cv = pd.read_csv("guadex_tw/outputs/tables/cv_metrics.csv")
pd.set_option("display.width", 240)
pd.set_option("display.max_columns", 30)
p = cv[cv["test_station"].astype(str).str.contains("POOLED", na=False)]
print(p[["model_name", "scope", "n_test", "nse", "rmse_c", "mae_c", "bias_c", "r2", "notes"]]
      .to_string(index=False))
print()
print("models per scope:")
print(cv.groupby(["scope", "model_name"]).size())
print()
co = json.load(open("guadex_tw/models/model_coefficients.json", encoding="utf-8"))
for scope, v in co.items():
    print("===", scope, "n", v.get("n"), "sites", v.get("n_sites"), "area", v.get("use_area"))
    for k in ("ta", "z", "ta_z", "Intercept"):
        if k in v["params"]:
            print("  {:8s} b={:+.4f} se={:.4f} p={:.3g}".format(
                k, v["params"][k], v["bse"][k], v["pvalues"][k]))
    print("  LRT ta_z", v.get("lr_test_ta_z"), "LRT z", v.get("lr_test_z"))
    print("  AIC month", round(v.get("aic_month", 0), 1), "harmonic",
          round(v.get("aic_harmonic", 0), 1), "resid_sd", round(v.get("resid_sd", 0), 3))
    print("  summer-winter amp", v.get("summer_minus_winter_pred_c"))
    if "mixedlm" in v:
        print("  mixedlm ICC", v["mixedlm"].get("icc"), "converged",
              v["mixedlm"].get("converged"),
              "fixed", {k: round(x, 3) for k, x in
                        (v["mixedlm"].get("fixed_effects") or {}).items()
                        if k in ("Intercept", "ta", "z", "ta_z")})
