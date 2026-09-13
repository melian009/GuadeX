"""Attempt the EEA Waterbase (WISE-SoE) water-quality ICM direct download.

Uses the no-auth EEA discodata SQL API (Microsoft SQL Server / T-SQL):
station metadata from Waterbase_S_WISE_SpatialObject_DerivedData and
temperature observations from Waterbase_T_WISE6_DisaggregatedData.
Records SUCCESS / PARTIAL / FAILED with the exact error if it fails.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.parse
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C

BASE = "https://discodata.eea.europa.eu/sql"
SPATIAL = "[WISE_SOE].[latest].[Waterbase_S_WISE_SpatialObject_DerivedData]"
OBS = "[WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData]"
NOTES = "WISE-SoE Waterbase ICM (latest revision on discodata)"
LICENSE = "EEA open data / CC-BY"

raw_dir = C.RAW / "waterbase"
raw_dir.mkdir(parents=True, exist_ok=True)


def api(sql: str, page: int = 1, hits: int = 5000, retries: int = 4):
    url = BASE + "?query=" + urllib.parse.quote(sql) + f"&p={page}&nrOfHits={hits}"
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=300)
            j = r.json()
            if "errors" in j:
                last = j["errors"]
                time.sleep(2 * (attempt + 1))
                continue
            return j.get("results", [])
        except Exception as exc:  # noqa: BLE001
            last = repr(exc)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"discodata query failed after {retries} retries: {last}")


def paginate(sql: str, hits: int = 5000, max_pages: int = 400, label: str = "") -> pd.DataFrame:
    rows, page = [], 1
    while page <= max_pages:
        res = api(sql, page=page, hits=hits)
        rows.extend(res)
        print(f"  [{label}] page {page}: {len(res)} rows (total {len(rows)})", flush=True)
        if len(res) < hits:
            break
        page += 1
    return pd.DataFrame(rows)


def main() -> None:
    t0 = time.time()
    status, err = "SUCCESS", ""
    try:
        print("Fetching ES monitoring-site metadata ...")
        spat_sql = (
            f"select monitoringSiteIdentifier as site_id, monitoringSiteName as site_name, "
            f"waterBodyIdentifier as wb_id, waterBodyName as wb_name, rbdIdentifier as rbd_id, "
            f"rbdName as rbd_name, thematicIdIdentifier as thematic_id, reservoir as reservoir, "
            f"surfaceWaterBodyTypeCode as swb_type, lon as lon, lat as lat "
            f"from {SPATIAL} where countryCode='ES' and monitoringSiteIdentifier is not null"
        )
        spat = paginate(spat_sql, hits=5000, label="spatial")
        spat = spat.drop_duplicates(subset=["site_id"])
        spat.to_csv(C.INTERIM / "waterbase_spatial_spain.csv", index=False)
        print(f"  {len(spat)} Spanish sites")

        print("Fetching ES river temperature observations ...")
        obs_sql = (
            f"select UID as uid, monitoringSiteIdentifier as site_id, "
            f"phenomenonTimeSamplingDate as obs_date, resultObservedValue as obs_value, "
            f"resultUom as unit, resultQualityObservedValueBelowLOQ as below_loq, "
            f"parameterWaterBodyCategory as category, resultObservationStatus as obs_status "
            f"from {OBS} where countryCode='ES' and parameterWaterBodyCategory='RW' "
            f"and observedPropertyDeterminandCode='EEA_3121-01-5' "
            f"and phenomenonTimeReferenceYear>=1990"
        )
        obs = paginate(obs_sql, hits=10000, label="temp")
        obs = obs.drop_duplicates(subset=["uid"]) if not obs.empty else obs
        obs.to_csv(C.INTERIM / "waterbase_obs_raw_spain.csv", index=False)
        print(f"  {len(obs)} Spanish river temperature rows")

        if obs.empty:
            status = "PARTIAL"
            err = "query returned zero temperature rows"
            out = pd.DataFrame()
        else:
            df = obs.merge(spat, on="site_id", how="left")
            df["obs_date"] = pd.to_datetime(df["obs_date"], errors="coerce")
            df["obs_value"] = pd.to_numeric(df["obs_value"], errors="coerce")
            df = df.dropna(subset=["obs_date", "obs_value", "lat", "lon"])
            df["unit"] = "Deg C"
            df["param_code"] = "TEMP"
            df["param_name"] = "Water Temperature"
            df["source_param_code"] = "EEA_3121-01-5"
            df["source"] = "WATERBASE_DIRECT"
            df["site_country"] = "Spain"
            df["drainage_region_name"] = df["rbd_name"]
            df["detection_limit_flag"] = df["below_loq"].map(
                {True: "below_loq", False: None}
            )
            df["site_class"] = df["reservoir"].map(
                {True: "reservoir", False: "river"}
            ).fillna("unknown")
            df["obs_id"] = "WB_" + df["uid"].astype(str)
            cols = [
                "obs_id", "site_id", "site_name", "site_country", "lat", "lon",
                "obs_date", "obs_value", "unit", "param_code", "param_name",
                "source_param_code", "source", "drainage_region_name", "wb_id", "wb_name",
                "rbd_id", "site_class", "detection_limit_flag", "category", "obs_status",
            ]
            out = df[cols].copy()
            out.to_csv(C.INTERIM / "waterbase_temp_spain.csv", index=False)
            gq = out[out["rbd_id"].astype(str).str.startswith("ES050")]
            gq.to_csv(C.INTERIM / "waterbase_temp_guadalquivir.csv", index=False)
            print(f"  Spain river Tw rows: {len(out)}, sites: {out['site_id'].nunique()}")
            print(f"  Guadalquivir (ES050) rows: {len(gq)}, "
                  f"sites: {gq['site_id'].nunique() if len(gq) else 0}")
            C.record(
                "Waterbase direct (WISE-SoE discodata)",
                url="https://discodata.eea.europa.eu/sql?query=...WISE6_DisaggregatedData",
                path=str(C.INTERIM / "waterbase_temp_spain.csv"), license=LICENSE,
                version="WISE_SOE.latest (WISE6_DisaggregatedData)", status=status,
                notes=f"{len(out)} ES river Tw obs, {out['site_id'].nunique()} sites; "
                      f"{len(gq)} in Guadalquivir ES050",
            )
    except Exception as exc:  # noqa: BLE001
        status, err = "FAILED", repr(exc)
        print(f"Waterbase direct FAILED: {exc!r}")
        C.record(
            "Waterbase direct (WISE-SoE discodata)",
            url="https://discodata.eea.europa.eu/sql", status="FAILED", error=err,
            license=LICENSE, version="WISE_SOE.latest",
            notes="EEA discodata SQL API attempt",
        )
    finally:
        (C.LOGS / "waterbase_attempt.json").write_text(
            json.dumps({"status": status, "error": err,
                        "seconds": round(time.time() - t0, 1)}, indent=2),
            encoding="utf-8",
        )
        print(f"Waterbase attempt status: {status} ({round(time.time()-t0,1)} s)")


if __name__ == "__main__":
    main()
