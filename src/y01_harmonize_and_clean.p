#!/usr/bin/env python3
"""
Phase 1: Harmonization, Quality Control, and Master Dataset Creation (Sukkur 1981-2025)
ECMWF ERA5-Land (27.70N, 68.80E) vs NASA POWER MERRA-2 (27.72N, 68.82E)
"""
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = REPO_ROOT / "data" / "raw"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

POWER_FILE = RAW_DIR / "POWER_Point_Daily_19810101_20251231_027d72N_068d82E_LST.csv"
ERA5_FILE = RAW_DIR / "reanalysis-era5-land-timeseries-sfc-2m-temperaturetxkzu847.csv"


def main():
    with open(POWER_FILE, "r", encoding="utf-8") as f:
        lines = f.readlines()
    header_idx = next(i + 1 for i, line in enumerate(lines) if "-END HEADER-" in line)

    df_power = pd.read_csv(POWER_FILE, skiprows=header_idx)
    for col in ["T2M", "T2M_MAX", "T2M_MIN"]:
        df_power[col] = pd.to_numeric(df_power[col], errors="coerce")
        df_power.loc[df_power[col] <= -900, col] = np.nan

    df_power["Date"] = pd.to_datetime(
        df_power["YEAR"].astype(int).astype(str) + "-" + df_power["DOY"].astype(int).astype(str),
        format="%Y-%j"
    )
    df_power = df_power.rename(columns={
        "YEAR": "Year", "DOY": "Day_of_Year",
        "T2M": "POWER_Tmean_C", "T2M_MAX": "POWER_Tmax_C", "T2M_MIN": "POWER_Tmin_C"
    })
    df_power["POWER_DTR_C"] = (df_power["POWER_Tmax_C"] - df_power["POWER_Tmin_C"]).round(3)

    df_era5_hr = pd.read_csv(ERA5_FILE)
    df_era5_hr["valid_time_utc"] = pd.to_datetime(df_era5_hr["valid_time"])
    df_era5_hr["valid_time_pkt"] = df_era5_hr["valid_time_utc"] + pd.Timedelta(hours=5)
    df_era5_hr["t2m_C"] = df_era5_hr["t2m"] - 273.15
    df_era5_hr["Date"] = df_era5_hr["valid_time_pkt"].dt.floor("D")

    df_era5_hr_filtered = df_era5_hr[
        (df_era5_hr["Date"] >= pd.Timestamp("1981-01-01")) &
        (df_era5_hr["Date"] <= pd.Timestamp("2025-12-31"))
    ].copy()

    df_era5_daily = df_era5_hr_filtered.groupby("Date").agg(
        ERA5_Hours_Count=("t2m_C", "count"),
        ERA5_Tmean_C=("t2m_C", "mean"),
        ERA5_Tmax_C=("t2m_C", "max"),
        ERA5_Tmin_C=("t2m_C", "min"),
        ERA5_Tstd_C=("t2m_C", "std")
    ).reset_index()
    df_era5_daily["ERA5_DTR_C"] = df_era5_daily["ERA5_Tmax_C"] - df_era5_daily["ERA5_Tmin_C"]

    master_daily = pd.merge(df_power, df_era5_daily, on="Date", how="inner").sort_values("Date").reset_index(drop=True)
    master_daily["Month"] = master_daily["Date"].dt.month
    master_daily["Day"] = master_daily["Date"].dt.day

    master_daily["Consensus_Tmean_C"] = ((master_daily["POWER_Tmean_C"] + master_daily["ERA5_Tmean_C"]) / 2.0).round(3)
    master_daily["Consensus_Tmax_C"] = ((master_daily["POWER_Tmax_C"] + master_daily["ERA5_Tmax_C"]) / 2.0).round(3)
    master_daily["Consensus_Tmin_C"] = ((master_daily["POWER_Tmin_C"] + master_daily["ERA5_Tmin_C"]) / 2.0).round(3)
    master_daily["Consensus_DTR_C"] = ((master_daily["POWER_DTR_C"] + master_daily["ERA5_DTR_C"]) / 2.0).round(3)

    master_daily["Diff_Tmean_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmean_C"] - master_daily["ERA5_Tmean_C"]).round(3)
    master_daily["Diff_Tmax_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmax_C"] - master_daily["ERA5_Tmax_C"]).round(3)
    master_daily["Diff_Tmin_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmin_C"] - master_daily["ERA5_Tmin_C"]).round(3)

    master_daily["ERA5_HotDay_40C"] = (master_daily["ERA5_Tmax_C"] >= 40.0).astype(int)
    master_daily["ERA5_ExtremeHeat_45C"] = (master_daily["ERA5_Tmax_C"] >= 45.0).astype(int)
    master_daily["ERA5_WarmNight_30C"] = (master_daily["ERA5_Tmin_C"] >= 30.0).astype(int)
    master_daily["POWER_HotDay_40C"] = (master_daily["POWER_Tmax_C"] >= 40.0).astype(int)
    master_daily["POWER_ExtremeHeat_45C"] = (master_daily["POWER_Tmax_C"] >= 45.0).astype(int)
    master_daily["POWER_WarmNight_30C"] = (master_daily["POWER_Tmin_C"] >= 30.0).astype(int)

    for col in ["ERA5_Tmean_C", "ERA5_Tmax_C", "ERA5_Tmin_C", "ERA5_Tstd_C", "ERA5_DTR_C",
                "POWER_Tmean_C", "POWER_Tmax_C", "POWER_Tmin_C", "POWER_DTR_C"]:
        master_daily[col] = master_daily[col].round(3)

    master_daily.to_csv(PROCESSED_DIR / "Sukkur_Climate_Master_Daily_1981_2025.csv", index=False)


if __name__ == "__main__":
    main()
