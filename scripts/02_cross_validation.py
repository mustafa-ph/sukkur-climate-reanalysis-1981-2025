#!/usr/bin/env python3
"""
Phase 2: Cross-Validation & Microclimate Analysis (1981–2025)

Description: This script checks how well the ERA5-Land and NASA POWER datasets align across daily, monthly, and yearly timeframes. 

I've also measured the "microclimate offset" between the two models. 
Because ERA5-Land has a tighter grid (~9 km), it picks up the cooling effect from the Indus River and local canal irrigation. 
NASA POWER (~50 km) is much broader, so it blends our irrigated areas with the hotter, dry desert air from the surrounding Thar/Rohri margins.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

# Setup paths depending on whether I'm running this locally or in the VM
BASE_DIR = Path("/working_dir/c_546465af28aa3c43")
if (BASE_DIR / "artifacts/file_generation/ttl=63d/output/Sukkur_Climate_Master_Daily_1981_2025.csv").exists():
    PROCESSED_DIR = BASE_DIR / "artifacts/file_generation/ttl=63d/output"
    TABLES_OUT_DIR = BASE_DIR / "artifacts/file_generation/ttl=63d/output"
else:
    REPO_ROOT = Path(__file__).resolve().parent.parent
    PROCESSED_DIR = REPO_ROOT / "data" / "processed"
    TABLES_OUT_DIR = REPO_ROOT / "results" / "tables"

TABLES_OUT_DIR.mkdir(parents=True, exist_ok=True)
DAILY_MASTER_FILE = PROCESSED_DIR / "Sukkur_Climate_Master_Daily_1981_2025.csv"


def willmott_d(x_ref, y_cmp):
    """Compute Willmott's Index of Agreement (d) between reference (ERA5-Land) and comparison (NASA POWER)."""
    x_bar = np.mean(x_ref)
    num = np.sum((y_cmp - x_ref) ** 2)
    den = np.sum((np.abs(y_cmp - x_bar) + np.abs(x_ref - x_bar)) ** 2)
    return float(1.0 - (num / den)) if den != 0 else np.nan


def compute_validation_metrics(df, era5_col, power_col, scale_label, subset_label, var_label):
    """Compute complete cross-validation metrics between ERA5-Land (ref) and NASA POWER (cmp)."""
    valid = df[[era5_col, power_col]].dropna()
    x = valid[era5_col].values  # ERA5-Land (~9 km)
    y = valid[power_col].values # NASA POWER (~50 km)
    n = len(valid)

    diff = y - x  # POWER minus ERA5-Land
    r, p_val = stats.pearsonr(x, y)
    r2 = r ** 2
    mbe = float(np.mean(diff))
    sd_diff = float(np.std(diff, ddof=1))
    mae = float(np.mean(np.abs(diff)))
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    crmse = float(np.sqrt(np.mean(((y - np.mean(y)) - (x - np.mean(x))) ** 2)))
    d_idx = willmott_d(x, y)
    ols = stats.linregress(x, y)

    return {
        "Timescale": scale_label,
        "Subset": subset_label,
        "Variable": var_label,
        "N": int(n),
        "ERA5_Mean_C": round(float(np.mean(x)), 3),
        "ERA5_SD_C": round(float(np.std(x, ddof=1)), 3),
        "POWER_Mean_C": round(float(np.mean(y)), 3),
        "POWER_SD_C": round(float(np.std(y, ddof=1)), 3),
        "Pearson_r": round(float(r), 4),
        "R_Squared": round(float(r2), 4),
        "p_value": float(p_val),
        "MBE_POWER_minus_ERA5_C": round(mbe, 3),
        "Diff_SD_C": round(sd_diff, 3),
        "MAE_C": round(mae, 3),
        "RMSE_C": round(rmse, 3),
        "Centered_RMSE_C": round(crmse, 3),
        "Willmott_d": round(d_idx, 4),
        "OLS_Slope": round(float(ols.slope), 4),
        "OLS_Intercept_C": round(float(ols.intercept), 3),
    }


def main():
    print(f"Loading Harmonized Daily Master Dataset: {DAILY_MASTER_FILE}")
    df_daily = pd.read_csv(DAILY_MASTER_FILE, parse_dates=["Date"])

    # 1. Roll up daily data into monthly averages
    df_monthly = df_daily.groupby(["Year", "Month"], as_index=False).agg(
        ERA5_Tmean_C=("ERA5_Tmean_C", "mean"),
        ERA5_Tmax_C=("ERA5_Tmax_C", "mean"),
        ERA5_Tmin_C=("ERA5_Tmin_C", "mean"),
        ERA5_DTR_C=("ERA5_DTR_C", "mean"),
        POWER_Tmean_C=("POWER_Tmean_C", "mean"),
        POWER_Tmax_C=("POWER_Tmax_C", "mean"),
        POWER_Tmin_C=("POWER_Tmin_C", "mean"),
        POWER_DTR_C=("POWER_DTR_C", "mean"),
    )

    # 2. Construct Annual Series (45 years: 1981 to 2025)
    df_annual = df_daily.groupby("Year", as_index=False).agg(
        ERA5_Tmean_C=("ERA5_Tmean_C", "mean"),
        ERA5_Tmax_C=("ERA5_Tmax_C", "mean"),
        ERA5_Tmin_C=("ERA5_Tmin_C", "mean"),
        ERA5_DTR_C=("ERA5_DTR_C", "mean"),
        POWER_Tmean_C=("POWER_Tmean_C", "mean"),
        POWER_Tmax_C=("POWER_Tmax_C", "mean"),
        POWER_Tmin_C=("POWER_Tmin_C", "mean"),
        POWER_DTR_C=("POWER_DTR_C", "mean"),
    )

    var_pairs = [
        ("ERA5_Tmean_C", "POWER_Tmean_C", "Tmean (2m Daily Mean)"),
        ("ERA5_Tmax_C", "POWER_Tmax_C", "Tmax (2m Daily Maximum)"),
        ("ERA5_Tmin_C", "POWER_Tmin_C", "Tmin (2m Daily Minimum)"),
        ("ERA5_DTR_C", "POWER_DTR_C", "DTR (Diurnal Temp. Range)"),
    ]

    summary_rows = []
    # Calculate validation metrics across all time scales
    for e_col, p_col, v_label in var_pairs:
        summary_rows.append(compute_validation_metrics(df_daily, e_col, p_col, "Daily", "All Days (1981-2025)", v_label))
        summary_rows.append(compute_validation_metrics(df_monthly, e_col, p_col, "Monthly", "All Months (1981-2025)", v_label))
        summary_rows.append(compute_validation_metrics(df_annual, e_col, p_col, "Annual", "All Years (1981-2025)", v_label))

    # Seasonal cross-validation on Daily data (Microclimate Offset Analysis)
    seasons_order = [
        "Winter (DJF)",
        "Pre-Monsoon / Spring-Summer (MAM)",
        "Summer Monsoon (JJAS)",
        "Post-Monsoon / Autumn (ON)"
    ]
    for s_name in seasons_order:
        df_s = df_daily[df_daily["Season"] == s_name]
        for e_col, p_col, v_label in var_pairs:
            summary_rows.append(compute_validation_metrics(df_s, e_col, p_col, "Seasonal (Daily)", s_name, v_label))

    df_summary = pd.DataFrame(summary_rows)
    summary_csv_path = TABLES_OUT_DIR / "Sukkur_Phase2_CrossValidation_Metrics_1981_2025.csv"
    df_summary.to_csv(summary_csv_path, index=False)
    print(f"Saved Cross-Validation Metrics Table: {summary_csv_path}")

    # 3. Monthly Microclimate Offset Profile (Jan-Dec)
    month_names = {
        1: "January", 2: "February", 3: "March", 4: "April",
        5: "May", 6: "June", 7: "July", 8: "August",
        9: "September", 10: "October", 11: "November", 12: "December"
    }
    monthly_profile_rows = []
    for m in range(1, 13):
        sub = df_daily[df_daily["Month"] == m]
        r_tmean, _ = stats.pearsonr(sub["ERA5_Tmean_C"], sub["POWER_Tmean_C"])
        r_tmax, _ = stats.pearsonr(sub["ERA5_Tmax_C"], sub["POWER_Tmax_C"])
        r_tmin, _ = stats.pearsonr(sub["ERA5_Tmin_C"], sub["POWER_Tmin_C"])

        diff_tmean = sub["POWER_Tmean_C"] - sub["ERA5_Tmean_C"]
        diff_tmax = sub["POWER_Tmax_C"] - sub["ERA5_Tmax_C"]
        diff_tmin = sub["POWER_Tmin_C"] - sub["ERA5_Tmin_C"]
        diff_dtr = sub["POWER_DTR_C"] - sub["ERA5_DTR_C"]

        monthly_profile_rows.append({
            "Month_Num": m,
            "Month_Name": month_names[m],
            "Days_Count": int(len(sub)),
            "ERA5_Tmean_C": round(float(sub["ERA5_Tmean_C"].mean()), 3),
            "POWER_Tmean_C": round(float(sub["POWER_Tmean_C"].mean()), 3),
            "MBE_Tmean_C": round(float(diff_tmean.mean()), 3),
            "RMSE_Tmean_C": round(float(np.sqrt(np.mean(diff_tmean ** 2))), 3),
            "r_Tmean": round(float(r_tmean), 4),

            "ERA5_Tmax_C": round(float(sub["ERA5_Tmax_C"].mean()), 3),
            "POWER_Tmax_C": round(float(sub["POWER_Tmax_C"].mean()), 3),
            "MBE_Tmax_C": round(float(diff_tmax.mean()), 3),
            "RMSE_Tmax_C": round(float(np.sqrt(np.mean(diff_tmax ** 2))), 3),
            "r_Tmax": round(float(r_tmax), 4),

            "ERA5_Tmin_C": round(float(sub["ERA5_Tmin_C"].mean()), 3),
            "POWER_Tmin_C": round(float(sub["POWER_Tmin_C"].mean()), 3),
            "MBE_Tmin_C": round(float(diff_tmin.mean()), 3),
            "RMSE_Tmin_C": round(float(np.sqrt(np.mean(diff_tmin ** 2))), 3),
            "r_Tmin": round(float(r_tmin), 4),

            "ERA5_DTR_C": round(float(sub["ERA5_DTR_C"].mean()), 3),
            "POWER_DTR_C": round(float(sub["POWER_DTR_C"].mean()), 3),
            "MBE_DTR_C": round(float(diff_dtr.mean()), 3),
        })

    df_monthly_profile = pd.DataFrame(monthly_profile_rows)
    monthly_csv_path = TABLES_OUT_DIR / "Sukkur_Phase2_Monthly_Microclimate_Offset_1981_2025.csv"
    df_monthly_profile.to_csv(monthly_csv_path, index=False)
    print(f"Saved Monthly Microclimate Offset Table: {monthly_csv_path}")

    # Print formatted results for inspection
    print("\n================ MULTI-TIMESCALE VALIDATION SUMMARY ================")
    sub_all = df_summary[df_summary["Timescale"].isin(["Daily", "Monthly", "Annual"])]
    print(sub_all[["Timescale", "Variable", "N", "ERA5_Mean_C", "POWER_Mean_C", "Pearson_r", "R_Squared", "MBE_POWER_minus_ERA5_C", "MAE_C", "RMSE_C", "Willmott_d"]].to_string(index=False))

    print("\n================ SEASONAL MICROCLIMATE OFFSET SUMMARY ================")
    sub_seas = df_summary[df_summary["Timescale"] == "Seasonal (Daily)"]
    print(sub_seas[["Subset", "Variable", "N", "ERA5_Mean_C", "POWER_Mean_C", "MBE_POWER_minus_ERA5_C", "MAE_C", "RMSE_C", "Pearson_r", "R_Squared"]].to_string(index=False))

    print("\n================ MONTHLY MICROCLIMATE OFFSET PROFILE ================")
    print(df_monthly_profile[["Month_Name", "ERA5_Tmax_C", "POWER_Tmax_C", "MBE_Tmax_C", "ERA5_Tmin_C", "POWER_Tmin_C", "MBE_Tmin_C", "ERA5_DTR_C", "POWER_DTR_C", "MBE_DTR_C"]].to_string(index=False))


if __name__ == "__main__":
    main()
