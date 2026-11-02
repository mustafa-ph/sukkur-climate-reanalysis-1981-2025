#!/usr/bin/env python3
"""
Harmonization, Quality Control, and Climatological Trend Analysis of 2-Meter Air Temperature
over Sukkur, Sindh, Pakistan (1981–2025): A Comparative Reanalysis of ECMWF ERA5-Land and NASA POWER (MERRA-2)

Author: Climate Research Pipeline
Location: Sukkur, Sindh, Pakistan (ERA5-Land grid: 27.70°N, 68.80°E; NASA POWER point: 27.72°N, 68.82°E)
Time Zone Standard: Pakistan Standard Time (PKT, UTC+5) / Local Solar Time (LST)
Period: January 1, 1981 – December 31, 2025 (45 complete years; 16,436 days)
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats

# Paths
BASE_DIR = "/working_dir/c_546465af28aa3c43"
UPLOAD_DIR = os.path.join(BASE_DIR, "user_uploaded_files/ttl=63d/user_uploaded")
OUTPUT_DIR = os.path.join(BASE_DIR, "artifacts/file_generation/ttl=63d/output")
INTER_DIR = os.path.join(BASE_DIR, "artifacts/file_generation/ttl=63d/intermediate")

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(INTER_DIR, exist_ok=True)

POWER_FILE = os.path.join(UPLOAD_DIR, "POWER_Point_Daily_19810101_20251231_027d72N_068d82E_LST.csv")
ERA5_FILE = os.path.join(UPLOAD_DIR, "reanalysis-era5-land-timeseries-sfc-2m-temperaturetxkzu847.csv")


def sens_slope_and_mk(x, y):
    """Compute Mann-Kendall trend test and Theil-Sen slope estimator."""
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_c = np.asarray(x[mask], dtype=float)
    y_c = np.asarray(y[mask], dtype=float)
    n = len(y_c)
    if n < 3:
        return {"slope_yr": np.nan, "slope_dec": np.nan, "p_value": np.nan, "z_score": np.nan, "tau": np.nan}

    # Theil-Sen slope
    res_sen = stats.theilslopes(y_c, x_c, 0.95)
    slope_yr = float(res_sen[0])
    intercept = float(res_sen[1])
    ci_low_yr = float(res_sen[2])
    ci_high_yr = float(res_sen[3])

    # Mann-Kendall S statistic
    s = 0
    for k in range(n - 1):
        s += np.sum(np.sign(y_c[k + 1:] - y_c[k]))

    # Variance with tie correction
    unique, counts = np.unique(y_c, return_counts=True)
    ties = counts[counts > 1]
    tie_term = np.sum(ties * (ties - 1) * (2 * ties + 5)) if len(ties) > 0 else 0
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    p_val = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    tau, _ = stats.kendalltau(x_c, y_c)
    ols = stats.linregress(x_c, y_c)

    return {
        "slope_yr": round(slope_yr, 5),
        "slope_dec": round(slope_yr * 10.0, 4),
        "ci_low_dec": round(ci_low_yr * 10.0, 4),
        "ci_high_dec": round(ci_high_yr * 10.0, 4),
        "ols_slope_dec": round(float(ols.slope) * 10.0, 4),
        "ols_p_value": round(float(ols.pvalue), 5),
        "mk_z": round(float(z), 4),
        "mk_p": round(float(p_val), 5),
        "kendall_tau": round(float(tau), 4),
        "significant_05": bool(p_val < 0.05),
    }


def max_spell_length(bool_series):
    """Calculate maximum consecutive True days and number of spells >= 3 and >= 5 days."""
    arr = np.asarray(bool_series, dtype=int)
    if len(arr) == 0 or np.sum(arr) == 0:
        return 0, 0, 0
    # Find run lengths
    padded = np.concatenate(([0], arr, [0]))
    diffs = np.diff(padded)
    starts = np.where(diffs == 1)[0]
    ends = np.where(diffs == -1)[0]
    lengths = ends - starts
    max_len = int(np.max(lengths)) if len(lengths) > 0 else 0
    spells_3d = int(np.sum(lengths >= 3))
    spells_5d = int(np.sum(lengths >= 5))
    return max_len, spells_3d, spells_5d


def main():
    print("1. Loading and cleaning NASA POWER daily dataset...")
    # Skip header lines up to -END HEADER-
    with open(POWER_FILE, "r") as f:
        lines = f.readlines()
    header_idx = 0
    for i, line in enumerate(lines):
        if "-END HEADER-" in line:
            header_idx = i + 1
            break

    df_power = pd.read_csv(POWER_FILE, skiprows=header_idx)
    # Replace -999 missing value flags with NaN
    for col in ["T2M", "T2M_MAX", "T2M_MIN"]:
        df_power[col] = pd.to_numeric(df_power[col], errors="coerce")
        df_power.loc[df_power[col] <= -900, col] = np.nan

    # Construct calendar Date from YEAR and DOY
    df_power["Date"] = pd.to_datetime(
        df_power["YEAR"].astype(int).astype(str) + "-" + df_power["DOY"].astype(int).astype(str),
        format="%Y-%j"
    )
    df_power = df_power.rename(columns={
        "YEAR": "Year",
        "DOY": "Day_of_Year",
        "T2M": "POWER_Tmean_C",
        "T2M_MAX": "POWER_Tmax_C",
        "T2M_MIN": "POWER_Tmin_C"
    })
    df_power["POWER_DTR_C"] = (df_power["POWER_Tmax_C"] - df_power["POWER_Tmin_C"]).round(3)

    print("2. Loading and cleaning ERA5-Land hourly dataset...")
    df_era5_hr = pd.read_csv(ERA5_FILE)
    df_era5_hr["valid_time_utc"] = pd.to_datetime(df_era5_hr["valid_time"])
    # Convert UTC to Pakistan Standard Time (PKT = UTC + 5 hours) to align with NASA POWER LST (~UTC+4.59h)
    df_era5_hr["valid_time_pkt"] = df_era5_hr["valid_time_utc"] + pd.Timedelta(hours=5)
    # Convert Kelvin to Celsius
    df_era5_hr["t2m_C"] = df_era5_hr["t2m"] - 273.15
    df_era5_hr["Date"] = df_era5_hr["valid_time_pkt"].dt.floor("D")

    # Filter to 1981-01-01 through 2025-12-31 in PKT
    start_date = pd.Timestamp("1981-01-01")
    end_date = pd.Timestamp("2025-12-31")
    df_era5_hr_filtered = df_era5_hr[(df_era5_hr["Date"] >= start_date) & (df_era5_hr["Date"] <= end_date)].copy()

    # Aggregate ERA5-Land hourly to daily statistics in PKT
    df_era5_daily = df_era5_hr_filtered.groupby("Date").agg(
        ERA5_Hours_Count=("t2m_C", "count"),
        ERA5_Tmean_C=("t2m_C", "mean"),
        ERA5_Tmax_C=("t2m_C", "max"),
        ERA5_Tmin_C=("t2m_C", "min"),
        ERA5_Tstd_C=("t2m_C", "std")
    ).reset_index()

    df_era5_daily["ERA5_DTR_C"] = df_era5_daily["ERA5_Tmax_C"] - df_era5_daily["ERA5_Tmin_C"]

    print("3. Merging into Harmonized Daily Master Dataset...")
    master_daily = pd.merge(df_power, df_era5_daily, on="Date", how="inner")
    master_daily = master_daily.sort_values("Date").reset_index(drop=True)

    # Calendar metadata columns
    master_daily["Month"] = master_daily["Date"].dt.month
    master_daily["Day"] = master_daily["Date"].dt.day

    def get_season(m):
        if m in [12, 1, 2]:
            return "Winter (DJF)"
        elif m in [3, 4, 5]:
            return "Pre-Monsoon / Spring-Summer (MAM)"
        elif m in [6, 7, 8, 9]:
            return "Summer Monsoon (JJAS)"
        else:
            return "Post-Monsoon / Autumn (ON)"

    master_daily["Season"] = master_daily["Month"].apply(get_season)

    # Ensemble / Consensus metrics and inter-dataset differences (POWER minus ERA5-Land)
    master_daily["Consensus_Tmean_C"] = ((master_daily["POWER_Tmean_C"] + master_daily["ERA5_Tmean_C"]) / 2.0).round(3)
    master_daily["Consensus_Tmax_C"] = ((master_daily["POWER_Tmax_C"] + master_daily["ERA5_Tmax_C"]) / 2.0).round(3)
    master_daily["Consensus_Tmin_C"] = ((master_daily["POWER_Tmin_C"] + master_daily["ERA5_Tmin_C"]) / 2.0).round(3)
    master_daily["Consensus_DTR_C"] = ((master_daily["POWER_DTR_C"] + master_daily["ERA5_DTR_C"]) / 2.0).round(3)

    master_daily["Diff_Tmean_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmean_C"] - master_daily["ERA5_Tmean_C"]).round(3)
    master_daily["Diff_Tmax_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmax_C"] - master_daily["ERA5_Tmax_C"]).round(3)
    master_daily["Diff_Tmin_POWER_minus_ERA5_C"] = (master_daily["POWER_Tmin_C"] - master_daily["ERA5_Tmin_C"]).round(3)

    # Extreme heat flags
    master_daily["ERA5_HotDay_40C"] = (master_daily["ERA5_Tmax_C"] >= 40.0).astype(int)
    master_daily["ERA5_ExtremeHeat_45C"] = (master_daily["ERA5_Tmax_C"] >= 45.0).astype(int)
    master_daily["ERA5_WarmNight_30C"] = (master_daily["ERA5_Tmin_C"] >= 30.0).astype(int)

    master_daily["POWER_HotDay_40C"] = (master_daily["POWER_Tmax_C"] >= 40.0).astype(int)
    master_daily["POWER_ExtremeHeat_45C"] = (master_daily["POWER_Tmax_C"] >= 45.0).astype(int)
    master_daily["POWER_WarmNight_30C"] = (master_daily["POWER_Tmin_C"] >= 30.0).astype(int)

    # Round ERA5 columns for clean CSV storage
    for col in ["ERA5_Tmean_C", "ERA5_Tmax_C", "ERA5_Tmin_C", "ERA5_Tstd_C", "ERA5_DTR_C",
                "POWER_Tmean_C", "POWER_Tmax_C", "POWER_Tmin_C", "POWER_DTR_C"]:
        master_daily[col] = master_daily[col].round(3)

    # Reorder columns logically
    cols_order = [
        "Date", "Year", "Month", "Day", "Day_of_Year", "Season",
        "ERA5_Tmean_C", "ERA5_Tmax_C", "ERA5_Tmin_C", "ERA5_DTR_C", "ERA5_Tstd_C", "ERA5_Hours_Count",
        "POWER_Tmean_C", "POWER_Tmax_C", "POWER_Tmin_C", "POWER_DTR_C",
        "Consensus_Tmean_C", "Consensus_Tmax_C", "Consensus_Tmin_C", "Consensus_DTR_C",
        "Diff_Tmean_POWER_minus_ERA5_C", "Diff_Tmax_POWER_minus_ERA5_C", "Diff_Tmin_POWER_minus_ERA5_C",
        "ERA5_HotDay_40C", "ERA5_ExtremeHeat_45C", "ERA5_WarmNight_30C",
        "POWER_HotDay_40C", "POWER_ExtremeHeat_45C", "POWER_WarmNight_30C"
    ]
    master_daily = master_daily[cols_order]

    daily_out_path = os.path.join(OUTPUT_DIR, "Sukkur_Climate_Master_Daily_1981_2025.csv")
    master_daily.to_csv(daily_out_path, index=False)
    print(f"Saved Daily Master Dataset: {daily_out_path} ({len(master_daily)} rows)")

    print("4. Computing Annual Master Dataset & Extreme Heat Metrics...")
    # 90th percentile baselines (1981-2010 WMO standard reference period)
    ref_mask = (master_daily["Year"] >= 1981) & (master_daily["Year"] <= 2010)
    era5_tx90_thresh = master_daily.loc[ref_mask, "ERA5_Tmax_C"].quantile(0.90)
    era5_tn90_thresh = master_daily.loc[ref_mask, "ERA5_Tmin_C"].quantile(0.90)
    power_tx90_thresh = master_daily.loc[ref_mask, "POWER_Tmax_C"].quantile(0.90)
    power_tn90_thresh = master_daily.loc[ref_mask, "POWER_Tmin_C"].quantile(0.90)

    annual_records = []
    for yr, grp in master_daily.groupby("Year"):
        grp = grp.sort_values("Day_of_Year")
        # Onset and cessation of >= 40C (Summer arrival / retreat)
        e40 = grp[grp["ERA5_Tmax_C"] >= 40.0]
        p40 = grp[grp["POWER_Tmax_C"] >= 40.0]

        era5_first_40_doy = int(e40["Day_of_Year"].min()) if len(e40) > 0 else np.nan
        era5_last_40_doy = int(e40["Day_of_Year"].max()) if len(e40) > 0 else np.nan
        era5_season_span_40 = int(era5_last_40_doy - era5_first_40_doy + 1) if len(e40) > 0 else 0

        power_first_40_doy = int(p40["Day_of_Year"].min()) if len(p40) > 0 else np.nan
        power_last_40_doy = int(p40["Day_of_Year"].max()) if len(p40) > 0 else np.nan
        power_season_span_40 = int(power_last_40_doy - power_first_40_doy + 1) if len(p40) > 0 else 0

        # Sustained onset: first day of first 3-consecutive-day spell >= 40C
        def first_3day_spell_doy(df_sub, col_thresh):
            bools = (df_sub[col_thresh] >= 40.0).values
            doys = df_sub["Day_of_Year"].values
            for idx in range(len(bools) - 2):
                if bools[idx] and bools[idx + 1] and bools[idx + 2]:
                    return int(doys[idx])
            return np.nan

        era5_sust_onset_40 = first_3day_spell_doy(grp, "ERA5_Tmax_C")
        power_sust_onset_40 = first_3day_spell_doy(grp, "POWER_Tmax_C")

        # Heatwave spell metrics
        e_max_spell_40, e_spells3_40, e_spells5_40 = max_spell_length(grp["ERA5_Tmax_C"] >= 40.0)
        e_max_spell_45, e_spells3_45, e_spells5_45 = max_spell_length(grp["ERA5_Tmax_C"] >= 45.0)
        p_max_spell_40, p_spells3_40, p_spells5_40 = max_spell_length(grp["POWER_Tmax_C"] >= 40.0)
        p_max_spell_45, p_spells3_45, p_spells5_45 = max_spell_length(grp["POWER_Tmax_C"] >= 45.0)

        rec = {
            "Year": int(yr),
            "Days_Count": int(len(grp)),
            "ERA5_Tmean_C": round(grp["ERA5_Tmean_C"].mean(), 3),
            "ERA5_Tmax_Mean_C": round(grp["ERA5_Tmax_C"].mean(), 3),
            "ERA5_Tmin_Mean_C": round(grp["ERA5_Tmin_C"].mean(), 3),
            "ERA5_DTR_Mean_C": round(grp["ERA5_DTR_C"].mean(), 3),
            "ERA5_TXx_AbsMax_C": round(grp["ERA5_Tmax_C"].max(), 3),
            "ERA5_TNn_AbsMin_C": round(grp["ERA5_Tmin_C"].min(), 3),
            "ERA5_HotDays_40C": int((grp["ERA5_Tmax_C"] >= 40.0).sum()),
            "ERA5_ExtremeDays_45C": int((grp["ERA5_Tmax_C"] >= 45.0).sum()),
            "ERA5_WarmNights_30C": int((grp["ERA5_Tmin_C"] >= 30.0).sum()),
            "ERA5_TX90p_Days": int((grp["ERA5_Tmax_C"] > era5_tx90_thresh).sum()),
            "ERA5_TN90p_Days": int((grp["ERA5_Tmin_C"] > era5_tn90_thresh).sum()),
            "ERA5_First_40C_DOY": era5_first_40_doy,
            "ERA5_Sustained3d_40C_DOY": era5_sust_onset_40,
            "ERA5_Last_40C_DOY": era5_last_40_doy,
            "ERA5_ThermalSummer_Span_Days": era5_season_span_40,
            "ERA5_MaxSpell_40C_Days": e_max_spell_40,
            "ERA5_Spells3d_40C_Count": e_spells3_40,
            "ERA5_MaxSpell_45C_Days": e_max_spell_45,
            "ERA5_Spells3d_45C_Count": e_spells3_45,

            "POWER_Tmean_C": round(grp["POWER_Tmean_C"].mean(), 3),
            "POWER_Tmax_Mean_C": round(grp["POWER_Tmax_C"].mean(), 3),
            "POWER_Tmin_Mean_C": round(grp["POWER_Tmin_C"].mean(), 3),
            "POWER_DTR_Mean_C": round(grp["POWER_DTR_C"].mean(), 3),
            "POWER_TXx_AbsMax_C": round(grp["POWER_Tmax_C"].max(), 3),
            "POWER_TNn_AbsMin_C": round(grp["POWER_Tmin_C"].min(), 3),
            "POWER_HotDays_40C": int((grp["POWER_Tmax_C"] >= 40.0).sum()),
            "POWER_ExtremeDays_45C": int((grp["POWER_Tmax_C"] >= 45.0).sum()),
            "POWER_WarmNights_30C": int((grp["POWER_Tmin_C"] >= 30.0).sum()),
            "POWER_TX90p_Days": int((grp["POWER_Tmax_C"] > power_tx90_thresh).sum()),
            "POWER_TN90p_Days": int((grp["POWER_Tmin_C"] > power_tn90_thresh).sum()),
            "POWER_First_40C_DOY": power_first_40_doy,
            "POWER_Sustained3d_40C_DOY": power_sust_onset_40,
            "POWER_Last_40C_DOY": power_last_40_doy,
            "POWER_ThermalSummer_Span_Days": power_season_span_40,
            "POWER_MaxSpell_40C_Days": p_max_spell_40,
            "POWER_Spells3d_40C_Count": p_spells3_40,
            "POWER_MaxSpell_45C_Days": p_max_spell_45,
            "POWER_Spells3d_45C_Count": p_spells3_45,

            "Consensus_Tmean_C": round(grp["Consensus_Tmean_C"].mean(), 3),
            "Consensus_Tmax_Mean_C": round(grp["Consensus_Tmax_C"].mean(), 3),
            "Consensus_Tmin_Mean_C": round(grp["Consensus_Tmin_C"].mean(), 3),
        }
        annual_records.append(rec)

    master_annual = pd.DataFrame(annual_records)
    annual_out_path = os.path.join(OUTPUT_DIR, "Sukkur_Climate_Master_Annual_Summary_1981_2025.csv")
    master_annual.to_csv(annual_out_path, index=False)
    print(f"Saved Annual Summary Dataset: {annual_out_path} ({len(master_annual)} rows)")

    # Print key summary diagnostics
    print("\n--- KEY STATISTICAL DIAGNOSTICS (1981-2025) ---")
    years = master_annual["Year"].values
    for col in [
        "ERA5_Tmean_C", "ERA5_Tmax_Mean_C", "ERA5_Tmin_Mean_C", "ERA5_HotDays_40C", "ERA5_ExtremeDays_45C",
        "ERA5_WarmNights_30C", "ERA5_First_40C_DOY", "ERA5_Last_40C_DOY", "ERA5_ThermalSummer_Span_Days",
        "POWER_Tmean_C", "POWER_Tmax_Mean_C", "POWER_Tmin_Mean_C", "POWER_HotDays_40C", "POWER_ExtremeDays_45C",
        "POWER_WarmNights_30C", "POWER_First_40C_DOY", "POWER_Last_40C_DOY", "POWER_ThermalSummer_Span_Days",
        "Consensus_Tmean_C", "Consensus_Tmax_Mean_C", "Consensus_Tmin_Mean_C"
    ]:
        res = sens_slope_and_mk(years, master_annual[col].values)
        mean_val = master_annual[col].mean()
        print(f"{col:30s}: Mean={mean_val:7.2f} | SenSlope={res['slope_dec']:+7.3f}/dec | MK_Z={res['mk_z']:+6.2f} | p={res['mk_p']:.4f}")

    # Inter-dataset comparison stats
    for pair in [("POWER_Tmean_C", "ERA5_Tmean_C"), ("POWER_Tmax_C", "ERA5_Tmax_C"), ("POWER_Tmin_C", "ERA5_Tmin_C")]:
        p_col, e_col = pair
        diff = master_daily[p_col] - master_daily[e_col]
        r, _ = stats.pearsonr(master_daily[p_col], master_daily[e_col])
        rmse = np.sqrt(np.mean(diff ** 2))
        mae = np.mean(np.abs(diff))
        mbe = np.mean(diff)
        print(f"Daily {p_col} vs {e_col}: r={r:.4f}, MBE={mbe:+.3f}C, MAE={mae:.3f}C, RMSE={rmse:.3f}C")


if __name__ == "__main__":
    main()