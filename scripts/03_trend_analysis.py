#!/usr/bin/env python3
"""
Phase 3: Core Statistical Trend Analysis (1981–2025)
Study: Shifting Summer Temperatures, Earlier Spring Heat Arrival, and Rising Warm Nights
       in Sukkur, Pakistan (1981–2025): Insights from ERA5-Land and NASA POWER Reanalyses
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
TABLES_OUT_DIR = REPO_ROOT / "data" / "results"
TABLES_OUT_DIR.mkdir(parents=True, exist_ok=True)

DAILY_MASTER_FILE = PROCESSED_DIR / "Sukkur_Climate_Master_Daily_1981_2025.csv"
ANNUAL_MASTER_FILE = PROCESSED_DIR / "Sukkur_Climate_Master_Annual_Summary_1981_2025.csv"


def mk_and_sens_slope(years, values, scale_factor=10.0):
    mask = ~np.isnan(years) & ~np.isnan(values)
    x = np.asarray(years[mask], dtype=float)
    y = np.asarray(values[mask], dtype=float)
    n = len(y)
    if n < 5:
        return None

    sen_res = stats.theilslopes(y, x, 0.95)
    sen_slope_yr = float(sen_res[0])
    sen_ci_low_yr = float(sen_res[2])
    sen_ci_high_yr = float(sen_res[3])

    s = 0
    for k in range(n - 1):
        s += np.sum(np.sign(y[k + 1:] - y[k]))

    unique, counts = np.unique(y, return_counts=True)
    ties = counts[counts > 1]
    tie_term = np.sum(ties * (ties - 1) * (2 * ties + 5)) if len(ties) > 0 else 0
    var_s = (n * (n - 1) * (2 * n + 5) - tie_term) / 18.0

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    mk_p = 2.0 * (1.0 - stats.norm.cdf(abs(z)))
    tau, _ = stats.kendalltau(x, y)
    ols = stats.linregress(x, y)

    early_mask = (x >= 1981) & (x <= 1995)
    late_mask = (x >= 2011) & (x <= 2025)
    mean_early_15 = float(np.mean(y[early_mask])) if np.any(early_mask) else np.nan
    mean_late_15 = float(np.mean(y[late_mask])) if np.any(late_mask) else np.nan
    delta_15yr = mean_late_15 - mean_early_15

    if mk_p < 0.001:
        sig_flag = "*** (p < 0.001)"
    elif mk_p < 0.01:
        sig_flag = "** (p < 0.01)"
    elif mk_p < 0.05:
        sig_flag = "* (p < 0.05)"
    elif mk_p < 0.10:
        sig_flag = "+ (p < 0.10)"
    else:
        sig_flag = "ns (p >= 0.10)"

    return {
        "N_Years": int(n),
        "Mean_1981_2025": round(float(np.mean(y)), 3),
        "SD_1981_2025": round(float(np.std(y, ddof=1)), 3),
        "Min_Val": round(float(np.min(y)), 3),
        "Min_Year": int(x[np.argmin(y)]),
        "Max_Val": round(float(np.max(y)), 3),
        "Max_Year": int(x[np.argmax(y)]),
        "Mean_1981_1995": round(mean_early_15, 3),
        "Mean_2011_2025": round(mean_late_15, 3),
        "Regime_Shift_Late_minus_Early": round(delta_15yr, 3),
        "Sen_Slope_per_Decade": round(sen_slope_yr * scale_factor, 4),
        "Sen_CI95_Low_per_Decade": round(sen_ci_low_yr * scale_factor, 4),
        "Sen_CI95_High_per_Decade": round(sen_ci_high_yr * scale_factor, 4),
        "Total_45yr_Sen_Change": round(sen_slope_yr * 44.0, 3),
        "OLS_Slope_per_Decade": round(float(ols.slope) * scale_factor, 4),
        "OLS_p_value": round(float(ols.pvalue), 5),
        "MK_S": int(s),
        "MK_Z": round(float(z), 4),
        "MK_p_value": round(float(mk_p), 5),
        "Kendall_Tau": round(float(tau), 4),
        "Significance": sig_flag,
    }


def main():
    df_daily = pd.read_csv(DAILY_MASTER_FILE, parse_dates=["Date"])
    df_annual = pd.read_csv(ANNUAL_MASTER_FILE)

    extra_annual = []
    for yr, grp in df_daily.groupby("Year"):
        grp = grp.sort_values("Day_of_Year")
        e38 = grp[grp["ERA5_Tmax_C"] >= 38.0]
        p38 = grp[grp["POWER_Tmax_C"] >= 38.0]
        c38 = grp[grp["Consensus_Tmax_C"] >= 38.0]
        c40 = grp[grp["Consensus_Tmax_C"] >= 40.0]

        m_jjas = grp[grp["Month"].isin([6, 7, 8, 9])]
        m_mam = grp[grp["Month"].isin([3, 4, 5])]
        m_djf = grp[grp["Month"].isin([12, 1, 2])]
        m_on = grp[grp["Month"].isin([10, 11])]
        m_may_aug = grp[grp["Month"].isin([5, 6, 7, 8])]
        m_mar_apr = grp[grp["Month"].isin([3, 4])]

        extra_annual.append({
            "Year": int(yr),
            "ERA5_First_38C_DOY": int(e38["Day_of_Year"].min()) if len(e38) > 0 else np.nan,
            "POWER_First_38C_DOY": int(p38["Day_of_Year"].min()) if len(p38) > 0 else np.nan,
            "Consensus_First_38C_DOY": int(c38["Day_of_Year"].min()) if len(c38) > 0 else np.nan,
            "Consensus_First_40C_DOY": int(c40["Day_of_Year"].min()) if len(c40) > 0 else np.nan,
            "Consensus_Last_40C_DOY": int(c40["Day_of_Year"].max()) if len(c40) > 0 else np.nan,
            "Consensus_HotDays_40C": int((grp["Consensus_Tmax_C"] >= 40.0).sum()),
            "Consensus_ExtremeDays_45C": int((grp["Consensus_Tmax_C"] >= 45.0).sum()),
            "Consensus_WarmNights_30C": int((grp["Consensus_Tmin_C"] >= 30.0).sum()),
            "Consensus_DTR_Mean_C": round(float(grp["Consensus_DTR_C"].mean()), 3),
            "ERA5_SummerMayAug_Tmean_C": round(float(m_may_aug["ERA5_Tmean_C"].mean()), 3),
            "ERA5_SummerMayAug_Tmax_C": round(float(m_may_aug["ERA5_Tmax_C"].mean()), 3),
            "ERA5_SummerMayAug_Tmin_C": round(float(m_may_aug["ERA5_Tmin_C"].mean()), 3),
            "POWER_SummerMayAug_Tmean_C": round(float(m_may_aug["POWER_Tmean_C"].mean()), 3),
            "POWER_SummerMayAug_Tmax_C": round(float(m_may_aug["POWER_Tmax_C"].mean()), 3),
            "POWER_SummerMayAug_Tmin_C": round(float(m_may_aug["POWER_Tmin_C"].mean()), 3),
            "Consensus_SummerMayAug_Tmean_C": round(float(m_may_aug["Consensus_Tmean_C"].mean()), 3),
            "Consensus_SummerMayAug_Tmax_C": round(float(m_may_aug["Consensus_Tmax_C"].mean()), 3),
            "Consensus_SummerMayAug_Tmin_C": round(float(m_may_aug["Consensus_Tmin_C"].mean()), 3),
            "ERA5_SpringMarApr_Tmax_C": round(float(m_mar_apr["ERA5_Tmax_C"].mean()), 3),
            "ERA5_SpringMarApr_Tmin_C": round(float(m_mar_apr["ERA5_Tmin_C"].mean()), 3),
            "POWER_SpringMarApr_Tmax_C": round(float(m_mar_apr["POWER_Tmax_C"].mean()), 3),
            "POWER_SpringMarApr_Tmin_C": round(float(m_mar_apr["POWER_Tmin_C"].mean()), 3),
            "Consensus_SpringMarApr_Tmax_C": round(float(m_mar_apr["Consensus_Tmax_C"].mean()), 3),
            "Consensus_SpringMarApr_Tmin_C": round(float(m_mar_apr["Consensus_Tmin_C"].mean()), 3),
            "ERA5_DJF_Tmean_C": round(float(m_djf["ERA5_Tmean_C"].mean()), 3),
            "POWER_DJF_Tmean_C": round(float(m_djf["POWER_Tmean_C"].mean()), 3),
            "Consensus_DJF_Tmean_C": round(float(m_djf["Consensus_Tmean_C"].mean()), 3),
            "ERA5_MAM_Tmean_C": round(float(m_mam["ERA5_Tmean_C"].mean()), 3),
            "ERA5_MAM_Tmin_C": round(float(m_mam["ERA5_Tmin_C"].mean()), 3),
            "POWER_MAM_Tmean_C": round(float(m_mam["POWER_Tmean_C"].mean()), 3),
            "POWER_MAM_Tmin_C": round(float(m_mam["POWER_Tmin_C"].mean()), 3),
            "Consensus_MAM_Tmean_C": round(float(m_mam["Consensus_Tmean_C"].mean()), 3),
            "Consensus_MAM_Tmin_C": round(float(m_mam["Consensus_Tmin_C"].mean()), 3),
            "ERA5_JJAS_Tmean_C": round(float(m_jjas["ERA5_Tmean_C"].mean()), 3),
            "ERA5_JJAS_Tmax_C": round(float(m_jjas["ERA5_Tmax_C"].mean()), 3),
            "ERA5_JJAS_Tmin_C": round(float(m_jjas["ERA5_Tmin_C"].mean()), 3),
            "POWER_JJAS_Tmean_C": round(float(m_jjas["POWER_Tmean_C"].mean()), 3),
            "POWER_JJAS_Tmax_C": round(float(m_jjas["POWER_Tmax_C"].mean()), 3),
            "POWER_JJAS_Tmin_C": round(float(m_jjas["POWER_Tmin_C"].mean()), 3),
            "Consensus_JJAS_Tmean_C": round(float(m_jjas["Consensus_Tmean_C"].mean()), 3),
            "Consensus_JJAS_Tmax_C": round(float(m_jjas["Consensus_Tmax_C"].mean()), 3),
            "Consensus_JJAS_Tmin_C": round(float(m_jjas["Consensus_Tmin_C"].mean()), 3),
            "ERA5_ON_Tmean_C": round(float(m_on["ERA5_Tmean_C"].mean()), 3),
            "POWER_ON_Tmean_C": round(float(m_on["POWER_Tmean_C"].mean()), 3),
            "Consensus_ON_Tmean_C": round(float(m_on["Consensus_Tmean_C"].mean()), 3),
        })

    df_extra = pd.DataFrame(extra_annual)
    df_full_ann = pd.merge(df_annual, df_extra, on="Year", how="inner")

    metrics_spec = [
        ("Annual Temperatures", "ERA5-Land (~9 km)", "ERA5_Tmean_C", "Annual Mean Tmean (°C)", "°C/decade"),
        ("Annual Temperatures", "NASA POWER (~50 km)", "POWER_Tmean_C", "Annual Mean Tmean (°C)", "°C/decade"),
        ("Annual Temperatures", "Consensus Ensemble", "Consensus_Tmean_C", "Annual Mean Tmean (°C)", "°C/decade"),
        ("Annual Temperatures", "ERA5-Land (~9 km)", "ERA5_Tmax_Mean_C", "Annual Mean Tmax (°C)", "°C/decade"),
        ("Annual Temperatures", "NASA POWER (~50 km)", "POWER_Tmax_Mean_C", "Annual Mean Tmax (°C)", "°C/decade"),
        ("Annual Temperatures", "Consensus Ensemble", "Consensus_Tmax_Mean_C", "Annual Mean Tmax (°C)", "°C/decade"),
        ("Annual Temperatures", "ERA5-Land (~9 km)", "ERA5_Tmin_Mean_C", "Annual Mean Tmin (°C)", "°C/decade"),
        ("Annual Temperatures", "NASA POWER (~50 km)", "POWER_Tmin_Mean_C", "Annual Mean Tmin (°C)", "°C/decade"),
        ("Annual Temperatures", "Consensus Ensemble", "Consensus_Tmin_Mean_C", "Annual Mean Tmin (°C)", "°C/decade"),
        ("Annual Temperatures", "ERA5-Land (~9 km)", "ERA5_DTR_Mean_C", "Annual Mean DTR (°C)", "°C/decade"),
        ("Annual Temperatures", "NASA POWER (~50 km)", "POWER_DTR_Mean_C", "Annual Mean DTR (°C)", "°C/decade"),
        ("Annual Temperatures", "Consensus Ensemble", "Consensus_DTR_Mean_C", "Annual Mean DTR (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_SummerMayAug_Tmean_C", "Core Summer (May–Aug) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_SummerMayAug_Tmean_C", "Core Summer (May–Aug) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_SummerMayAug_Tmean_C", "Core Summer (May–Aug) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_SummerMayAug_Tmax_C", "Core Summer (May–Aug) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_SummerMayAug_Tmax_C", "Core Summer (May–Aug) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_SummerMayAug_Tmax_C", "Core Summer (May–Aug) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_SummerMayAug_Tmin_C", "Core Summer (May–Aug) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_SummerMayAug_Tmin_C", "Core Summer (May–Aug) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_SummerMayAug_Tmin_C", "Core Summer (May–Aug) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_SpringMarApr_Tmax_C", "Spring Onset (Mar–Apr) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_SpringMarApr_Tmax_C", "Spring Onset (Mar–Apr) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_SpringMarApr_Tmax_C", "Spring Onset (Mar–Apr) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_SpringMarApr_Tmin_C", "Spring Onset (Mar–Apr) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_SpringMarApr_Tmin_C", "Spring Onset (Mar–Apr) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_SpringMarApr_Tmin_C", "Spring Onset (Mar–Apr) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_MAM_Tmean_C", "Pre-Monsoon (MAM) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_MAM_Tmean_C", "Pre-Monsoon (MAM) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_MAM_Tmean_C", "Pre-Monsoon (MAM) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_MAM_Tmin_C", "Pre-Monsoon (MAM) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_MAM_Tmin_C", "Pre-Monsoon (MAM) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_MAM_Tmin_C", "Pre-Monsoon (MAM) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_JJAS_Tmean_C", "Summer Monsoon (JJAS) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_JJAS_Tmean_C", "Summer Monsoon (JJAS) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_JJAS_Tmean_C", "Summer Monsoon (JJAS) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_JJAS_Tmax_C", "Summer Monsoon (JJAS) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_JJAS_Tmax_C", "Summer Monsoon (JJAS) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_JJAS_Tmax_C", "Summer Monsoon (JJAS) Tmax (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_JJAS_Tmin_C", "Summer Monsoon (JJAS) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_JJAS_Tmin_C", "Summer Monsoon (JJAS) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_JJAS_Tmin_C", "Summer Monsoon (JJAS) Tmin (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_DJF_Tmean_C", "Winter (DJF) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_DJF_Tmean_C", "Winter (DJF) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_DJF_Tmean_C", "Winter (DJF) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "ERA5-Land (~9 km)", "ERA5_ON_Tmean_C", "Post-Monsoon (ON) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "NASA POWER (~50 km)", "POWER_ON_Tmean_C", "Post-Monsoon (ON) Tmean (°C)", "°C/decade"),
        ("Seasonal Temperatures", "Consensus Ensemble", "Consensus_ON_Tmean_C", "Post-Monsoon (ON) Tmean (°C)", "°C/decade"),
        ("Angle 1: Heat Arrival & Timing", "ERA5-Land (~9 km)", "ERA5_First_38C_DOY", "First Tmax >= 38°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "NASA POWER (~50 km)", "POWER_First_38C_DOY", "First Tmax >= 38°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "Consensus Ensemble", "Consensus_First_38C_DOY", "First Tmax >= 38°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "ERA5-Land (~9 km)", "ERA5_First_40C_DOY", "First Tmax >= 40°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "NASA POWER (~50 km)", "POWER_First_40C_DOY", "First Tmax >= 40°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "Consensus Ensemble", "Consensus_First_40C_DOY", "First Tmax >= 40°C Arrival (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "ERA5-Land (~9 km)", "ERA5_Sustained3d_40C_DOY", "Sustained 3-Day >= 40°C Onset (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "NASA POWER (~50 km)", "POWER_Sustained3d_40C_DOY", "Sustained 3-Day >= 40°C Onset (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "ERA5-Land (~9 km)", "ERA5_Last_40C_DOY", "Last Tmax >= 40°C Retreat (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "NASA POWER (~50 km)", "POWER_Last_40C_DOY", "Last Tmax >= 40°C Retreat (DOY)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "ERA5-Land (~9 km)", "ERA5_ThermalSummer_Span_Days", "Thermal Summer Span (Days between 1st & Last >= 40°C)", "days/decade"),
        ("Angle 1: Heat Arrival & Timing", "NASA POWER (~50 km)", "POWER_ThermalSummer_Span_Days", "Thermal Summer Span (Days between 1st & Last >= 40°C)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_HotDays_40C", "Hot Days (Tmax >= 40°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_HotDays_40C", "Hot Days (Tmax >= 40°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "Consensus Ensemble", "Consensus_HotDays_40C", "Hot Days (Tmax >= 40°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_ExtremeDays_45C", "Extreme Heat Days (Tmax >= 45°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_ExtremeDays_45C", "Extreme Heat Days (Tmax >= 45°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "Consensus Ensemble", "Consensus_ExtremeDays_45C", "Extreme Heat Days (Tmax >= 45°C, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_WarmNights_30C", "Tropical Warm Nights (Tmin >= 30°C, nights/yr)", "nights/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_WarmNights_30C", "Tropical Warm Nights (Tmin >= 30°C, nights/yr)", "nights/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "Consensus Ensemble", "Consensus_WarmNights_30C", "Tropical Warm Nights (Tmin >= 30°C, nights/yr)", "nights/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_TX90p_Days", "Warm Days TX90p (>90th Pct Tmax, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_TX90p_Days", "Warm Days TX90p (>90th Pct Tmax, days/yr)", "days/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_TN90p_Days", "Warm Nights TN90p (>90th Pct Tmin, nights/yr)", "nights/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_TN90p_Days", "Warm Nights TN90p (>90th Pct Tmin, nights/yr)", "nights/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_TXx_AbsMax_C", "Annual Absolute Maximum Tmax (TXx, °C)", "°C/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_TXx_AbsMax_C", "Annual Absolute Maximum Tmax (TXx, °C)", "°C/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "ERA5-Land (~9 km)", "ERA5_TNn_AbsMin_C", "Annual Absolute Minimum Tmin (TNn, °C)", "°C/decade"),
        ("Angle 2: Extreme Heat & Warm Nights", "NASA POWER (~50 km)", "POWER_TNn_AbsMin_C", "Annual Absolute Minimum Tmin (TNn, °C)", "°C/decade"),
    ]

    years_arr = df_full_ann["Year"].values
    trend_rows = []
    for cat, ds_name, col_name, label, unit in metrics_spec:
        res = mk_and_sens_slope(years_arr, df_full_ann[col_name].values, scale_factor=10.0)
        if res is not None:
            row = {"Category": cat, "Metric": label, "Dataset": ds_name, "Column_Code": col_name, "Slope_Unit": unit}
            row.update(res)
            trend_rows.append(row)

    pd.DataFrame(trend_rows).to_csv(TABLES_OUT_DIR / "Sukkur_Phase3_Core_Trend_Analysis_1981_2025.csv", index=False)

    month_names = {1: "January", 2: "February", 3: "March", 4: "April", 5: "May", 6: "June",
                   7: "July", 8: "August", 9: "September", 10: "October", 11: "November", 12: "December"}
    df_monthly = df_daily.groupby(["Year", "Month"], as_index=False).agg(
        ERA5_Tmean_C=("ERA5_Tmean_C", "mean"), ERA5_Tmax_C=("ERA5_Tmax_C", "mean"),
        ERA5_Tmin_C=("ERA5_Tmin_C", "mean"), ERA5_DTR_C=("ERA5_DTR_C", "mean"),
        POWER_Tmean_C=("POWER_Tmean_C", "mean"), POWER_Tmax_C=("POWER_Tmax_C", "mean"),
        POWER_Tmin_C=("POWER_Tmin_C", "mean"), POWER_DTR_C=("POWER_DTR_C", "mean"),
        Consensus_Tmean_C=("Consensus_Tmean_C", "mean"), Consensus_Tmax_C=("Consensus_Tmax_C", "mean"),
        Consensus_Tmin_C=("Consensus_Tmin_C", "mean"), Consensus_DTR_C=("Consensus_DTR_C", "mean"),
    )

    monthly_trend_rows = []
    for m in range(1, 13):
        sub_m = df_monthly[df_monthly["Month"] == m].sort_values("Year")
        yrs_m = sub_m["Year"].values
        for ds_label, prefix in [("ERA5-Land (~9 km)", "ERA5"), ("NASA POWER (~50 km)", "POWER"), ("Consensus Ensemble", "Consensus")]:
            for var_short in ["Tmean", "Tmax", "Tmin", "DTR"]:
                col_c = f"{prefix}_{var_short}_C"
                res_m = mk_and_sens_slope(yrs_m, sub_m[col_c].values, scale_factor=10.0)
                monthly_trend_rows.append({
                    "Month_Num": m, "Month_Name": month_names[m], "Dataset": ds_label, "Variable": var_short,
                    "Mean_1981_2025_C": res_m["Mean_1981_2025"], "Mean_1981_1995_C": res_m["Mean_1981_1995"],
                    "Mean_2011_2025_C": res_m["Mean_2011_2025"], "Regime_Shift_C": res_m["Regime_Shift_Late_minus_Early"],
                    "Sen_Slope_C_per_Decade": res_m["Sen_Slope_per_Decade"], "Sen_CI95_Low": res_m["Sen_CI95_Low_per_Decade"],
                    "Sen_CI95_High": res_m["Sen_CI95_High_per_Decade"], "Total_45yr_Change_C": res_m["Total_45yr_Sen_Change"],
                    "MK_Z": res_m["MK_Z"], "MK_p_value": res_m["MK_p_value"], "Kendall_Tau": res_m["Kendall_Tau"],
                    "Significance": res_m["Significance"],
                })

    pd.DataFrame(monthly_trend_rows).to_csv(TABLES_OUT_DIR / "Sukkur_Phase3_Monthly_Trends_1981_2025.csv", index=False)

if __name__ == "__main__":
    main()
