#!/usr/bin/env python3
"""
Phase 4: Plotting the Final Figures (1981–2025)

Description: This script generates the 8 final charts for my paper. 
I set them up to export at high resolution (300 DPI) so they are ready for publication. 
They export in the exact order they appear in the manuscript:

- Fig1_CrossValidation_Daily_Scatter_Density.png
- Fig2_Monthly_Climatology_Microclimate_Offset.png
- Fig3_Annual_Tmin_Warming_Trends_1981_2025.png
- Fig4_Annual_Tmax_Evaporative_Buffering_1981_2025.png
- Fig5_Diurnal_Temperature_Range_Compression_1981_2025.png
- Fig6_Monthly_Warming_Asymmetry_Profile.png
- Fig7_Spring_Heat_Onset_DOY_1981_2025.png
- Fig8_Tropical_Warm_Nights_Surge_1981_2025.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.titlesize": 12.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 11,
    "axes.labelweight": "bold",
    "xtick.labelsize": 9.8,
    "ytick.labelsize": 9.8,
    "legend.fontsize": 9.2,
    "figure.titlesize": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.28,
    "grid.linestyle": "--",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.22,
})

REPO_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
RESULTS_DIR = REPO_ROOT / "data" / "results"
OUT_DIR = REPO_ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DAILY_FILE = PROCESSED_DIR / "Sukkur_Climate_Master_Daily_1981_2025.csv"
ANNUAL_FILE = PROCESSED_DIR / "Sukkur_Climate_Master_Annual_Summary_1981_2025.csv"
MONTHLY_OFFSET_FILE = RESULTS_DIR / "Sukkur_Phase2_Monthly_Microclimate_Offset_1981_2025.csv"
MONTHLY_TRENDS_FILE = RESULTS_DIR / "Sukkur_Phase3_Monthly_Trends_1981_2025.csv"

COL_ERA5 = "#0072B2"
COL_POWER = "#D55E00"
COL_CONS = "#111111"


def sen_line(years, vals):
    mask = ~np.isnan(years) & ~np.isnan(vals)
    x = np.asarray(years[mask], dtype=float)
    y = np.asarray(vals[mask], dtype=float)
    res = stats.theilslopes(y, x, 0.95)
    slope = float(res[0])
    intercept = float(res[1])
    return slope, intercept, intercept + slope * years


def main():
    df_daily = pd.read_csv(DAILY_FILE, parse_dates=["Date"])
    df_annual = pd.read_csv(ANNUAL_FILE)
    df_m_off = pd.read_csv(MONTHLY_OFFSET_FILE)
    df_m_trd = pd.read_csv(MONTHLY_TRENDS_FILE)

    df_annual["Consensus_DTR_Mean_C"] = ((df_annual["ERA5_DTR_Mean_C"] + df_annual["POWER_DTR_Mean_C"]) / 2.0).round(3)
    c_daily = df_daily.groupby("Year").agg(
        Consensus_WarmNights_30C=("Consensus_Tmin_C", lambda s: int((s >= 30.0).sum())),
        Consensus_First_40C_DOY=("Day_of_Year", lambda s: int(s[df_daily.loc[s.index, "Consensus_Tmax_C"] >= 40.0].min()) if (df_daily.loc[s.index, "Consensus_Tmax_C"] >= 40.0).any() else np.nan)
    ).reset_index()
    df_annual = pd.merge(df_annual, c_daily, on="Year", how="left")
    years = df_annual["Year"].values

    # FIGURE 1
    fig, ax = plt.subplots(figsize=(9.2, 7.2))
    ax.scatter(df_daily["ERA5_Tmax_C"], df_daily["POWER_Tmax_C"], s=8, alpha=0.12, color="#D55E00", label="_nolegend_", rasterized=True)
    ax.scatter(df_daily["ERA5_Tmean_C"], df_daily["POWER_Tmean_C"], s=8, alpha=0.12, color="#009E73", label="_nolegend_", rasterized=True)
    ax.scatter(df_daily["ERA5_Tmin_C"], df_daily["POWER_Tmin_C"], s=8, alpha=0.12, color="#0072B2", label="_nolegend_", rasterized=True)
    lims = [0, 53]
    ax.plot(lims, lims, "k--", lw=1.6, alpha=0.75, label="1:1 Exact Agreement Line ($y = x$)")
    for e_col, p_col, col_hex, name_lbl in [
        ("ERA5_Tmax_C", "POWER_Tmax_C", "#B22222", "$T_{\\max}$ ($r = 0.982, R^2 = 0.965$, MBE = +3.07°C)"),
        ("ERA5_Tmean_C", "POWER_Tmean_C", "#006D5B", "$T_{\\mathrm{mean}}$ ($r = 0.988, R^2 = 0.976$, MBE = +1.41°C)"),
        ("ERA5_Tmin_C", "POWER_Tmin_C", "#004C99", "$T_{\\min}$ ($r = 0.975, R^2 = 0.951$, MBE = +0.23°C)"),
    ]:
        ols = stats.linregress(df_daily[e_col], df_daily[p_col])
        x_sub = np.linspace(df_daily[e_col].min(), df_daily[e_col].max(), 100)
        ax.plot(x_sub, ols.intercept + ols.slope * x_sub, color=col_hex, lw=2.4, label=name_lbl)
    ax.annotate(
        "Daytime Desert-Grid Warm Offset ($T_{\\max}$):\nMBE = +3.07°C (Centered RMSE = 1.48°C)\nERA5-Land (9 km) resolves Indus River & canal cooling",
        xy=(39.5, 43.2), xytext=(3.0, 32.5),
        arrowprops=dict(arrowstyle="->", color="#B22222", lw=1.5, connectionstyle="arc3,rad=-0.12"),
        fontsize=9.2, bbox=dict(boxstyle="round,pad=0.4", fc="#FFF5F0", ec="#D55E00", lw=1.0)
    )
    ax.annotate(
        "Nighttime Convergence ($T_{\\min}$):\nNear-zero bias (MBE = +0.23°C)\nWillmott $d = 0.986$ (Daily), $0.996$ (Monthly)",
        xy=(16.0, 16.0), xytext=(24.5, 6.5),
        arrowprops=dict(arrowstyle="->", color="#004C99", lw=1.5),
        fontsize=9.2, bbox=dict(boxstyle="round,pad=0.4", fc="#F0F8FF", ec="#0072B2", lw=1.0)
    )
    ax.set_xlim(0, 53)
    ax.set_ylim(0, 56)
    ax.set_xlabel("ECMWF ERA5-Land (~9 km Grid) Daily 2m Temperature (°C)")
    ax.set_ylabel("NASA POWER / MERRA-2 (~50 km Grid) Daily 2m Temperature (°C)", labelpad=8)
    ax.set_title("Daily Cross-Validation of 2m Air Temperatures over Sukkur (1981–2025, N = 16,436 Days)", pad=14)
    ax.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.96, fontsize=9.2)
    fig.savefig(OUT_DIR / "Fig1_CrossValidation_Daily_Scatter_Density.png")
    plt.close(fig)

    # FIGURE 2
    fig, ax1 = plt.subplots(figsize=(9.6, 6.4))
    months_x = df_m_off["Month_Num"].values
    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    ax1.fill_between(months_x, df_m_off["ERA5_Tmax_C"], df_m_off["POWER_Tmax_C"],
                     color="#FDB863", alpha=0.35, label="Daytime Evaporative Cooling Offset ($\\Delta T_{\\max}$: +2.09 to +4.09°C)")
    ax1.plot(months_x, df_m_off["POWER_Tmax_C"], "o-", color=COL_POWER, lw=2.4, ms=6.5, label="NASA POWER $T_{\\max}$ (~50 km Regional Desert Grid)")
    ax1.plot(months_x, df_m_off["ERA5_Tmax_C"], "s-", color=COL_ERA5, lw=2.4, ms=6.5, label="ERA5-Land $T_{\\max}$ (~9 km Indus River / Barrage Grid)")
    ax1.plot(months_x, df_m_off["POWER_Tmin_C"], "o--", color="#E66101", lw=1.9, ms=5.5, alpha=0.85, label="NASA POWER $T_{\\min}$ (~50 km Grid)")
    ax1.plot(months_x, df_m_off["ERA5_Tmin_C"], "s--", color="#2B83BA", lw=1.9, ms=5.5, alpha=0.85, label="ERA5-Land $T_{\\min}$ (~9 km Grid)")
    for m_idx, row in df_m_off.iterrows():
        m_num = row["Month_Num"]
        mid_y = (row["ERA5_Tmax_C"] + row["POWER_Tmax_C"]) / 2.0
        ax1.text(m_num, mid_y, f"+{row['MBE_Tmax_C']:.2f}°C", ha="center", va="center",
                 fontsize=8.2, fontweight="bold", color="#7F2704",
                 bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.80))
    ax1.annotate(
        "Peak Kharif Irrigation & Monsoon Discharge (Jul–Sep):\nMaximum latent heat cooling in 9 km Indus grid\n(Aug $\\Delta T_{\\max} = +4.09^\\circ$C, Sep $\\Delta T_{\\max} = +3.87^\\circ$C)",
        xy=(8, 39.8), xytext=(4.3, 11.2),
        arrowprops=dict(arrowstyle="->", color="#7F2704", lw=1.5, connectionstyle="arc3,rad=0.15"),
        fontsize=9.2, bbox=dict(boxstyle="round,pad=0.4", fc="#FFF8EB", ec="#D55E00", lw=1.0)
    )
    ax1.set_xticks(months_x)
    ax1.set_xticklabels(month_labels)
    ax1.set_ylim(6, 54.5)
    ax1.set_xlabel("Month of Year (1981–2025 Climatological Mean)")
    ax1.set_ylabel("2-Meter Air Temperature (°C)", labelpad=8)
    ax1.set_title("Annual Cycle & Indus Oasis Evaporative Cooling Offset over Sukkur (1981–2025)", pad=14)
    ax1.legend(loc="upper left", ncol=2, frameon=True, facecolor="white", framealpha=0.96, fontsize=8.8)
    fig.savefig(OUT_DIR / "Fig2_Monthly_Climatology_Microclimate_Offset.png")
    plt.close(fig)

    # FIGURE 3
    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    s_p, i_p, fit_p = sen_line(years, df_annual["POWER_Tmin_Mean_C"].values)
    s_e, i_e, fit_e = sen_line(years, df_annual["ERA5_Tmin_Mean_C"].values)
    s_c, i_c, fit_c = sen_line(years, df_annual["Consensus_Tmin_Mean_C"].values)
    ax.plot(years, df_annual["POWER_Tmin_Mean_C"], "o-", color=COL_POWER, lw=1.5, ms=4.5, alpha=0.8, label="NASA POWER Annual $T_{\\min}$ (Mean: 21.15°C)")
    ax.plot(years, fit_p, "--", color=COL_POWER, lw=2.3, label=f"NASA POWER Sen's Slope: +{s_p*10:.3f}°C/decade ($Z = +3.37, p = 0.0008$, ***)")
    ax.plot(years, df_annual["ERA5_Tmin_Mean_C"], "s-", color=COL_ERA5, lw=1.5, ms=4.5, alpha=0.8, label="ERA5-Land Annual $T_{\\min}$ (Mean: 20.92°C)")
    ax.plot(years, fit_e, "--", color=COL_ERA5, lw=2.0, label=f"ERA5-Land Sen's Slope: +{s_e*10:.3f}°C/decade ($Z = +1.48, p = 0.140$, ns)")
    ax.plot(years, df_annual["Consensus_Tmin_Mean_C"], "-", color=COL_CONS, lw=2.4, label=f"Consensus Ensemble $T_{{\\min}}$ (+{s_c*10:.3f}°C/decade, $Z = +2.68, p = 0.0073$, **)")
    ax.axvspan(2011, 2025, color="#FEE0D2", alpha=0.32, label="Recent Warm Regime (2011–2025: Consensus 21.14°C vs. 1981–1995: 20.83°C)")
    ax.set_xlim(1980, 2026)
    ax.set_ylim(19.8, 22.95)
    ax.set_xlabel("Year (1981–2025)")
    ax.set_ylabel("Annual Mean Minimum Temperature $T_{\\min}$ (°C)", labelpad=8)
    ax.set_title("Sustained Nocturnal Warming ($T_{\\min}$) over Sukkur (1981–2025)", pad=14)
    ax.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.96, fontsize=8.8)
    fig.savefig(OUT_DIR / "Fig3_Annual_Tmin_Warming_Trends_1981_2025.png")
    plt.close(fig)

    # FIGURE 4
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    s_p, i_p, fit_p = sen_line(years, df_annual["POWER_Tmax_Mean_C"].values)
    s_e, i_e, fit_e = sen_line(years, df_annual["ERA5_Tmax_Mean_C"].values)
    s_c, i_c, fit_c = sen_line(years, df_annual["Consensus_Tmax_Mean_C"].values)
    ax.plot(years, df_annual["POWER_Tmax_Mean_C"], "o-", color=COL_POWER, lw=1.5, ms=4.5, alpha=0.8, label="NASA POWER Annual $T_{\\max}$ (Mean: 37.08°C)")
    ax.plot(years, fit_p, "--", color=COL_POWER, lw=2.2, label=f"NASA POWER Sen's Slope: {s_p*10:+.3f}°C/decade ($Z = -1.32, p = 0.187$, ns)")
    ax.plot(years, df_annual["Consensus_Tmax_Mean_C"], "-", color=COL_CONS, lw=2.3, label=f"Consensus Ensemble $T_{{\\max}}$ (Mean: 35.55°C | Slope: {s_c*10:+.3f}°C/decade, $p = 0.035$, *)")
    ax.plot(years, df_annual["ERA5_Tmax_Mean_C"], "s-", color=COL_ERA5, lw=1.5, ms=4.5, alpha=0.8, label="ERA5-Land Annual $T_{\\max}$ (Mean: 34.02°C)")
    ax.plot(years, fit_e, "--", color=COL_ERA5, lw=2.2, label=f"ERA5-Land Sen's Slope: {s_e*10:+.3f}°C/decade ($Z = -2.53, p = 0.011$, *)")
    ax.annotate(
        "2022 Mega-Flood & Monsoon Anomaly:\nRecord July–August rainfall suppressed\nannual $T_{\\max}$ (35.00°C POWER / 33.24°C ERA5)",
        xy=(2022, 33.24), xytext=(2004, 31.95),
        arrowprops=dict(arrowstyle="->", color="#0072B2", lw=1.5),
        fontsize=9.0, bbox=dict(boxstyle="round,pad=0.35", fc="#F0F8FF", ec="#0072B2", lw=1.0)
    )
    ax.set_xlim(1980, 2026)
    ax.set_ylim(31.6, 40.5)
    ax.set_xlabel("Year (1981–2025)")
    ax.set_ylabel("Annual Mean Maximum Temperature $T_{\\max}$ (°C)", labelpad=8)
    ax.set_title("Daytime Maximum Temperature ($T_{\\max}$) & Irrigation Evaporative Buffering (1981–2025)", pad=14)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.96, fontsize=8.8)
    fig.savefig(OUT_DIR / "Fig4_Annual_Tmax_Evaporative_Buffering_1981_2025.png")
    plt.close(fig)

    # FIGURE 5
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    s_p, i_p, fit_p = sen_line(years, df_annual["POWER_DTR_Mean_C"].values)
    s_e, i_e, fit_e = sen_line(years, df_annual["ERA5_DTR_Mean_C"].values)
    s_c, i_c, fit_c = sen_line(years, df_annual["Consensus_DTR_Mean_C"].values)
    ax.plot(years, df_annual["POWER_DTR_Mean_C"], "o-", color=COL_POWER, lw=1.6, ms=4.5, alpha=0.85, label="NASA POWER Annual DTR (Mean: 15.94°C)")
    ax.plot(years, fit_p, "--", color=COL_POWER, lw=2.3, label=f"NASA POWER Sen's Slope: {s_p*10:+.3f}°C/decade ($Z = -4.49, p = 0.00001$, ***)")
    ax.plot(years, df_annual["Consensus_DTR_Mean_C"], "-", color=COL_CONS, lw=2.4, label=f"Consensus Ensemble DTR (Mean: 14.52°C | Slope: {s_c*10:+.3f}°C/decade, $p = 0.00002$, ***)")
    ax.plot(years, df_annual["ERA5_DTR_Mean_C"], "s-", color=COL_ERA5, lw=1.6, ms=4.5, alpha=0.85, label="ERA5-Land Annual DTR (Mean: 13.10°C)")
    ax.plot(years, fit_e, "--", color=COL_ERA5, lw=2.3, label=f"ERA5-Land Sen's Slope: {s_e*10:+.3f}°C/decade ($Z = -4.15, p = 0.00003$, ***)")
    ax.annotate(
        "Strong Diurnal Asymmetry Signal ($p < 0.0001$):\nDTR narrowed by $-1.12^\\circ$C (ERA5) to $-1.45^\\circ$C (POWER)\nover 45 years as nights warmed faster than days",
        xy=(2019, 12.23), xytext=(1983, 11.45),
        arrowprops=dict(arrowstyle="->", color="#111111", lw=1.5),
        fontsize=9.0, bbox=dict(boxstyle="round,pad=0.35", fc="#F7F7F7", ec="#555555", lw=1.0)
    )
    ax.set_xlim(1980, 2026)
    ax.set_ylim(11.1, 19.2)
    ax.set_xlabel("Year (1981–2025)")
    ax.set_ylabel("Annual Diurnal Temp. Range $\\mathrm{DTR} = T_{\\max} - T_{\\min}$ (°C)", labelpad=8)
    ax.set_title("Rapid Compression of Diurnal Temperature Range (DTR) over Sukkur (1981–2025)", pad=14)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.96, fontsize=8.8)
    fig.savefig(OUT_DIR / "Fig5_Diurnal_Temperature_Range_Compression_1981_2025.png")
    plt.close(fig)

    # FIGURE 6 Monthly Warming Asymmetry Profile
    fig, ax = plt.subplots(figsize=(10.2, 6.5))
    sub_cons_max = df_m_trd[(df_m_trd["Dataset"] == "Consensus Ensemble") & (df_m_trd["Variable"] == "Tmax")].sort_values("Month_Num")
    sub_cons_min = df_m_trd[(df_m_trd["Dataset"] == "Consensus Ensemble") & (df_m_trd["Variable"] == "Tmin")].sort_values("Month_Num")
    sub_pow_min = df_m_trd[(df_m_trd["Dataset"] == "NASA POWER (~50 km)") & (df_m_trd["Variable"] == "Tmin")].sort_values("Month_Num")
    sub_pow_max = df_m_trd[(df_m_trd["Dataset"] == "NASA POWER (~50 km)") & (df_m_trd["Variable"] == "Tmax")].sort_values("Month_Num")
    x = np.arange(1, 13)
    w = 0.20
    b1 = ax.bar(x - 1.5*w, sub_pow_min["Sen_Slope_C_per_Decade"], width=w, color="#D55E00", alpha=0.9, label="NASA POWER $T_{\\min}$ Slope (°C/decade)")
    b2 = ax.bar(x - 0.5*w, sub_cons_min["Sen_Slope_C_per_Decade"], width=w, color="#FDBF6F", edgecolor="#B22222", lw=0.8, label="Consensus Ensemble $T_{\\min}$ Slope (°C/decade)")
    b3 = ax.bar(x + 0.5*w, sub_pow_max["Sen_Slope_C_per_Decade"], width=w, color="#6BAED6", edgecolor="#08519C", lw=0.8, label="NASA POWER $T_{\\max}$ Slope (°C/decade)")
    b4 = ax.bar(x + 1.5*w, sub_cons_max["Sen_Slope_C_per_Decade"], width=w, color="#0072B2", alpha=0.9, label="Consensus Ensemble $T_{\\max}$ Slope (°C/decade)")
    ax.axhline(0, color="black", lw=1.2)
    for bars, df_sub in [(b1, sub_pow_min), (b2, sub_cons_min), (b3, sub_pow_max), (b4, sub_cons_max)]:
        for bar, (_, r_row) in zip(bars, df_sub.iterrows()):
            pval = r_row["MK_p_value"]
            val = r_row["Sen_Slope_C_per_Decade"]
            if pval < 0.05:
                star = "***" if pval < 0.001 else ("**" if pval < 0.01 else "*")
                y_pos = val + 0.018 if val >= 0 else val - 0.042
                ax.text(bar.get_x() + bar.get_width() / 2.0, y_pos, star, ha="center", va="bottom" if val >= 0 else "top", fontsize=8.8, fontweight="bold", color="#111111")
    ax.annotate(
        "Late-Winter / Spring Emergence (Feb–Apr):\nStrongest positive warming across both\n$T_{\\max}$ (+0.38 to +0.46°C/dec) & $T_{\\min}$ (+0.33 to +0.40°C/dec)",
        xy=(2.1, 0.455), xytext=(4.3, 0.47),
        arrowprops=dict(arrowstyle="->", color="#B22222", lw=1.4),
        fontsize=9.0, bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F0", ec="#D55E00", lw=1.0)
    )
    ax.annotate(
        "June Pre-Monsoon / Kharif Irrigation Split:\nDaytime $T_{\\max}$ drops (-0.54°C/dec, ***) while\nNighttime $T_{\\min}$ rises (+0.29°C/dec, *** in POWER)",
        xy=(6.3, -0.545), xytext=(6.9, -0.52),
        arrowprops=dict(arrowstyle="->", color="#0072B2", lw=1.4),
        fontsize=9.0, bbox=dict(boxstyle="round,pad=0.35", fc="#F0F8FF", ec="#0072B2", lw=1.0)
    )
    ax.set_xticks(x)
    ax.set_xticklabels(month_labels)
    ax.set_ylim(-0.76, 0.68)
    ax.set_xlabel("Month of Year (1981–2025)")
    ax.set_ylabel("Theil–Sen Trend Slope (°C / decade)", labelpad=8)
    ax.set_title("Monthly Warming Asymmetry ($T_{\\min}$ vs. $T_{\\max}$ Sen's Slopes) over Sukkur (1981–2025)", pad=14)
    ax.legend(loc="lower left", ncol=2, frameon=True, facecolor="white", framealpha=0.96, fontsize=8.8)
    fig.savefig(OUT_DIR / "Fig6_Monthly_Warming_Asymmetry_Profile.png")
    plt.close(fig)

    # FIGURE 7 (Sequential Manuscript Order: Earlier Spring Arrival of 40°C Heat DOY)
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    s_p1, _, fit_p1 = sen_line(years, df_annual["POWER_First_40C_DOY"].values)
    s_p3, _, fit_p3 = sen_line(years, df_annual["POWER_Sustained3d_40C_DOY"].values)
    s_c1, _, fit_c1 = sen_line(years, df_annual["Consensus_First_40C_DOY"].values)
    ax.plot(years, df_annual["POWER_First_40C_DOY"], "o-", color=COL_POWER, lw=1.6, ms=5.0, alpha=0.85, label="NASA POWER First $T_{\\max} \\geq 40^\\circ$C Day (1981–1995 Mean: Apr 5 → 2011–2025: Mar 27)")
    ax.plot(years, fit_p1, "--", color=COL_POWER, lw=2.5, label=f"NASA POWER First $\\geq 40^\\circ$C Trend: {s_p1*10:+.2f} days/decade ($Z = -2.16, p = 0.031$, *)")
    ax.plot(years, df_annual["POWER_Sustained3d_40C_DOY"], "^-", color="#8C2D04", lw=1.4, ms=4.5, alpha=0.75, label="NASA POWER Sustained 3-Day $\\geq 40^\\circ$C Onset (Apr 11 → Apr 2)")
    ax.plot(years, fit_p3, ":", color="#8C2D04", lw=2.2, label=f"NASA POWER Sustained 3-Day Trend: {s_p3*10:+.2f} days/decade ($Z = -2.03, p = 0.043$, *)")
    ax.plot(years, df_annual["Consensus_First_40C_DOY"], "s-", color=COL_CONS, lw=1.8, ms=4.5, alpha=0.85, label=f"Consensus Ensemble First $\\geq 40^\\circ$C Day ({s_c1*10:+.2f} days/decade, $Z = -1.87, p = 0.061$, +)")
    doy_ticks = [60, 70, 80, 91, 100, 110, 121, 130, 140]
    doy_labels = ["DOY 60\n(Mar 1)", "DOY 70\n(Mar 11)", "DOY 80\n(Mar 21)", "DOY 91\n(Apr 1)", "DOY 100\n(Apr 10)", "DOY 110\n(Apr 20)", "DOY 121\n(May 1)", "DOY 130\n(May 10)", "DOY 140\n(May 20)"]
    ax.set_yticks(doy_ticks)
    ax.set_yticklabels(doy_labels)
    ax.set_xlim(1980, 2026)
    ax.set_ylim(58, 146)
    ax.set_xlabel("Year (1981–2025)")
    ax.set_ylabel("Spring Arrival Date (Day of Year & Calendar Date)", labelpad=8)
    ax.set_title("Earlier Spring Arrival of 40°C Heat over Sukkur (1981–2025)", pad=14)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.96, fontsize=8.6)
    fig.savefig(OUT_DIR / "Fig7_Spring_Heat_Onset_DOY_1981_2025.png")
    plt.close(fig)

    # FIGURE 8 (Sequential Manuscript Order: Escalation of Tropical Warm Nights Tmin >= 30°C)
    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    s_pw, _, fit_pw = sen_line(years, df_annual["POWER_WarmNights_30C"].values)
    s_cw, _, fit_cw = sen_line(years, df_annual["Consensus_WarmNights_30C"].values)
    s_ew, _, fit_ew = sen_line(years, df_annual["ERA5_WarmNights_30C"].values)
    ax.plot(years, df_annual["POWER_WarmNights_30C"], "o-", color=COL_POWER, lw=1.6, ms=5.0, alpha=0.85, label="NASA POWER Warm Nights $T_{\\min} \\geq 30^\\circ$C (1981–1995: 42.0 → 2011–2025: 59.3 nights/yr)")
    ax.plot(years, fit_pw, "--", color=COL_POWER, lw=2.5, label=f"NASA POWER Sen's Slope: +{s_pw*10:.2f} nights/decade ($Z = +3.46, p = 0.0005$, ***)")
    ax.plot(years, df_annual["Consensus_WarmNights_30C"], "-", color=COL_CONS, lw=2.3, label=f"Consensus Ensemble Warm Nights (30.7 → 42.2 nights/yr | +{s_cw*10:.2f}/decade, $p = 0.0028$, **)")
    ax.plot(years, df_annual["ERA5_WarmNights_30C"], "s-", color=COL_ERA5, lw=1.5, ms=4.5, alpha=0.8, label=f"ERA5-Land Warm Nights $T_{{\\min}} \\geq 30^\\circ$C (Mean: 27.0 nights/yr | +{s_ew*10:.2f}/decade, ns)")
    ax.annotate(
        "41% Surge in Nocturnal Heat Stress:\nNASA POWER warm nights ($\\geq 30^\\circ$C) rose by\n+17.3 nights/yr between 1981–1995 and 2011–2025",
        xy=(2017, 84), xytext=(1998, 4.0),
        arrowprops=dict(arrowstyle="->", color=COL_POWER, lw=1.5, connectionstyle="arc3,rad=0.18"),
        fontsize=9.0, bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5F0", ec="#D55E00", lw=1.0)
    )
    ax.set_xlim(1980, 2026)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Year (1981–2025)")
    ax.set_ylabel("Warm Nights ($T_{\\min} \\geq 30^\\circ$C, nights/year)", labelpad=8)
    ax.set_title("Escalating Frequency of Sultry Tropical Nights ($T_{\\min} \\geq 30^\\circ$C) in Sukkur (1981–2025)", pad=14)
    ax.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.96, fontsize=8.7)
    fig.savefig(OUT_DIR / "Fig8_Tropical_Warm_Nights_Surge_1981_2025.png")
    plt.close(fig)
    print("All 8 figures generated in sequential order (Fig1 to Fig8).")

if __name__ == "__main__":
    main()
