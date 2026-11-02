# Shifting Summer Temperatures, Earlier Spring Heat Arrival, and Rising Warm Nights in Sukkur, Pakistan (1981–2025): Insights from ERA5-Land and NASA POWER Reanalyses

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![Study Period](https://img.shields.io/badge/Period-1981--2025%20(45%20Years)-orange.svg)](#)

This repository contains the harmonized datasets, statistical validation and trend-detection tables, high-resolution (300 DPI) figures, and reproducible Python pipeline for the 45-year (1981–2025, N = 16,436 days) climatological investigation of 2-meter air temperatures over **Sukkur, Sindh, Pakistan** (27.70°N, 68.80°E), comparing **ECMWF ERA5-Land** (0.1° × 0.1°, ~9 km grid) and **NASA POWER / MERRA-2** (0.5° × 0.625°, ~50 km grid).

---

## Key Empirical Findings

1. **Multi-Timescale Phase Coherence & The Indus Oasis Evaporative Offset:**
   - **High Phase Agreement:** Hourly-aggregated ERA5-Land (UTC+5) and daily NASA POWER exhibit strong temporal synchronicity across daily (r = 0.975–0.988, cRMSE = 1.38–1.84°C) and monthly (r = 0.992–0.995, Willmott's d = 0.951–0.996) scales.
   - **Daytime Oasis Cooling vs. Desert Sensible Heating:** While nocturnal Tmin shows near-zero systematic bias (MBE = +0.233°C), daytime Tmax in NASA POWER is systematically warmer by +3.069°C—lowest during the January Sukkur Barrage canal closure (+2.09°C) and peaking during July–September Kharif irrigation and monsoon river discharge (+3.48°C to +4.09°C).

2. **Diurnal Warming Asymmetry & DTR Compression:**
   - **Nocturnal Warming (Tmin):** Annual mean Tmin rose significantly at +0.214°C/decade (Z = +3.37, p = 0.00076***) in NASA POWER and +0.139°C/decade (Z = +2.68, p = 0.00735**) in the Consensus Ensemble.
   - **Summer Daytime Buffering (Tmax):** Core summer (May–August) Tmax moderated by −0.226 to −0.274°C/decade (p < 0.006) under intensified irrigation and monsoon moisture.
   - **DTR Compression:** Annual Diurnal Temperature Range (DTR = Tmax − Tmin) compressed rapidly across all series (−0.255°C/decade in ERA5-Land, −0.330°C/decade in NASA POWER, −0.293°C/decade in Consensus; all p ≤ 0.00003***).

3. **Earlier Spring Arrival of Extreme Heat (DOY):**
   - **February–April Warming Hotspot:** Late winter and spring represent the fastest-warming window of the calendar year (+0.33 to +0.46°C/decade, p < 0.05).
   - **Advancing 40°C Onset:** In NASA POWER, the first annual occurrence of Tmax ≥ 40°C advanced by −3.00 days/decade (Z = −2.16, p = 0.0311*, shifting from April 5 in 1981–1995 to March 27 in 2011–2025), and the first sustained 3-day ≥ 40°C spell advanced by −2.92 days/decade (Z = −2.03, p = 0.0427*, from April 11 to April 2).

4. **41% Surge in Sultry Tropical Warm Nights (Tmin ≥ 30°C):**
   - **Nocturnal Heat Stress Escalation:** Annual Tropical Warm Nights (Tmin ≥ 30°C) surged by +41.1% in NASA POWER from 42.0 nights/year (1981–1995) to 59.3 nights/year (2011–2025) at +6.20 nights/decade (Z = +3.46, p = 0.00053***), and in the Consensus Ensemble from 30.7 to 42.2 nights/year at +4.36 nights/decade (Z = +2.99, p = 0.00282**).

---

## Repository Structure

```text
sukkur-climate-reanalysis-1981-2025/
├── README.md
├── CITATION.cff
├── LICENSE
├── requirements.txt
├── scripts/
│   ├── 01_harmonize_and_clean.py       # Phase 1: UTC+5 alignment, QC & master dataset creation
│   ├── 02_cross_validation.py          # Phase 2: Multi-timescale cross-validation & oasis offset
│   ├── 03_trend_analysis.py            # Phase 3: Mann-Kendall & Theil-Sen trend detection
│   └── 04_generate_figures.py          # Phase 4: 300-DPI publication figure generation (Fig 1-8)
├── data/
│   ├── processed/
│   │   ├── Sukkur_Climate_Master_Daily_1981_2025.csv           # 16,436 daily records (1981-2025)
│   │   └── Sukkur_Climate_Master_Annual_Summary_1981_2025.csv  # 45 annual summary records
│   └── results/
│       ├── Sukkur_Phase2_CrossValidation_Metrics_1981_2025.csv
│       ├── Sukkur_Phase2_Monthly_Microclimate_Offset_1981_2025.csv
│       ├── Sukkur_Phase3_Core_Trend_Analysis_1981_2025.csv
│       ├── Sukkur_Phase3_Monthly_Trends_1981_2025.csv
│       └── Sukkur_Phase4_Decadal_Regime_Summary_1981_2025.csv
└── figures/
    ├── Fig1_CrossValidation_Daily_Scatter_Density.png
    ├── Fig2_Monthly_Climatology_Microclimate_Offset.png
    ├── Fig3_Annual_Tmin_Warming_Trends_1981_2025.png
    ├── Fig4_Annual_Tmax_Evaporative_Buffering_1981_2025.png
    ├── Fig5_Diurnal_Temperature_Range_Compression_1981_2025.png
    ├── Fig6_Monthly_Warming_Asymmetry_Profile.png
    ├── Fig7_Spring_Heat_Onset_DOY_1981_2025.png
    └── Fig8_Tropical_Warm_Nights_Surge_1981_2025.png
```

---

## Reproducing the Analysis

```bash
pip install -r requirements.txt
python3 scripts/01_harmonize_and_clean.py
python3 scripts/02_cross_validation.py
python3 scripts/03_trend_analysis.py
python3 scripts/04_generate_figures.py
```
