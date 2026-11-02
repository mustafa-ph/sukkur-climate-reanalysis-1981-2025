# Shifting Summer Temperatures, Earlier Spring Heat Arrival, and Rising Warm Nights in Sukkur, Pakistan (1981–2025): Insights from ERA5-Land and NASA POWER Reanalyses

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![Study Period](https://img.shields.io/badge/Period-1981--2025%20(45%20Years)-orange.svg)](#)

This repository holds all the data, charts, and Python code for a 45-year (1981–2025) temperature study I did on Sukkur, Pakistan. I compared two major climate models (ERA5-Land and NASA POWER) to track local climate changes, and here are my key findings:

* **The oasis cooling effect:** Both datasets match up well overall, but NASA POWER records daytime highs about 3°C warmer. This happens because ERA5-Land better captures the local cooling effect from the Indus River and seasonal irrigation, while NASA treats the area more like a dry desert.
* **Shrinking temperature gaps:** Our nights are getting steadily warmer, but our core summer days have actually cooled slightly due to that farming moisture. Because of this, the gap between our daytime highs and nighttime lows is shrinking rapidly.
* **Earlier extreme heat:** Late winter and spring are warming faster than any other time of year. I found that our first 40°C day is showing up about three days earlier every decade, shifting from early April in the 1980s to late March today.
* **A surge in hot nights:** We are experiencing way more nocturnal heat stress. Sticky, tropical nights (where the temperature never drops below 30°C) have jumped by 41%, climbing from about 42 nights a year to nearly 60.

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
