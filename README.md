# Sukkur Climate Reanalysis & Extreme Heat Trend Analysis (1981–2025)

A multi-reanalysis climate assessment of 2-meter air temperature regimes, seasonal shifts, and extreme heatwave dynamics over **Sukkur, Sindh, Pakistan (1981–2025)**, harmonizing high-resolution **ECMWF ERA5-Land** (~9 km) and **NASA POWER / MERRA-2** (~50 km) datasets.

## Project Roadmap & Workflow Phases

- [x] **Phase 1: Data Cleaning, Temporal Alignment & Master Dataset Creation** (`src/01_harmonize_and_clean.py`)
- [ ] **Phase 2: Dataset Cross-Validation (Proving Data Reliability)** (`src/02_cross_validation.py`)
- [ ] **Phase 3: Core Statistical Trend Analysis (Mann–Kendall & Theil–Sen Slopes)** (`src/03_trend_analysis.py`)
- [ ] **Phase 4: Publication-Ready Visualizations & Summary Tables** (`src/04_generate_figures.py`)
- [ ] **Phase 5: Literature Framing & Step-by-Step Manuscript Writing** (`manuscript/`)

## Repository Structure

- `data/raw/` — Raw NASA POWER daily CSV and ECMWF ERA5-Land hourly CSV datasets (1981–2025)
- `data/processed/` — Cleaned & harmonized daily (`Sukkur_Climate_Master_Daily_1981_2025.csv`, 16,436 days) and annual (`Sukkur_Climate_Master_Annual_Summary_1981_2025.csv`, 45 years) master datasets
- `src/` — Python scripts for data harmonization, cross-validation, trend analysis, and visualization
- `results/tables/` — Statistical cross-validation and trend summary tables
- `results/figures/` — High-resolution publication figures
- `manuscript/` — Paper draft and LaTeX source files
