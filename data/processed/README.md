# Processed Datasets (1981–2025)

Running `src/01_harmonize_and_clean.py` generates my two main dataset files:

- `Sukkur_Climate_Master_Daily_1981_2025.csv` — All 16,436 daily records, perfectly aligned to our local Pakistan Standard Time (PKT, UTC+5).
- `Sukkur_Climate_Master_Annual_Summary_1981_2025.csv` — A year-by-year breakdown of the 45 years. I included specific metrics like extreme heat days (>=40°C and >=45°C), warm nights (>=30°C), heatwave spell durations, and the exact day of the year (DOY) that our summer heat started and retreated.
