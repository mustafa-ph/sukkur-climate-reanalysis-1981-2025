# Raw Reanalysis Datasets (1981–2025)

Place the two raw downloaded files in this directory before running `python src/01_harmonize_and_clean.py`:
1. `POWER_Point_Daily_19810101_20251231_027d72N_068d82E_LST.csv` — NASA POWER (MERRA-2) daily 2m air temperature (`T2M`, `T2M_MAX`, `T2M_MIN`) at 27.72°N, 68.82°E (~50 km grid).
2. `reanalysis-era5-land-timeseries-sfc-2m-temperaturetxkzu847.csv` — ECMWF ERA5-Land hourly 2m air temperature (`t2m` in Kelvin, UTC) at 27.70°N, 68.80°E (~9 km grid).
