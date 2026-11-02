# Raw Datasets (1981–2025)

The raw hourly ERA5-Land file is around 26 MB, which is just over GitHub's file size limit. Because of this, I couldn't include the raw input files in this repository. 

If you just want to see the final data, don't worry—I've already included my fully cleaned and harmonized datasets in the [`../processed/`](../processed/) folder. 

However, if you want to run my cleaning script (`python3 scripts/01_harmonize_and_clean.py`) from scratch, you will need to download the raw files yourself and drop them into this `data/raw/` directory:

1. **NASA POWER Daily Series (`POWER_Point_Daily_19810101_20251231_027d72N_068d82E_LST.csv`)**
   - **Where to get it:** [NASA POWER Data Access Viewer](https://power.larc.nasa.gov/data-access-viewer/)
   - **Settings I used:** Agroclimatology / Daily (`1981-01-01` to `2025-12-31`, Local Solar Time)
   - **Coordinates:** Latitude `27.72°N`, Longitude `68.82°E`
   - **Variables needed:** `T2M` (Mean Temp), `T2M_MAX` (Max Temp), and `T2M_MIN` (Min Temp)

2. **ECMWF ERA5-Land Hourly Series (`reanalysis-era5-land-timeseries-sfc-2m-temperaturetxkzu847.csv`)**
   - **Where to get it:** [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land)
   - **Settings I used:** Hourly data from `1981-01-01 00:00 UTC` to `2025-12-31 23:00 UTC`
   - **Coordinates:** Latitude `27.70°N`, Longitude `68.80°E`
   - **Variable needed:** `2m_temperature`
