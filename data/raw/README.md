# Raw Reanalysis Datasets (1981–2025)

Because the raw hourly ECMWF ERA5-Land time series exceeds GitHub's 25 MB web-upload limit (~26 MB uncompressed), raw input files are retrieved directly from their official institutional archives. Note that the fully harmonized, quality-controlled daily (`N = 16,436 days`) and annual (`N = 45 years`) master datasets are already provided in [`../processed/`](../processed/).

To re-run `python3 scripts/01_harmonize_and_clean.py` from scratch, place the following two raw files in this directory (`data/raw/`):

1. **NASA POWER Agroclimatology Daily Series (`POWER_Point_Daily_19810101_20251231_027d72N_068d82E_LST.csv`)**
   - **Source Portal:** [NASA POWER Data Access Viewer](https://power.larc.nasa.gov/data-access-viewer/) (NASA Langley Research Center / GMAO MERRA-2)
   - **Community / Temporal Level:** Agroclimatology / Daily (`1981-01-01` to `2025-12-31`, Local Solar Time)
   - **Point Coordinates:** Latitude `27.72°N`, Longitude `68.82°E` (~50 km grid; elevation `62.79 m`)
   - **Parameters:** `T2M` (2m Mean Temperature, °C), `T2M_MAX` (2m Maximum Temperature, °C), `T2M_MIN` (2m Minimum Temperature, °C)

2. **ECMWF ERA5-Land Hourly Series (`reanalysis-era5-land-timeseries-sfc-2m-temperaturetxkzu847.csv`)**
   - **Source Portal:** [Copernicus Climate Change Service (C3S) Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land) (`DOI: 10.24381/cds.e2161bac`)
   - **Temporal Coverage:** Hourly (`1981-01-01 00:00 UTC` to `2025-12-31 23:00 UTC`, `394,488 hours`)
   - **Grid Node Coordinates:** Latitude `27.70°N`, Longitude `68.80°E` (`0.1° × 0.1°`, ~9 km grid)
   - **Variable:** `2m_temperature` (`t2m`, reported in Kelvin)
