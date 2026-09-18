"""
ERA5 monthly means download for Palermo, Trento and Frankfurt
1995-2025: 2m temperature and total precipitation

Dataset: reanalysis-era5-single-levels-monthly-means
Resolution: ERA5 ~0.25 degrees (~28 km)
Download size: tiny (one file for the full 30-year period)

Setup
-----
1. Register at https://cds.climate.copernicus.eu/
2. Create ~/.cdsapirc:
   url: https://cds.climate.copernicus.eu/api
   key: <your-api-key>

Output
------
era5_monthly/era5_monthly_1995_2025.nc
"""

import cdsapi
import os

OUTDIR = "era5_monthly"
os.makedirs(OUTDIR, exist_ok=True)

OUTFILE = os.path.join(OUTDIR, "era5_monthly_1995_2025.nc")

# Bounding box covering all three cities [North, West, South, East]
AREA = [52, 6, 36, 16]

YEARS  = [str(y) for y in range(1995, 2026)]
MONTHS = [f"{m:02d}" for m in range(1, 13)]


def main():
    if os.path.exists(OUTFILE):
        print(f"File already exists: {OUTFILE}")
        print("Delete it to re-download.")
        return

    client = cdsapi.Client()

    print("Requesting ERA5 monthly means (1995-2025)...", flush=True)

    client.retrieve(
        "reanalysis-era5-single-levels-monthly-means",
        {
            "product_type": "monthly_averaged_reanalysis",
            "variable": [
                "2m_temperature",
                "total_precipitation",
            ],
            "year":        YEARS,
            "month":       MONTHS,
            "time":        "00:00",
            "area":        AREA,
            "data_format": "netcdf",
        },
        OUTFILE,
    )

    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
