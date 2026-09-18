"""
ERA5 monthly means processing for Palermo, Trento and Frankfurt.

Variable notes
--------------
2m_temperature (t2m):
    Monthly mean in Kelvin. Convert to Celsius by subtracting 273.15.
    Annual mean = mean of 12 monthly values.

total_precipitation (tp):
    Monthly mean of the hourly accumulation rate in m/hour.
    Monthly total in mm:
        monthly_total_mm = mean_rate_m_per_hr x 24 x days_in_month x 1000
    Annual total = sum of 12 monthly totals.

Output (in era5_processed/)
----------------------------
monthly_climate.csv   monthly temperature and precipitation per city
annual_climate.csv    annual mean temperature and annual total precipitation
cv_summary.csv        coefficient of variation across years per city
"""

import os
import zipfile
import calendar
import numpy as np
import pandas as pd
import xarray as xr

CITIES = {
    "Palermo":   {"lat": 38.115, "lon": 13.361},
    "Trento":    {"lat": 46.067, "lon": 11.121},
    "Frankfurt": {"lat": 50.110, "lon":  8.682},
}

INFILE  = os.path.join("era5_monthly", "era5_monthly_1995_2025.nc")
OUTDIR  = "era5_processed"
KELVIN  = 273.15
os.makedirs(OUTDIR, exist_ok=True)


# ── File handling ─────────────────────────────────────────────────────────────
def resolve_and_open(path):
    """
    The new CDS API returns a zip archive containing one NetCDF file
    per variable type (averaged vs accumulated). Extract ALL .nc files,
    open each, and merge into a single Dataset.
    Returns a merged xarray Dataset.
    """
    if zipfile.is_zipfile(path):
        print(f"  {path} is a zip archive.")
        extract_dir = os.path.dirname(path)

        with zipfile.ZipFile(path, "r") as zf:
            nc_names = [f for f in zf.namelist() if f.endswith(".nc")]
            print(f"  Contains {len(nc_names)} NetCDF file(s): {nc_names}")
            zf.extractall(extract_dir)

        nc_paths = [os.path.join(extract_dir, f) for f in nc_names]

    else:
        nc_paths = [path]

    datasets = []
    for p in nc_paths:
        print(f"  Opening {os.path.basename(p)}...")
        for engine in ("h5netcdf", "netcdf4", "scipy"):
            try:
                ds = xr.open_dataset(p, engine=engine)
                print(f"    Variables: {list(ds.data_vars)}  |  engine: {engine}")
                datasets.append(ds)
                break
            except Exception:
                continue
        else:
            raise IOError(
                f"Could not open {p} with any available engine. "
                "Try: pip install h5netcdf netCDF4 scipy"
            )

    if len(datasets) == 1:
        return datasets[0]

    # Merge all variable files into one Dataset
    # Use join="override" to handle any minor coordinate differences
    merged = xr.merge(datasets, join="override")
    print(f"\n  Merged variables: {list(merged.data_vars)}")
    return merged


# ── Grid extraction ───────────────────────────────────────────────────────────
def extract_nearest(ds, lat, lon):
    lat_dim = "latitude"  if "latitude"  in ds.dims else "lat"
    lon_dim = "longitude" if "longitude" in ds.dims else "lon"
    return ds.sel({lat_dim: lat, lon_dim: lon}, method="nearest")


def first_var(ds, candidates):
    for c in candidates:
        for v in ds.data_vars:
            if c in v.lower():
                return v
    raise KeyError(
        f"None of {candidates} found in dataset variables: "
        f"{list(ds.data_vars)}"
    )


def days_in_month(year, month):
    return calendar.monthrange(year, month)[1]


# ── Processing ────────────────────────────────────────────────────────────────
def main():
    if not os.path.exists(INFILE):
        raise FileNotFoundError(
            f"{INFILE} not found. Run 01_download_era5_monthly.py first."
        )

    print(f"Loading {INFILE}...\n")
    ds = resolve_and_open(INFILE)

    # Handle both 'time' and 'valid_time' dimension names
    time_dim = "valid_time" if "valid_time" in ds.dims else "time"
    if time_dim == "valid_time":
        ds = ds.rename({"valid_time": "time"})

    print(f"\n  Final dimensions: {dict(ds.sizes)}")
    print(f"  Final variables:  {list(ds.data_vars)}\n")

    t_var = first_var(ds, ["t2m", "temperature"])
    p_var = first_var(ds, ["tp",  "precipitation"])

    records = []

    for city_name, coords in CITIES.items():
        city_ds = extract_nearest(ds, coords["lat"], coords["lon"])

        times       = pd.to_datetime(city_ds.time.values)
        temp_k      = city_ds[t_var].values.ravel()
        precip_rate = city_ds[p_var].values.ravel()   # m/hr

        for date, t, p_rate in zip(times, temp_k, precip_rate):
            year  = date.year
            month = date.month
            ndays = days_in_month(year, month)

            temp_c    = float(t) - KELVIN
            # monthly total mm = mean hourly rate (m/hr) x 24hr x n_days x 1000
            precip_mm = max(0.0, float(p_rate) * ndays * 1000.0)

            records.append({
                "city":        city_name,
                "year":        year,
                "month":       month,
                "temp_mean_c": round(temp_c, 3),
                "precip_mm":   round(precip_mm, 2),
            })

    ds.close()

    monthly = (
        pd.DataFrame(records)
        .sort_values(["city", "year", "month"])
        .reset_index(drop=True)
    )
    monthly.to_csv(os.path.join(OUTDIR, "monthly_climate.csv"), index=False)
    print(f"Monthly data saved: {len(monthly)} rows")

    # ── Annual aggregation ────────────────────────────────────────────────────
    annual = (
        monthly
        .groupby(["city", "year"])
        .agg(
            temp_annual_mean    = ("temp_mean_c", "mean"),
            precip_annual_total = ("precip_mm",   "sum"),
        )
        .reset_index()
    )
    annual.to_csv(os.path.join(OUTDIR, "annual_climate.csv"), index=False)
    print("Annual data saved")

    # ── Coefficient of variation ──────────────────────────────────────────────
    def cv(s):
        return (s.std() / s.mean()) * 100

    cv_df = (
        annual
        .groupby("city")
        .agg(
            temp_mean_overall   = ("temp_annual_mean",    "mean"),
            temp_cv_pct         = ("temp_annual_mean",    cv),
            precip_mean_overall = ("precip_annual_total", "mean"),
            precip_cv_pct       = ("precip_annual_total", cv),
        )
        .reset_index()
    )
    cv_df.to_csv(os.path.join(OUTDIR, "cv_summary.csv"), index=False)
    print("CV summary saved\n")
    print(cv_df.to_string(index=False, float_format="{:.2f}".format))


if __name__ == "__main__":
    main()
