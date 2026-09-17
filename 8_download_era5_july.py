import cdsapi

c = cdsapi.Client()

years = [str(y) for y in range(2016, 2026)]  # edit range as needed

for yr in years:
  c.retrieve(
    "derived-era5-single-levels-daily-statistics",
    {
      "product_type": "reanalysis",
      "variable": "2m_temperature",
      "year": yr,
      "month": "07",
      "day": [
        "01","02","03","04","05","06","07","08","09","10",
        "11","12","13","14","15","16","17","18","19","20",
        "21","22","23","24","25","26","27","28","29","30","31"
      ],
      "daily_statistic": "daily_mean",
      "time_zone": "UTC+00:00",
      "frequency": "1_hourly",
      # Europe: N, W, S, E
      "area": [72, -25, 34, 45],
      "data_format": "netcdf",
      "download_format": "unarchived"
    },
    f"ERA5_dailymean_t2m_Europe_{yr}07.nc"
  )