# Data and Analysis Code

Data and analysis code for the manuscript "Diurnal temperature fluctuations outweigh population origin in shaping thermal performance and range predictions in European *Aedes albopictus*", testing whether diurnal temperature fluctuations or population origin more strongly shape larval and adult performance across a 1,200 km latitudinal gradient (Frankfurt–Rhine-Main, Germany; Trento and Palermo, Italy).

F1 offspring from three established invasive populations were reared at five mean temperatures (13, 18, 23, 28 and 33°C) under constant and fluctuating (±5°C stepped diel cycle) thermal regimes. Larval survival, larva-to-adult development time, adult wing length and adult longevity were measured, and the resulting thermal performance relationships were projected across Europe using ERA5 reanalysis data.

## Repository Contents
### Scripts

Scripts are numbered in execution order. `0_common.R` must be sourced first.

| Script | Purpose |
|---|---|
| `0_common.R` | Shared setup: package loading, data import, factor level definitions, colour palettes and plotting theme. |
| `1_survival_development-time_models.R` | Binomial GLMM for larva-to-adult survival; Poisson GLMM for development time. |
| `2_wing_length_model.R` | Gaussian LMM for adult wing length with sex × temperature and sex × regime interactions. |
| `3_longevity.R` | Poisson GLMM for adult longevity (constant-temperature treatment only), with Weibull AFT sensitivity check. |
| `4_briere.R` | Brière-1 thermal performance curve fitting for constant- and fluctuating-temperature development rates; rate summation bias quantification with bootstrap confidence intervals. |
| `5_person_period_competing_risks.R` | Person-period expansion and cause-specific discrete-time hazard GLMMs (cloglog link) for larval/pupal death and adult emergence as competing outcomes. |
| `6_adult_longevity_cumulative_incidence.R` | Discrete-time hazard model and cumulative incidence curves for adult mortality; pairwise population comparisons. |
| `7_population_growth_model_adapted.R` | Integration of trait functions into a temperature-dependent net reproductive rate $R_0(T)$ model with bootstrap uncertainty. |
| `8_download_era5_july.py` | Retrieval of ERA5 single-level daily mean 2 m air temperature for July 2016–2025 via the Copernicus Climate Data Store API. |
| `9_spatial_mapping.R` | Application of constant- and fluctuating-temperature Brière parameters to ERA5 data; continental-scale development rate and bias maps. |
| `10a_figure_2.R` | Main text Figure 2: temperature-dependent life-history traits, thermal performance curves and $R_0(T)$. |
| `10b_figure_3.R` | Main text Figure 3: cause-specific cumulative incidence of larval/pupal death, adult emergence and adult mortality. |

## Requirements

### R

**R ≥ 4.5.0** with the following packages:

- `glmmTMB`
- `survival`
- `survminer`
- `minpack.lm`
- `emmeans`
- `DHARMa`
- `terra`
- `sf`
- `tidyverse`
- `patchwork`

### Python

**Python ≥ 3.9** with:

- `cdsapi`

The `cdsapi` package is required for `8_download_era5_july.py`.

Running the ERA5 download script requires:

1. A registered **Copernicus Climate Data Store (CDS)** account.
2. A local `.cdsapirc` credentials file configured with your CDS API credentials.
