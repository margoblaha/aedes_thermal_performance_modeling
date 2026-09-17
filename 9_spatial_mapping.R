# Ae. albopictus L2A development rate -- CT vs FT -- Europe, July 2016-2025 (10-yr mean)

# 0. Setup ----
library(terra)
library(lubridate)
library(ggplot2)
library(tidyterra)   # geom_spatraster()
library(sf)
library(rnaturalearth)
library(scales)      # scales::squish(), for the delta-map legend
library(patchwork)   # 3-panel composite figure
library(ggtext)

outdir <- "maps"

# Climate data -- 10 years of July daily-mean temperature ----
# Expects one file per year, named as produced by the download script:
# ERA5_dailymean_t2m_Europe_{YYYY}07.nc
data_dir <- "/home/margo/Dati/Experiments/maps/"
years    <- 2016:2025   # edit to match whichever 10 years you downloaded

# 1. Briere development-rate function ----
briere <- function(temp, a, Tmin, Tmax) {
  rate <- a * temp * (temp - Tmin) * sqrt(pmax(Tmax - temp, 0))
  rate[is.na(temp)] <- NA
  rate[!is.na(temp) & (temp < Tmin | temp > Tmax)] <- 0
  rate
}

# 2. L2A dev_rate parameters (Briere fit), CT vs FT ----
# (pulled straight from params_CT$l2a$dev_rate / params_FT$l2a$dev_rate
dev_rate_CT <- c(a = 0.0000680, Tmin = 8.79, Tmax = 37.7)
dev_rate_FT <- c(a = 0.0000486, Tmin = 10.3, Tmax = 38.0)

dev_rate_l2a <- function(temp, p) briere(temp, p["a"], p["Tmin"], p["Tmax"])

# 3. Per-year July mean development rate, then average across years ----
# For each year: load daily temps, compute daily dev_rate (CT & FT), take
# that year's July mean (same additive-rate logic as before). Store each
# year's mean layer, then average the 10 yearly means. Because every year
# contributes the same number of July days (31), the mean of yearly means
# is identical to pooling all ~310 days directly -- but doing it per-year
# means each file is processed and dropped in turn, which keeps memory
# use flat instead of holding 10 years of daily rasters at once.
yearly_CT <- vector("list", length(years))
yearly_FT <- vector("list", length(years))

for (i in seq_along(years)) {
  yr <- years[i]
  message("Processing ", yr, " ...")

  f <- file.path(data_dir, sprintf("ERA5_dailymean_t2m_Europe_%d07.nc", yr))
  daily_t <- rast(f)
  daily_t <- app(daily_t, function(x) x - 273.15)  # Kelvin -> Celsius

  july_dates <- seq(as.Date(sprintf("%d-07-01", yr)), as.Date(sprintf("%d-07-31", yr)), by = "day")
  names(daily_t) <- as.character(july_dates)
  time(daily_t)  <- july_dates

  r_d_CT <- app(daily_t, dev_rate_l2a, p = dev_rate_CT)
  r_d_FT <- app(daily_t, dev_rate_l2a, p = dev_rate_FT)

  yearly_CT[[i]] <- app(r_d_CT, mean, na.rm = TRUE)  # this year's July mean
  yearly_FT[[i]] <- app(r_d_FT, mean, na.rm = TRUE)

  rm(daily_t, r_d_CT, r_d_FT); gc()
}

r_l2a_mth_CT <- app(rast(yearly_CT), mean, na.rm = TRUE)  # 10-year mean of July means
r_l2a_mth_FT <- app(rast(yearly_FT), mean, na.rm = TRUE)

names(r_l2a_mth_CT) <- "dev_rate"
names(r_l2a_mth_FT) <- "dev_rate"

# 4. Difference: CT - FT ----
# Unlike the survival/abundance ratio in the full pipeline, dev_rate is
# never raised to a power -- briere() just clamps to exactly 0 outside
# [Tmin, Tmax] -- so it doesn't suffer the near-zero blowup that forced
# the log2fc + floor-threshold treatment there. A plain difference is
# well-behaved everywhere; no floor/log transform needed.
diff_l2a <- 100 * (r_l2a_mth_CT - r_l2a_mth_FT) / r_l2a_mth_FT
names(diff_l2a) <- "diff"

# 5. Shared plotting setup ----
e <- as.vector(ext(r_l2a_mth_CT))  # xmin, xmax, ymin, ymax

land <- ne_countries(scale = "medium", returnclass = "sf")

base_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = ggtext::element_textbox_simple(
        fill = "grey90",
        box.colour = "grey10",
        linewidth = 0.5,         # border width
        halign = 0.5,
        padding = margin(2, 2, 2, 2),
      ),
      legend.position = "bottom"
    )
}

# 6. CT / FT maps, calibrated to a shared scale ----
# CT and FT are put on the SAME fill scale so the two panels are directly
# comparable by eye -- calibrated to CT's range. Cells
# where CT exceeds FT's max are squished to the top color rather than
# turned NA (oob = scales::squish): they still register as "at or above
# FT's maximum" instead of disappearing from the map.

# If you'd rather see CT's true ceiling, switch ft_range below to
# range(c(values(r_l2a_mth_CT), values(r_l2a_mth_FT)), na.rm = TRUE).
ct_range <- range(values(r_l2a_mth_CT), na.rm = TRUE)

plot_rate <- function(rast_layer, title_txt, limits) {
  ggplot() +
    geom_spatraster(data = mask(rast_layer, vect(land)), aes(fill = dev_rate)) +
    geom_sf(data = land, fill = NA, color = "grey30", linewidth = 0.2) +
    scale_fill_viridis_b(
      name = "L2A development rate \n(1/day)",
      limits = limits,
      n.breaks = 7,
      oob = scales::squish,
      na.value = "transparent"
    ) +
    coord_sf(xlim = e[c("xmin", "xmax")], ylim = e[c("ymin", "ymax")], expand = FALSE) +
    labs(title = title_txt, x = NULL, y = NULL) +
    base_theme()+
    theme(legend.key.width = unit(1.5, "cm"))
}

map_CT <- plot_rate(r_l2a_mth_CT, "Constant-temperature parameters", ct_range)
map_FT <- plot_rate(r_l2a_mth_FT, "Fluctuating-temperature parameters", ct_range)

# 7. Delta map: CT - FT ----
# Robust (percentile-based) limit + squish, same pattern as the full
# script's delta map -- here it's just guarding against a handful of
# genuine outlier pixels, not doing floor-masking work.
scale_limit <- as.numeric(quantile(abs(values(diff_l2a)), probs = 0.90, na.rm = TRUE))
delta_breaks <- seq(-125, 125, by = 25)
delta_scale <- function() {
  scale_fill_stepsn(
    name = "\u0394 development rate\n[(CT − FT)/FT]*100",
    colours = rev(RColorBrewer::brewer.pal(9, "RdBu")),
    limits = c(-scale_limit, scale_limit),
    breaks = delta_breaks,
    oob = scales::squish,
    na.value = "transparent"
  )
}

map_delta <- ggplot() +
  geom_spatraster(data = mask(diff_l2a, land), aes(fill = diff)) +
  geom_sf(data = land, fill = NA, color = "grey30", linewidth = 0.2) +
  # scale_fill_distiller(
  #   name = "\u0394 development rate\n[(CT \u2212 FT)/FT]*100",
  #   palette = "RdBu",
  #   limits = c(-scale_limit, scale_limit),
  #   oob = scales::squish,
  #   na.value = "transparent"
  # ) +
  delta_scale() +
  coord_sf(xlim = e[c("xmin", "xmax")], ylim = e[c("ymin", "ymax")], expand = FALSE) +
  labs(
    title = "% difference",
    x = NULL, y = NULL
  ) +
  base_theme()+
  theme(legend.key.width = unit(1.5, "cm"))

# 8. 3-panel composite ----
map_composite <- (map_CT | map_FT | map_delta)

map_composite

ggsave(file.path(outdir, "l2a_dev_rate_3panel_jul2016-2025.png"), map_composite,
       width = 16, height = 6, dpi = 300)
