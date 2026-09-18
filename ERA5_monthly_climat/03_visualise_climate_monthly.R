################################################################################
# Figure 1: Climate characterisation of study populations + experimental design
#
# Panels:
#   A  Map of sampling locations (filled circle markers)
#   B  ERA5 mean seasonal temperature cycle (1995-2025)
#   C  ERA5 mean seasonal precipitation cycle (1995-2025)
#   D  ERA5 interannual coefficient of variation (temperature & precipitation)
#   E  Experimental thermal regime (constant vs fluctuating, 23°C example)
#
# Input:  era5_processed/annual_climate.csv
#         era5_processed/monthly_climate.csv
# Output: figures/figure1.pdf  +  .png
################################################################################

# ── 0. Packages ───────────────────────────────────────────────────────────────
pkgs <- c("ggplot2", "patchwork", "tidyverse", "ggrepel",
          "ggspatial", "rnaturalearth", "rnaturalearthdata", "scales")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))

dir.create("figures", showWarnings = FALSE)

# ── 1. Shared aesthetics ──────────────────────────────────────────────────────
CITY_LEVELS <- c("Frankfurt", "Trento", "Palermo")
CITY_COLS   <- c(Frankfurt = "#5E4BB5", Trento = "#2E86AB", Palermo = "#E07A5F")

# Treatment colours matching existing thermal performance figures
TREAT_COLS <- c(constant = "#7B3294", fluctuating = "#ffc425")

origins <- tibble::tribble(
  ~origin,       ~lon,   ~lat,
  "Frankfurt",   8.68,  50.12,
  "Trento",     11.12,  46.07,
  "Palermo",    13.36,  38.12
) %>%
  mutate(origin = factor(origin, levels = CITY_LEVELS))

theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey90", linewidth = 0.3),
      plot.tag           = element_text(face = "bold", size = 13),
      strip.background   = element_blank(),
      strip.text         = element_text(size = 10)
    )
}

# ── 2. Load ERA5 data ─────────────────────────────────────────────────────────
annual <- read_csv("era5_processed/annual_climate.csv", show_col_types = FALSE) %>%
  mutate(origin = factor(city, levels = CITY_LEVELS))

monthly_era5 <- read_csv("era5_processed/monthly_climate.csv",
                         show_col_types = FALSE) %>%
  mutate(origin = factor(city, levels = CITY_LEVELS))

cv_df <- annual %>%
  group_by(origin) %>%
  summarise(
    temp_mean_overall   = mean(temp_annual_mean),
    temp_cv_pct         = sd(temp_annual_mean)    / mean(temp_annual_mean)    * 100,
    precip_mean_overall = mean(precip_annual_total),
    precip_cv_pct       = sd(precip_annual_total) / mean(precip_annual_total) * 100,
    .groups = "drop"
  )

# ── 3. PANEL A — Map (filled circle markers) ──────────────────────────────────
europe_sf <- rnaturalearth::ne_countries(
  continent   = c("europe", "africa"),
  scale       = "medium",
  returnclass = "sf"
)

p_map <- ggplot() +
  geom_sf(data      = europe_sf,
          fill      = "grey93",
          colour    = "grey65",
          linewidth = 0.25) +
  geom_point(data = origins,
             aes(x = lon, y = lat, colour = origin, fill = origin),
             shape  = 21,
             size   = 4,
             stroke = 0.9) +
  ggrepel::geom_label_repel(
    data        = origins,
    aes(x = lon, y = lat, label = origin, colour = origin),
    fill        = "white",
    size        = 3.2,
    fontface    = "bold",
    box.padding = 0.6,
    show.legend = FALSE
  ) +
  ggspatial::annotation_scale(
    location   = "bl",
    width_hint = 0.22,
    text_cex   = 0.65
  ) +
  ggspatial::annotation_north_arrow(
    location    = "tl",
    which_north = "true",
    height      = unit(0.9, "cm"),
    width       = unit(0.9, "cm"),
    style       = ggspatial::north_arrow_fancy_orienteering(text_size = 8)
  ) +
  coord_sf(xlim = c(-5, 22), ylim = c(33, 57), expand = FALSE) +
  scale_colour_manual(values = CITY_COLS) +
  scale_fill_manual(values   = CITY_COLS) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", tag = "A") +
  theme_pub() +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.25))


# ── 4. PANEL B — Seasonal temperature cycle (ERA5 monthly means) ─────────────
monthly_avg <- monthly_era5 %>%
  group_by(origin, month) %>%
  summarise(
    temp_mean   = mean(temp_mean_c),
    temp_sd     = sd(temp_mean_c),
    precip_mean = mean(precip_mm),
    precip_sd   = sd(precip_mm),
    .groups     = "drop"
  ) %>%
  mutate(month_label = factor(month.abb[month], levels = month.abb))

p_seasonal_temp <- ggplot(
  monthly_avg,
  aes(x = month, y = temp_mean, colour = origin, fill = origin)
) +
  geom_ribbon(
    aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd),
    alpha = 0.15, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, stroke = 0.7) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values   = CITY_COLS, name = "Population") +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) +
  labs(
    x   = "Month",
    y   = "Mean temperature (°C)\n± 1 SD across years",
    tag = "B"
  ) +
  theme_pub()


# ── 5. PANEL C — Seasonal precipitation cycle (ERA5 monthly means) ───────────
p_seasonal_precip <- ggplot(
  monthly_avg,
  aes(x = month, y = precip_mean, colour = origin, fill = origin)
) +
  geom_ribbon(
    aes(ymin = pmax(0, precip_mean - precip_sd), ymax = precip_mean + precip_sd),
    alpha = 0.15, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, stroke = 0.7) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values   = CITY_COLS, name = "Population") +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) +
  labs(
    x   = "Month",
    y   = "Mean precipitation (mm)\n± 1 SD across years",
    tag = "C"
  ) +
  theme_pub() +
  theme(legend.position  = "none",
        legend.title     = element_text(face = "bold", size = 10),
        legend.text      = element_text(size = 9))


# ── 6. PANEL D — Interannual CV ───────────────────────────────────────────────
p_cv <- cv_df %>%
  pivot_longer(
    cols      = c(temp_cv_pct, precip_cv_pct),
    names_to  = "trait",
    values_to = "cv"
  ) %>%
  mutate(
    trait = recode(trait,
                   temp_cv_pct   = "Temperature CV (%)",
                   precip_cv_pct = "Precipitation CV (%)"),
    trait = factor(trait, levels = c("Temperature CV (%)", "Precipitation CV (%)"))
  ) %>%
  ggplot(aes(x = origin, y = cv, fill = origin)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", cv)),
            vjust = -0.5, size = 3, colour = "grey30") +
  scale_fill_manual(values = CITY_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(
    x   = NULL,
    y   = "Coefficient of variation (%)",
    tag = "D"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

# ── 7. Experimental thermal regime (23°C example, overlay) ─────────
# Constant vs fluctuating (±5°C) diel cycle at a single representative
# mean temperature (23°C).
#
# Fluctuating profile (24h cycle):
#   8h  plateau at mean + 5°C   (hours 0-8)
#   4h  linear descent to mean - 5°C   (hours 8-12)
#   8h  plateau at mean - 5°C   (hours 12-20)
#   4h  linear ascent back to mean + 5°C   (hours 20-24)

example_temp <- 23
amplitude    <- 5
high_t       <- example_temp + amplitude
low_t        <- example_temp - amplitude

hour_seq <- seq(0, 24, by = 0.1)

fluct_df <- tibble(hour = hour_seq) %>%
  mutate(
    temperature = case_when(
      hour >= 0  & hour <= 8  ~ high_t,
      hour > 8   & hour <= 12 ~ high_t + (low_t - high_t) * (hour - 8) / 4,
      hour > 12  & hour <= 20 ~ low_t,
      hour > 20  & hour <= 24 ~ low_t + (high_t - low_t) * (hour - 20) / 4
    ),
    regime = "fluctuating"
  )

const_df <- tibble(
  hour        = c(0, 24),
  temperature = example_temp,
  regime      = "constant"
)

p_regime <- ggplot() +
  geom_line(data = fluct_df,
            aes(x = hour, y = temperature, colour = regime),
            linewidth = 1.1) +
  geom_line(data = const_df,
            aes(x = hour, y = temperature, colour = regime),
            linewidth = 1.1) +
  scale_colour_manual(values = TREAT_COLS, name = "Treatment") +
  scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(low_t - 2, high_t + 2, 2)) +
  labs(
    x     = "Time of day (h)",
    y     = "Temperature (°C)"
    # tag   = "E"
  ) +
  theme_pub() +
  theme(legend.position = "bottom",
        text = element_text(size = 16))

# ── 8. Assemble and save ────────────────────────────────────────────────────────────────
top_row <- p_map | p_seasonal_temp  
bot_row <- p_seasonal_precip  | p_cv

fig1 <- (top_row / bot_row) +
  plot_layout(heights = c(1.2, 1)); fig1

ggsave("figures/figure1.png",
       fig1, width = 14, height = 9.5, units = "in", dpi = 300)

ggsave("figures/figureSM1.png",
       p_regime, width = 14, height = 9.5, units = "in", dpi = 300)

message("Saved: figures/figure1.pdf  and  figures/figure1.png")