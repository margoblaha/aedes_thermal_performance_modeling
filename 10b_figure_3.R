################################################################################
# figure_3.R
#
# Assembles the cumulative-incidence multi-panel figure. This script does NOT
# load any data itself — it rebuilds each panel from objects left in the
# environment by earlier scripts, then assembles them with patchwork.
#
# Run order required before this script:
#   1. person_period_competing_risks.R
#        -> p_death_long, p_eclose_long, CITY_COLS (via _common.R)
#   2. adult_longevity_cumulative_incidence.R
#        -> p_ci_long
#   3. longevity.R -> pred_adult_longevity_city, obs_city (city-stratified)
#
# Produces: figures/combined_figure_3.png
################################################################################

source("~/Dati/Experiments/_common.R")

suppressPackageStartupMessages({
  library(patchwork)
})

# ── Panel B: cumulative incidence of larval/pupal death ──────────────────────
plot_1a <- ggplot(p_death_long, aes(x = day, y = CI_death, colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_death_lo, ymax = CI_death_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(method_temp ~ avg_temp, labeller = labeller(avg_temp = function(x) paste0(x, " °C"))) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(labels = scales::label_percent(suffix = ""), limits = c(0, 1)) +
  labs(x = "Days from hatching", y = "Cumulative incidence \n  of larval and pupal death (%)", tag = "B") +
  base_theme_hazard() +
  theme(legend.position = "none")

# ── Panel A: cumulative incidence of adult emergence (eclosion) ──────────────
plot_1b <- ggplot(p_eclose_long, aes(x = day, y = CI_eclose, colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_eclose_lo, ymax = CI_eclose_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(method_temp ~ avg_temp, labeller = labeller(avg_temp = function(x) paste0(x, " °C"))) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(labels = scales::label_percent(suffix = ""), limits = c(0, 1)) +
  labs(x = "Days from hatching", y = "Cumulative incidence \n  of adult emergence (%)", tag = "A") +
  base_theme_hazard() +
  theme(legend.position = "none")

# ── Panel C: cumulative incidence of adult death across the life span ────────
plot_ci <- ggplot(p_ci_long, aes(x = day, y = CI_death, colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_death_lo, ymax = CI_death_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(sex ~ avg_temp, labeller = labeller(avg_temp = function(x) paste0(x, " °C"))) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  scale_y_continuous(labels = scales::label_percent(suffix = ""), limits = c(0, 1)) +
  labs(x = "Days from hatching", y = "Cumulative incidence \n of death across the life span (%)", tag = "C") +
  base_theme_hazard() +
  theme(legend.position = "none")

# ── Panel D: adult longevity by population ───────────────────────────────────
p_adult_longevity_pred <- ggplot() +
  scale_fill_manual(values = CITY_COLS, name = "Population", breaks = CITY_LEVELS) +
  scale_colour_manual(values = CITY_COLS, name = "Population", breaks = CITY_LEVELS) +
  # Dodged by city so the three populations don't overplot at the same temperature
  geom_boxplot(data = obs_city,
               aes(x = avg_temp, y = prop, group = interaction(city, avg_temp), colour = city),
               show.legend = FALSE, position = position_dodge2(width = 3, preserve = "single"),
               width = 1.2, outlier.shape = NA) +
  geom_jitter(data = obs_city, aes(x = avg_temp, y = prop, colour = city),
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 3),
              alpha = 0.35, size = 1.2, show.legend = FALSE) +
  geom_ribbon(data = pred_adult_longevity_city,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = city), alpha = 0.15, show.legend = FALSE) +
  geom_line(data = pred_adult_longevity_city, aes(x = avg_temp, y = predicted, colour = city),
            linewidth = 0.75, show.legend = FALSE) +
  facet_wrap(~ sex, nrow = 2) +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Temperature (°C)", y = "Adult longevity (days)", tag = "D") +
  guides(colour = guide_legend(override.aes = list(shape = NA, linetype = NA, fill = NA))) +
  base_theme_hazard() +
  theme(legend.position = "none")

# ── Assembly ──────────────────────────────────────────────────────────────────
bottom_row <- (plot_ci + p_adult_longevity_pred) +
  plot_layout(widths = c(1.7, 1))   # C gets 1.7x the width of D — tweak ratio as needed

figure_3 <- (plot_1b / plot_1a / bottom_row) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(figure_3)
ggsave("figures/combined_figure_3.png", figure_3, width = 9, height = 12, dpi = 300)
