################################################################################
# figure_2.R
#
# Assembles the main multi-panel trait figure. This script does NOT load any
# data itself — it rebuilds each panel (with combination-specific legend
# placement) from objects left in the environment by earlier scripts, then
# assembles them with patchwork.
#
# Run order required before this script:
#   1. trait_models.R                    -> pred_surv, obs_surv_mean (survival)
#                                            l2a_data, pred_l2a (L2A development time)
#   2. wing_length_model.R                -> wing_data, pred_grid
#   3. briere.R              -> curve_ci, ct_curve, ft_curve,
#                                            rs_curve, pooled_rate, palette
#   4. population_growth_model_adapted.R   -> curves, results
#
# Produces: figures/combined_figure_2_v2.png
################################################################################

source("~/Dati/Experiments/0_common.R")

suppressPackageStartupMessages({
  library(patchwork)
})

# ── Panel A: survival ────────────────────────────────────────────────────────
# No panel in this figure keeps a traditional (side/bottom) legend; panels B
# and D instead draw a small legend inside the panel itself (see their own
# theme() calls below). Panels A, C, and E have no legend at all.
pA <- ggplot() +
  geom_ribbon(data = pred_surv,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = method_temp, group = method_temp),
              alpha = 0.12) +
  geom_line(data = pred_surv,
            aes(x = avg_temp, y = predicted, colour = method_temp, group = method_temp),
            linewidth = 0.85) +
  geom_errorbar(data = obs_surv_mean,
                aes(x = avg_temp, ymin = mean - se, ymax = mean + se, colour = method_temp),
                width = 0.4, alpha = 0.8, linewidth = 0.5) +
  geom_point(data = obs_surv_mean, aes(x = avg_temp, y = mean, colour = method_temp), size = 2.0) +
  scale_colour_manual(
    values = METHOD_PALETTE, name = "Treatment",
    guide = guide_legend(override.aes = list(linetype = 1, shape = 16, size = 2, linewidth = 1))
  ) +
  scale_fill_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Average temperature (°C)", y = "Survival probability", tag = "A") +
  base_theme() +
  theme(legend.position = "none")

# ── Panel B: larva -> adult development time ─────────────────────────────────
# Draws its own small legend inside the panel (top-right).
p_l2a_pred <- ggplot() +
  geom_boxplot(data = l2a_data,
               aes(x = avg_temp, y = days_larva_adult, colour = method_temp, fill = method_temp,
                   group = interaction(avg_temp, method_temp)),
               alpha = 0.20, outlier.shape = NA, width = 0.8, position = position_dodge(1.2)) +
  geom_jitter(data = l2a_data,
              aes(x = avg_temp, y = days_larva_adult, colour = method_temp,
                  group = interaction(avg_temp, method_temp)),
              alpha = 0.25, size = 0.8,
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1.2)) +
  geom_ribbon(data = pred_l2a,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = method_temp, group = method_temp),
              alpha = 0.25) +
  geom_line(data = pred_l2a,
            aes(x = avg_temp, y = predicted, colour = method_temp, group = method_temp),
            linewidth = 1.0) +
  scale_colour_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_fill_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Average temperature (°C)", y = "Days (larva to adult)", tag = "B") +
  base_theme() +
  theme(
    legend.position         = "inside",
    legend.position.inside  = c(0.78, 0.80),   # top-right, inside the panel
    legend.background       = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
    legend.key              = element_rect(fill = "white", colour = NA),
    legend.title            = element_text(size = 11, face = "bold"),
    legend.text             = element_text(size = 11),
    legend.key.size         = unit(1, "cm")
  )

# ── Panel C: wing length ──────────────────────────────────────────────────────
p_wing_pred <- ggplot(wing_data, aes(x = avg_temp, y = wing_mm, fill = method_temp, colour = method_temp)) +
  geom_boxplot(aes(group = interaction(avg_temp, method_temp)),
               position = position_dodge2(width = 0.8, padding = 0.1),
               alpha = 0.3, outlier.shape = NA, width = 2) +
  geom_jitter(aes(colour = method_temp),
              position = position_jitterdodge(dodge.width = 2, jitter.width = 0.25),
              size = 1, alpha = 0.6, stroke = 0.5) +
  geom_ribbon(data = pred_grid,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = method_temp, group = method_temp),
              inherit.aes = FALSE, alpha = 0.25, colour = NA) +
  geom_line(data = pred_grid,
            aes(x = avg_temp, y = predicted, colour = method_temp, group = method_temp),
            inherit.aes = FALSE, linewidth = 1.0) +
  scale_fill_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_colour_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  facet_wrap(~ sex, scales = "fixed") +
  labs(x = "Average temperature (°C)", y = "Wing length (mm)", tag = "C") +
  base_theme() +
  theme(legend.position = "none")

# ── Panel D: Brière curves + Rate Summation prediction ───────────────────────
p_briere_rs <- ggplot() +
  geom_ribbon(data = curve_ci, aes(x = avg_temp, ymin = ct_lo, ymax = ct_hi, fill = "Constant"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = ct_curve, aes(x = avg_temp, y = predicted, colour = "Constant"),
            linewidth = 1, show.legend = FALSE) +
  geom_ribbon(data = curve_ci, aes(x = avg_temp, ymin = ft_lo, ymax = ft_hi, fill = "Fluctuating"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = ft_curve, aes(x = avg_temp, y = predicted, colour = "Fluctuating"),
            linewidth = 1, show.legend = FALSE) +
  geom_ribbon(data = curve_ci, aes(x = avg_temp, ymin = rs_lo, ymax = rs_hi, fill = "Rate Summation"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = rs_curve, aes(x = avg_temp, y = predicted, linetype = "Rate Summation"),
            colour = "grey40", linewidth = 1) +
  geom_errorbar(data = pooled_rate,
                aes(x = avg_temp, ymin = rate_l2a - se, ymax = rate_l2a + se,
                    colour = ifelse(method_temp == "constant", "Constant", "Fluctuating")),
                width = 0.4, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = pooled_rate,
             aes(x = avg_temp, y = rate_l2a,
                 colour = ifelse(method_temp == "constant", "Constant", "Fluctuating")),
             size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = palette, breaks = names(palette), name = NULL) +
  scale_fill_manual(values = palette, breaks = names(palette), name = NULL) +
  scale_linetype_manual(values = c("Rate Summation" = "dashed"), name = NULL) +
  scale_x_continuous(breaks = TEMP_BREAKS, limits = c(8, 40)) +
  labs(x = "Temperature (°C)", y = "Development rate (day⁻¹)", tag = "D") +
  base_theme() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.27, 0.90),   # top-right, inside the panel
    legend.background      = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
    legend.key             = element_rect(fill = "white", colour = NA),
    legend.title           = element_text(size = 11, face = "bold"),
    legend.text            = element_text(size = 11),
    legend.key.size        = unit(1, "cm")
  )

# ── Panel E: R0 curve ─────────────────────────────────────────────────────────
# No colour/fill aesthetic mapped, so no legend would appear anyway;
# legend.position = "none" is just explicit.
r0_plot <- ggplot(curves, aes(x = temperature)) +
  geom_ribbon(aes(ymin = pmax(0, R0_lo), ymax = R0_hi), fill = "#6A00A8", alpha = 0.20) +
  geom_line(aes(y = R0), linewidth = 1, colour = "#6A00A8") +
  geom_point(data = results, aes(x = temperature, y = R0),
             shape = 21, size = 2.5, stroke = 0.6, fill = "white", colour = "#6A00A8") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Temperature (°C)", y = expression(italic(R)[0](T)), tag = "E") +
  theme_classic(base_size = 13) + base_theme() +
  theme(legend.position = "none")

# ── Assembly ──────────────────────────────────────────────────────────────────
main_grid <- (pA + p_l2a_pred) / (p_wing_pred + p_briere_rs) / r0_plot

print(main_grid)
ggsave("figures/combined_figure_2_v2.png", main_grid, width = 9, height = 12, dpi = 300)

# ==============================================================================
# OPTIONAL — an earlier panel D used adult longevity (faceted by sex) instead
# of the Brière/Rate-Summation curves above. Kept for reference in case the
# figure layout is revisited.
# ==============================================================================

# p_adult_longevity_pred <- ggplot() +
#   scale_fill_manual(values = METHOD_PALETTE, limits = names(METHOD_PALETTE), drop = FALSE, name = "Treatment") +
#   scale_colour_manual(values = METHOD_PALETTE, limits = names(METHOD_PALETTE), drop = FALSE, name = "Treatment") +
#   geom_boxplot(data = obs, aes(x = avg_temp, y = prop, group = avg_temp, colour = "constant")) +
#   geom_jitter(data = obs, aes(x = avg_temp, y = prop),
#               colour = METHOD_PALETTE[["constant"]], width = 0.4, height = 0, alpha = 0.35, size = 1.2) +
#   geom_ribbon(data = pred_adult_longevity, aes(x = avg_temp, ymin = lwr, ymax = upr, fill = "constant"), alpha = 0.15) +
#   geom_line(data = pred_adult_longevity, aes(x = avg_temp, y = predicted, colour = "constant"), linewidth = 1) +
#   facet_wrap(~ sex) +
#   scale_x_continuous(breaks = TEMP_BREAKS) +
#   labs(x = "Temperature (°C)", y = "Adult longevity", tag = "D") +
#   base_theme() +
#   theme(legend.position = "none")
