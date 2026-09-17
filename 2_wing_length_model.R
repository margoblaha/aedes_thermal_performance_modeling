################################################################################
# wing_length_model.R
#
# Wing length (a proxy for adult body size) as a function of temperature,
# treatment, sex, and city.
#
# Layout: data loading and exploration -> candidate models and selection ->
# predictions and figure, mirroring the structure used in trait_models.R.
#
# Requires: _common.R; wing_data_all.RDS
# Produces: p_wing_pred (used by figure_2.R / longevity.R)
################################################################################

source("~/Dati/Experiments/0_common.R")

suppressPackageStartupMessages({
  library(tidyr)
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
  library(emmeans)
  library(broom.mixed)
  library(stringr)
})

# ==============================================================================
# 1. DATA LOADING AND EXPLORATION
# ==============================================================================

wing_data <- readRDS("wing_data_all.RDS") %>%
  filter(avg_temp != 25) %>%
  mutate(
    replic     = factor(str_extract(replic, "\\d+")),
    replic_uid = interaction(city, method_temp, avg_temp, replic, drop = TRUE)
  )

cat("\n=== Wing-length data ===\n")
cat("n rows:", nrow(wing_data), "\n")
cat("avg_temp unique:", sort(unique(wing_data$avg_temp)), "\n")
cat("method_temp unique:", levels(wing_data$method_temp), "\n")
cat("city unique:", levels(wing_data$city), "\n")
cat("wing_mm range:", round(range(wing_data$wing_mm), 4), "\n")

p_wing <- ggplot(wing_data, aes(x = avg_temp, y = wing_mm, fill = method_temp, colour = method_temp)) +
  geom_boxplot(aes(group = interaction(avg_temp, method_temp)),
               position = position_dodge2(width = 0.8, padding = 0.1),
               alpha = 0.3, outlier.shape = NA, width = 2) +
  geom_jitter(aes(colour = method_temp),
              position = position_jitterdodge(dodge.width = 2, jitter.width = 0.25),
              size = 1, alpha = 0.6, stroke = 0.5) +
  scale_fill_manual(values = METHOD_PALETTE) +
  scale_colour_manual(values = METHOD_PALETTE) +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  facet_grid(sex ~ city, scales = "fixed") +
  labs(x = "Average temperature (°C)", y = "Wing length (mm)", fill = "Treatment", colour = "Treatment") +
  base_theme() +
  theme(legend.position = "bottom")

print(p_wing)
ggsave("figures/wing_length_boxjitter.png", p_wing, width = 12, height = 6, dpi = 300)

# ==============================================================================
# 2. MODEL FITTING AND SELECTION
# ==============================================================================

m_wing_1 <- glmmTMB(
  wing_mm ~ city * sex + method_temp + avg_temp + I(avg_temp^2) + method_temp:avg_temp + (1 | replic_uid),
  data = wing_data, family = gaussian()
)
summary(m_wing_1)   # city not significant

m_wing_2 <- glmmTMB(
  wing_mm ~ sex + method_temp + avg_temp + I(avg_temp^2) + (1 | replic_uid),
  data = wing_data, family = gaussian()
)
summary(m_wing_2)

m_wing_3 <- glmmTMB(
  wing_mm ~ sex + method_temp + city * avg_temp + city * I(avg_temp^2) + (1 | replic_uid),
  data = wing_data, family = gaussian()
)
summary(m_wing_3)

m_wing_4 <- glmmTMB(
  wing_mm ~ sex + method_temp + city * avg_temp + (1 | replic_uid),
  data = wing_data, family = gaussian()
)
print(AIC(m_wing_1, m_wing_2, m_wing_3, m_wing_4))
print(anova(m_wing_1, m_wing_2, m_wing_3, m_wing_4))
summary(m_wing_4)

m_wing_5 <- glmmTMB(
  wing_mm ~ sex + method_temp + city + avg_temp + (1 | replic_uid),
  data = wing_data, family = gaussian()
)
print(AIC(m_wing_3, m_wing_5))
print(anova(m_wing_3, m_wing_5))
summary(m_wing_5)

m_wing_6 <- glmmTMB(
  wing_mm ~ sex + method_temp + avg_temp + (1 | replic_uid),
  data = wing_data, family = gaussian()   # city dropped entirely
)
print(AIC(m_wing_3, m_wing_4, m_wing_5, m_wing_6))
print(anova(m_wing_4, m_wing_6))
summary(m_wing_6)

m_wing_7 <- glmmTMB(
  wing_mm ~ sex * avg_temp + sex * method_temp + (1 | replic_uid),
  data = wing_data, family = gaussian()   # FINAL model used for predictions below
)


print(AIC(m_wing_3, m_wing_4, m_wing_5, m_wing_6, m_wing_7))
print(anova(m_wing_3, m_wing_4, m_wing_5, m_wing_6, m_wing_7))

res_wing <- simulateResiduals(m_wing_7, n = 1000)
plot(res_wing, main = "Wing length — Gaussian residuals (m_wing_7)")
testDispersion(res_wing)
testUniformity(res_wing)

res <- simulateResiduals(m_wing_7)
plotResiduals(res, wing_data$avg_temp)
plotResiduals(res, wing_data$city)
plotResiduals(res, wing_data$sex)

# ==============================================================================
# 3. PREDICTIONS AND FIGURE
#
# Uses m_wing_7 (final model: sex, avg_temp, method_temp; no city). A manual
# prediction grid is used instead of emmeans, which does not reliably average
# over the random-effects grouping factor for this model.
# ==============================================================================

pred_grid <- expand.grid(
  sex         = levels(wing_data$sex),
  method_temp = levels(wing_data$method_temp),
  avg_temp    = seq(min(wing_data$avg_temp), max(wing_data$avg_temp), length.out = 100),
  replic_uid  = NA   # NA so predict() returns population-level (fixed-effect only) predictions
)

pred_se <- predict(m_wing_7, newdata = pred_grid, se.fit = TRUE, re.form = NA, type = "response")
pred_grid$predicted <- pred_se$fit
pred_grid$lwr       <- pred_se$fit - 1.96 * pred_se$se.fit
pred_grid$upr       <- pred_se$fit + 1.96 * pred_se$se.fit

# NOTE: city is not in m_wing_7, so predicted lines would be identical across
# city panels; the figure below facets by sex only (matches the model).
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
  scale_fill_manual(values = METHOD_PALETTE) +
  scale_colour_manual(values = METHOD_PALETTE) +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  facet_wrap(~ sex, scales = "fixed") +
  labs(x = "Average temperature (°C)", y = "Wing length (mm)",
       fill = "Treatment", colour = "Treatment", tag = "C") +
  base_theme() +
  theme(legend.position = "bottom")

print(p_wing_pred)
ggsave("figures/wing_length_model_pred.png", p_wing_pred, width = 9, height = 6, dpi = 300)
