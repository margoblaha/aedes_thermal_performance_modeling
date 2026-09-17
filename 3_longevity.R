################################################################################
# longevity.R
#
# Adult longevity (constant temperature only) modelled two ways:
#   1. A Poisson/negative-binomial GLMM on longevity in days (consistent with
#      the count-based GLMM style used for the other traits in trait_models.R)
#   2. A Weibull accelerated-failure-time (AFT) survival model, as a
#      cross-check from a genuinely different model class
#
# Requires: _common.R; constant_data.RDS
# Produces: p_adult_longevity_pred (pooled), p_weibull_pred,
#           pred_adult_longevity_city / obs_city (by population — consumed by
#           figure_3.R panel D)
################################################################################

source("~/Dati/Experiments/0_common.R")

suppressPackageStartupMessages({
  library(tidyr)
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
  library(emmeans)
  library(survival)
})

const_data <- readRDS("experimental_data.RDS") %>%
  filter(method_temp == "constant")

# ==============================================================================
# 1. DATA EXPLORATION
# ==============================================================================

hist(const_data$adult_longevity)
print(range(const_data$adult_longevity, na.rm = TRUE))
boxplot(adult_longevity ~ avg_temp, data = const_data)
print(table(const_data$city, const_data$avg_temp))

adult_longevity_df <- const_data %>%
  select(city, sex, avg_temp, replic_uid, adult_longevity) %>%
  drop_na()

# ==============================================================================
# 2. MODEL FITTING AND SELECTION
# ==============================================================================

m_adult_longevity_poi <- glmmTMB(
  adult_longevity ~ city * avg_temp + sex + (1 | replic_uid),
  data = const_data, family = poisson
)

m_adult_longevity_poi2 <- glmmTMB(
  adult_longevity ~ city * avg_temp * sex + (1 | replic_uid),
  data = const_data, family = poisson
)

cat("--- LRT: full city x temp x sex interaction ---\n")
print(anova(m_adult_longevity_poi, m_adult_longevity_poi2))   # three-way term not supported

overdisp_fun(m_adult_longevity_poi)
summary(m_adult_longevity_poi)

m_adult_longevity_nb <- glmmTMB(
  adult_longevity ~ city * avg_temp + sex + (1 | replic_uid),
  data = const_data, family = nbinom2(link = "log")
)
overdisp_fun(m_adult_longevity_nb)

m_adult_longevity_quad <- glmmTMB(
  adult_longevity ~ city * (avg_temp + I(avg_temp^2)) + sex + (1 | replic_uid),
  data = const_data, family = nbinom2(link = "log")
)
overdisp_fun(m_adult_longevity_quad)

cat("--- AIC: Poisson vs NB vs NB-quadratic ---\n")
print(AIC(m_adult_longevity_poi, m_adult_longevity_nb, m_adult_longevity_quad))   # Poisson lowest
print(anova(m_adult_longevity_poi, m_adult_longevity_nb, m_adult_longevity_quad))

# Final model: m_adult_longevity_poi (lowest AIC)
res_adult_longevity <- simulateResiduals(m_adult_longevity_poi, n = 1000)
plot(res_adult_longevity, main = "Adult longevity — Poisson residuals")
testDispersion(res_adult_longevity)
testZeroInflation(res_adult_longevity)
plotResiduals(res_adult_longevity)

# Cross-check: diagnostics for the (non-selected) quadratic NB model
res_adult_longevity <- simulateResiduals(m_adult_longevity_quad, n = 1000)
plot(res_adult_longevity, main = "Adult longevity — NB-quadratic residuals")
testDispersion(res_adult_longevity)
testZeroInflation(res_adult_longevity)
summary(m_adult_longevity_quad)

# ==============================================================================
# 3. PREDICTIONS AND FIGURE (pooled across city)
# ==============================================================================

emm_adult_longevity <- emmeans(m_adult_longevity_poi, ~ sex + avg_temp, type = "response")
print(summary(emm_adult_longevity))

temp_seq_all <- seq(13, 33, by = 0.5)   # constant-temperature range only

pred_adult_longevity <- emmeans(m_adult_longevity_poi, ~ avg_temp + sex,
                                at = list(avg_temp = temp_seq_all), type = "response") %>%
  as.data.frame() %>%
  rename(predicted = rate, lwr = asymp.LCL, upr = asymp.UCL)

obs <- adult_longevity_df %>%
  group_by(city, sex, avg_temp, replic_uid) %>%
  summarise(prop = mean(adult_longevity, na.rm = TRUE), .groups = "drop")

p_adult_longevity_pred <- ggplot() +
  scale_fill_manual(values = c(constant = unname(METHOD_PALETTE["constant"])), name = "Treatment") +
  scale_colour_manual(values = c(constant = unname(METHOD_PALETTE["constant"])), name = "Treatment") +
  geom_boxplot(data = obs, aes(x = avg_temp, y = prop, group = avg_temp, colour = "constant")) +
  geom_jitter(data = obs, aes(x = avg_temp, y = prop),
              colour = unname(METHOD_PALETTE["constant"]), width = 0.4, height = 0, alpha = 0.35, size = 1.2) +
  geom_ribbon(data = pred_adult_longevity,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = "constant"), alpha = 0.15) +
  geom_line(data = pred_adult_longevity,
            aes(x = avg_temp, y = predicted, colour = "constant"), linewidth = 1) +
  facet_wrap(~ sex) +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Temperature (°C)", y = "Adult longevity", tag = "D") +
  base_theme()

print(p_adult_longevity_pred)
ggsave("figures/pred_adult_longevity_no_city.png", p_adult_longevity_pred, width = 9, height = 6, dpi = 300)

# NOTE: an earlier combined-figure assembly
# (pA + p_l2a_pred) / (p_wing_pred + p_adult_longevity_pred), using objects
# from trait_models.R and wing_length_model.R, has been superseded by the
# multi-panel layout built in figure_2.R.

# ==============================================================================
# 3b. PREDICTIONS AND FIGURE (by city — for figure_3.R panel D)
#
# Same model (m_adult_longevity_poi) as section 3, but city is now included
# in the emmeans formula so predictions and observed data are broken out by
# population rather than pooled.
# ==============================================================================

pred_adult_longevity_city <- emmeans(m_adult_longevity_poi, ~ avg_temp + city + sex,
                                     at = list(avg_temp = temp_seq_all), type = "response") %>%
  as.data.frame() %>%
  rename(predicted = rate, lwr = asymp.LCL, upr = asymp.UCL) %>%
  mutate(city = factor(city, levels = CITY_LEVELS))

obs_city <- adult_longevity_df %>%
  group_by(city, sex, avg_temp, replic_uid) %>%
  summarise(prop = mean(adult_longevity, na.rm = TRUE), .groups = "drop") %>%
  mutate(city = factor(city, levels = CITY_LEVELS))

# show.legend = FALSE on the boxplot/jitter keeps the shared bottom legend in
# figure_3.R (plot_layout(guides = "collect")) to a single "Population" key
# from the ribbon/line, rather than duplicating it with boxplot/jitter glyphs.
p_adult_longevity_by_city <- ggplot() +
  scale_fill_manual(values = CITY_COLS, name = "Population", breaks = CITY_LEVELS) +
  scale_colour_manual(values = CITY_COLS, name = "Population", breaks = CITY_LEVELS) +
  geom_boxplot(data = obs_city,
               aes(x = avg_temp, y = prop, group = interaction(city, avg_temp), colour = city),
               show.legend = FALSE, position = position_dodge2(width = 3, preserve = "single"),
               width = 1.2, outlier.shape = NA) +
  geom_jitter(data = obs_city, aes(x = avg_temp, y = prop, colour = city),
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 3),
              alpha = 0.35, size = 1.2, show.legend = FALSE) +
  geom_ribbon(data = pred_adult_longevity_city,
              aes(x = avg_temp, ymin = lwr, ymax = upr, fill = city), alpha = 0.15) +
  geom_line(data = pred_adult_longevity_city, aes(x = avg_temp, y = predicted, colour = city),
            linewidth = 0.75) +
  facet_wrap(~ sex, nrow = 2) +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Temperature (°C)", y = "Adult longevity (days)", tag = "D") +
  guides(colour = guide_legend(override.aes = list(shape = NA, linetype = NA, fill = NA))) +
  base_theme() +
  theme(legend.position = "none")

print(p_adult_longevity_by_city)
ggsave("figures/pred_adult_longevity_by_city.png", p_adult_longevity_by_city, width = 9, height = 8, dpi = 300)

# ==============================================================================
# 4. Additional check — Weibull survival model, faceted by city
#
# Cross-check of the GLMM result above using a genuinely different model
# class (parametric accelerated-failure-time survival), retaining city as an
# explicit facet rather than pooling it away.
# ==============================================================================

adult_longevity_surv <- const_data %>%
  select(city, avg_temp, replic_uid, sex, adult_longevity) %>%
  drop_na() %>%
  mutate(event = 1L) %>%       # all individuals died — no censoring
  filter(adult_longevity > 0)

cat("\n--- Weibull AFT input data ---\n")
cat("n rows:", nrow(adult_longevity_surv), "\n")
cat("adult_longevity range:", range(adult_longevity_surv$adult_longevity), "\n")
print(table(adult_longevity_surv$city, adult_longevity_surv$avg_temp))

obs_weibull <- adult_longevity_surv %>%
  group_by(city, avg_temp, replic_uid) %>%
  summarise(prop = mean(adult_longevity, na.rm = TRUE), .groups = "drop")

m_weibull <- survreg(
  Surv(adult_longevity, event) ~ city * avg_temp + sex,,
  data = adult_longevity_surv, dist = "weibull"
)
summary(m_weibull)