# Trait-level GL(M)Ms for three life-history traits (pooled across cities
# unless noted): larva-to-adult survival, larva-to-pupa development time, and
# larva-to-adult development time.

#   1. Data loading and raw-data exploration
#   2. MODEL FITTING AND SELECTION for all three traits (candidate models,
#      AIC/anova comparisons, overdispersion and DHARMa diagnostics)
#   3. PREDICTIONS AND PLOTS for all three traits (emmeans tables, ggplot
#      figures), built from the final model of each trait selected in part 2
#   4. Optional / alternative models explored but not carried forward

# Requires: _common.R; all_experimental_data.RDS
# Produces: pA (survival), p_l2a_pred / pred_l2a (L2A development time)

source("~/Dati/Experiments/0_common.R")

suppressPackageStartupMessages({
  library(tidyr)
  library(readODS)      # read_ods(), if wing-length data is pulled in later
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
  library(emmeans)
  library(broom.mixed)
  library(purrr)
  library(patchwork)
})

# ==============================================================================
# 1. DATA LOADING AND RAW-DATA EXPLORATION
# ==============================================================================

all_data  <- readRDS("experimental_data.RDS")
adulthood <- all_data %>% filter(!is.na(days_larva_adult))

cat("days_larva_adult range:", min(adulthood$days_larva_adult),
    "-", max(adulthood$days_larva_adult), "\n")

cat("=== Main dataset ===\n")
cat("Total individuals:", nrow(all_data), "\n")
cat("City x method_temp x avg_temp counts:\n")
print(all_data %>% count(city, method_temp, avg_temp, generation, experiment) %>% as.data.frame())

# Observed development-time distributions (raw data, not a model prediction)
p_constant <- adulthood %>%
  filter(method_temp == "constant") %>%
  ggplot(aes(x = days_larva_adult)) +
  geom_histogram(binwidth = 1, colour = METHOD_PALETTE["constant"], fill = METHOD_PALETTE["constant"]) +
  scale_x_continuous(limits = c(0, 115)) +
  facet_grid(city ~ factor(avg_temp)) +
  theme_bw() +
  labs(title = "Constant temperatures", x = "", y = "Count")

p_fluct <- adulthood %>%
  filter(method_temp == "fluctuating") %>%
  ggplot(aes(x = days_larva_adult)) +
  geom_histogram(binwidth = 1, colour = METHOD_PALETTE["fluctuating"], fill = METHOD_PALETTE["fluctuating"]) +
  scale_x_continuous(limits = c(0, 115)) +
  facet_grid(city ~ factor(avg_temp)) +
  theme_bw() +
  labs(title = "Fluctuating temperatures", x = "Larva-to-adult development time (days)", y = "Count")

print(p_constant / p_fluct)
ggsave("figures/observed_development_times.png", p_constant / p_fluct, width = 11, height = 8, dpi = 300)

# An alternative viability-TPC approach (rTPC::kamykowski_1985 fitted per
# city x method_temp via nls_multstart) was explored but not pursued in
# favour of the binomial GLMM survival model below.

# ==============================================================================
# 2. MODEL FITTING AND SELECTION
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 TRAIT 1 — Larva-to-adult SURVIVAL [binomial]
# ------------------------------------------------------------------------------

# Aggregate to replicate level: each row = one cup
rep_data <- all_data %>%
  group_by(city, method_temp, avg_temp, replic_uid) %>%
  summarise(n_total = n(), n_survived = sum(adult, na.rm = TRUE),
            n_failed = n_total - n_survived, .groups = "drop")

# Candidate models, increasingly simplified
m_surv_bin <- glmmTMB(
  cbind(n_survived, n_failed) ~ city * avg_temp + avg_temp * method_temp + (1 | replic_uid),
  data = rep_data, family = binomial
)

m_surv_bin2 <- glmmTMB(
  cbind(n_survived, n_failed) ~ city + avg_temp + method_temp + avg_temp:method_temp + (1 | replic_uid),
  data = rep_data, family = binomial
)

m_surv_bin3 <- glmmTMB(
  cbind(n_survived, n_failed) ~ poly(avg_temp, 2) + method_temp + poly(avg_temp, 2):method_temp + (1 | replic_uid),
  data = rep_data, family = binomial   # no city term
)

cat("\n--- Survival: AIC comparison across candidate models ---\n")
print(AIC(m_surv_bin, m_surv_bin2, m_surv_bin3))   # M3 (no city) preferred
print(tapply(rep_data$n_survived / rep_data$n_total,
             as.factor(rep_data$avg_temp):rep_data$method_temp, mean))

# Final model: same structure as m_surv_bin3, refit with raw (non-orthogonal)
# polynomial terms so fixed-effect coefficients are directly interpretable
# (needed for thermal_params() in the predictions section below).
m_surv_bin4 <- glmmTMB(
  cbind(n_survived, n_failed) ~ avg_temp + I(avg_temp^2) + method_temp
    + avg_temp:method_temp + I(avg_temp^2):method_temp + (1 | replic_uid),
  data = rep_data, family = binomial
)

overdisp_fun(m_surv_bin4)
summary(m_surv_bin4)   # no significant city effect

res4 <- simulateResiduals(m_surv_bin4, n = 1000)
plot(res4, main = "Survival binomial — final model (raw quadratic)")
testResiduals(res4)
testOutliers(res4)
testResiduals(m_surv_bin3)   # cross-check against the poly() parameterisation

# Separate, smaller check: constant-temperature-only survival curve, used for
# an "Aedes albopictus mechanism" cross-reference rather than for prediction.
rep_data_const <- all_data %>%
  filter(method_temp == "constant") %>%
  group_by(avg_temp, replic_uid) %>%
  summarise(n_total = n(), n_survived = sum(adult, na.rm = TRUE),
            n_failed = n_total - n_survived, .groups = "drop")

m_surv_bin_mechanism <- glmmTMB(
  cbind(n_survived, n_failed) ~ avg_temp + I(avg_temp^2) + (1 | replic_uid),
  data = rep_data_const, family = binomial
)
print(unname(coef(m_surv_bin_mechanism)))

# ------------------------------------------------------------------------------
# 2.2 TRAIT 2 — Larva-to-pupa DEVELOPMENT TIME [Poisson]
# ------------------------------------------------------------------------------

# Only individuals that reached pupation (non-NA days_larva_pupa)
# l2p_df <- all_data %>%
#   filter(!is.na(days_larva_pupa)) %>%
#   mutate(days_larva_pupa = as.integer(round(days_larva_pupa)))
# 
# cat("\n--- Development (L2P): descriptive stats ---\n")
# cat("n rows:", nrow(l2p_df), "\n")
# cat("days_larva_pupa range:", range(l2p_df$days_larva_pupa), "\n")
# cat("mean / variance:", round(mean(l2p_df$days_larva_pupa), 2), "/",
#     round(var(l2p_df$days_larva_pupa), 2),
#     " (ratio:", round(var(l2p_df$days_larva_pupa) / mean(l2p_df$days_larva_pupa), 2), ")\n")
# # ratio >> 1 => Poisson is likely overdispersed; a negative-binomial
# # alternative was tried and is kept in the optional section below.
# 
# m_l2p_poi <- glmmTMB(
#   days_larva_pupa ~ city * avg_temp + avg_temp * method_temp + (1 | replic_uid),
#   data = l2p_df, family = poisson(link = "log")
# )
# summary(m_l2p_poi)
# cat("\n--- Development (L2P): Poisson overdispersion check ---\n")
# overdisp_fun(m_l2p_poi)
# 
# # Quadratic temperature term
# m_l2p_poi_quad <- glmmTMB(
#   days_larva_pupa ~ city * poly(avg_temp, 2) + poly(avg_temp, 2) * method_temp + (1 | replic_uid),
#   data = l2p_df, family = poisson(link = "log")
# )
# print(AIC(m_l2p_poi_quad, m_l2p_poi))
# summary(m_l2p_poi_quad)
# 
# res_l2p <- simulateResiduals(m_l2p_poi_quad, n = 1000)
# plot(res_l2p, main = "Development time (L2P) — Poisson residuals")
# testDispersion(res_l2p)
# testZeroInflation(res_l2p)

# ------------------------------------------------------------------------------
# 2.3 TRAIT 3 — Larva-to-adult DEVELOPMENT TIME [Poisson / CMP / NB]
# ------------------------------------------------------------------------------

# Only individuals that reached adulthood (non-NA days_larva_adult)
l2a_data <- all_data %>%
  filter(!is.na(days_larva_adult)) %>%
  mutate(days_larva_adult = as.integer(round(days_larva_adult)))

cat("\n--- Development (L2A): descriptive stats ---\n")
cat("n rows:", nrow(l2a_data), "\n")
cat("days_larva_adult range:", range(l2a_data$days_larva_adult), "\n")
cat("mean / variance:", round(mean(l2a_data$days_larva_adult), 2), "/",
    round(var(l2a_data$days_larva_adult), 2),
    " (ratio:", round(var(l2a_data$days_larva_adult) / mean(l2a_data$days_larva_adult), 2), ")\n")
# ratio >> 1 => Poisson is likely overdispersed

m_l2a_poi <- glmmTMB(
  days_larva_adult ~ city * avg_temp + avg_temp * method_temp + (1 | replic_uid),
  data = l2a_data, family = poisson(link = "log")
)
summary(m_l2a_poi)

res_l2a <- simulateResiduals(m_l2a_poi, n = 1000)
plot(res_l2a, main = "Development time (L2A) — Poisson residuals")
testDispersion(res_l2a)   # fails
testZeroInflation(res_l2a)

# Quadratic temperature term, with an explicit city x poly(temp) check
m_l2a_poi_quad_1 <- glmmTMB(
  days_larva_adult ~ poly(avg_temp, 2) * method_temp + poly(avg_temp, 2) * city + (1 | replic_uid),
  data = l2a_data, family = poisson(link = "log")
)

# Preferred quadratic model (city dropped — no improvement over quad_1 above)
m_l2a_poi_quad <- glmmTMB(
  days_larva_adult ~ poly(avg_temp, 2) * method_temp + (1 | replic_uid),
  data = l2a_data, family = poisson(link = "log")
)
print(emmeans(m_l2a_poi_quad, ~ method_temp, at = list(avg_temp = 23), type = "response"))
summary(m_l2a_poi_quad)
anova(m_l2a_poi_quad)

overdisp_fun(m_l2a_poi_quad)
res_l2a <- simulateResiduals(m_l2a_poi_quad, n = 1000)
plot(res_l2a, main = "Development time (L2A) — Poisson + quadratic residuals")
testDispersion(res_l2a)   # still fails
testZeroInflation(res_l2a)

# Negative binomial (NB2) with the same quadratic x city structure.
# NOTE: this model shows convergence warnings and fails the DHARMa dispersion
m_l2a_nb <- glmmTMB(
  days_larva_adult ~ city * poly(avg_temp, 2) + poly(avg_temp, 2) * method_temp + (1 | replic_uid),
  data = l2a_data, family = nbinom2(link = "log")
)
summary(m_l2a_nb)

res_l2a <- simulateResiduals(m_l2a_nb, n = 1000)
plot(res_l2a, main = "Development time (L2A) — NB residuals")
testDispersion(res_l2a)
testZeroInflation(res_l2a)

# ==============================================================================
# 3. PREDICTIONS AND PLOTS
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 TRAIT 1 — Survival: predictions, figure, and thermal-limit estimates
# ------------------------------------------------------------------------------

temp_seq <- seq(8, 40, by = 0.3)

pred_surv <- emmeans(m_surv_bin4, ~ avg_temp * method_temp,
                     at = list(avg_temp = temp_seq), type = "response") %>%
  as.data.frame() %>%
  rename(predicted = prob, lwr = asymp.LCL, upr = asymp.UCL)

obs_surv <- rep_data %>% mutate(prop = n_survived / n_total)

obs_surv_mean <- obs_surv %>%
  group_by(method_temp, avg_temp) %>%
  summarise(mean = mean(prop), se = sd(prop) / sqrt(n()), .groups = "drop")

# Pooled across cities (matches the final, no-city model m_surv_bin4)
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
  scale_colour_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_fill_manual(values = METHOD_PALETTE, name = "Treatment") +
  scale_x_continuous(breaks = TEMP_BREAKS) +
  labs(x = "Average temperature (°C)", y = "Survival probability", tag = "A") +
  base_theme()

print(pA)
ggsave("figures/pred_surv_adult_no_city.png", pA, width = 9, height = 6, dpi = 300)

# --- Analytical thermal limits (Tmin/Topt/Tmax) from the survival GLMM -----
# Pooled version: m_surv_bin4 has no city term, so thermal_params() is
# evaluated per method_temp only (constant vs. fluctuating), not per city.
fe <- fixef(m_surv_bin4)$cond

thermal_params <- function(method_name, fe) {
  int <- fe["(Intercept)"]
  b1  <- fe["avg_temp"]
  b2  <- fe["I(avg_temp^2)"]
  if (method_name == "fluctuating") {
    b1 <- b1 + fe["avg_temp:method_tempfluctuating"]
    b2 <- b2 + fe["I(avg_temp^2):method_tempfluctuating"]
  }

  # Topt: vertex of the parabola in linear-predictor space
  Topt <- -b1 / (2 * b2)

  # Tmin / Tmax: where eta(T) = 0 (p = 0.5)
  disc <- b1^2 - 4 * b2 * int
  if (disc < 0) {
    Tmin <- NA; Tmax <- NA          # curve never crosses 50%
  } else {
    roots <- (-b1 + c(-1, 1) * sqrt(disc)) / (2 * b2)
    Tmin  <- min(roots)
    Tmax  <- max(roots)
  }

  tibble(method_temp = method_name, Tmin, Topt, Tmax)
}

thermal_limits_analytical <- map_dfr(METHOD_LEVELS, thermal_params, fe = fe)
print(thermal_limits_analytical)

# # ------------------------------------------------------------------------------
# # 3.2 TRAIT 2 — L2P: coefficient table and marginal means
# # ------------------------------------------------------------------------------
# 
# tab_l2p <- tidy(m_l2p_poi_quad, effects = "fixed", conf.int = TRUE) %>%
#   mutate(across(c(estimate, std.error, conf.low, conf.high), ~ round(.x, 4)),
#          p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
# print(tab_l2p)
# 
# emm_l2p <- emmeans(m_l2p_poi_quad, ~ city * avg_temp | method_temp, type = "response")
# print(summary(emm_l2p))

# ------------------------------------------------------------------------------
# 3.3 TRAIT 3 — L2A: coefficient table, marginal means, and prediction figure
#
# Final model: m_l2a_poi_quad (quadratic Poisson, pooled across city — no
# city term), used consistently below for the coefficient table, marginal
# means, and the prediction curve/figure.
# ------------------------------------------------------------------------------

tab_l2a <- tidy(m_l2a_poi_quad, effects = "fixed", conf.int = TRUE) %>%
  mutate(across(c(estimate, std.error, conf.low, conf.high), ~ round(.x, 4)),
         p.value = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
print(tab_l2a)

emm_l2a <- emmeans(m_l2a_poi_quad, ~ avg_temp | method_temp,
                   at = list(avg_temp = TEMP_BREAKS), type = "response")
print(summary(emm_l2a))

# --- Prediction curve + figure (from survival_l2p_and_predicted_plots.R) ----
temp_seq_l2a <- seq(13, 33, by = 0.5)   # observed constant/fluctuating range only

pred_l2a <- emmeans(m_l2a_poi_quad, ~ avg_temp * method_temp,
                    at = list(avg_temp = temp_seq_l2a), type = "response") %>%
  as.data.frame() %>%
  rename(predicted = rate, lwr = asymp.LCL, upr = asymp.UCL)

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
  base_theme()

print(p_l2a_pred)
ggsave("figures/pred_devtime_l2a_poi_quadratic.png", p_l2a_pred, width = 9, height = 6, dpi = 300)