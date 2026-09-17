################################################################################
# population_growth_model_adapted.R
#
# net reproductive rate R0(T)) built by combining the fitted trait models
# (development, survival, longevity, wing length) with literature-derived
# fecundity and egg-hatching functions. Constant-temperature data only.
#
# Requires: _common.R; all_combined_data.RDS
#           figures/thermal_dependence_of_r0.png,
#           figures/population_trajectory.png
################################################################################

source("~/Downloads/Scripts/0_common.R")

suppressPackageStartupMessages({
  library(tidyr)
  library(purrr)
  library(glmmTMB)
  library(minpack.lm)
  library(scales)
})

# ==============================================================================
# 0. DATA
# ==============================================================================

constant_data <- readRDS("experimental_data.RDS") %>%
  filter(method_temp == "constant") %>%
  mutate(city = as.character(city))

temperatures_obs <- sort(unique(constant_data$avg_temp))

cat("Temperatures observed:", paste(temperatures_obs, collapse = ", "), "\n")
cat("n individuals:", nrow(constant_data), "\n")

# ==============================================================================
# 1. LITERATURE — wing length to fecundity (Blackmore & Lord 2000)
# ==============================================================================

fecundity_from_wl <- function(wl) {
  10^(0.68 + 1.35 * log10(pmax(0.1, wl)))
}

# ==============================================================================
# 2. LITERATURE — egg hatching rate (Marini et al. 2020)
# ==============================================================================

hatch_lit <- tibble(
  temperature = c(10, 15, 25, 30),
  hatch_rate  = c(0.0, 0.20, 0.85, 0.60)
)
hatch_mod    <- lm(hatch_rate ~ temperature + I(temperature^2), data = hatch_lit)
hatch_rate_T <- function(temp) {
  pmax(0, pmin(1, predict(hatch_mod, newdata = data.frame(temperature = temp))))
}

gonotrophic_cycle_T <- function(temp) pmax(2, 15 - 0.3 * temp)

# ==============================================================================
# 3. DEVELOPMENT RATE — Brière-1 on pooled data (no city split)
# ==============================================================================

dev_time_days <- function(temp, a, Tmin, Tmax) {
  r <- briere1(temp, a, Tmin, Tmax)
  ifelse(r > 0, 1 / r, NA_real_)
}

# Zero-rate anchor at 13°C (Frankfurt and Palermo had zero adults there;
# pooling means we add one anchor row for the species as a whole)
ZERO_ANCHOR_TEMP <- 13

make_dev_data <- function(dat) {
  dat %>%
    filter(!is.na(days_larva_adult), days_larva_adult > 0) %>%
    group_by(avg_temp) %>%
    summarise(rate_l2a = mean(1 / days_larva_adult, na.rm = TRUE), n_initial = n(), .groups = "drop") %>%
    filter(avg_temp != ZERO_ANCHOR_TEMP) %>%
    bind_rows(tibble(avg_temp = ZERO_ANCHOR_TEMP, rate_l2a = 0, n_initial = 1L))
}

fit_briere_pooled <- function(dev_df, fallback = NULL) {
  pos_temps  <- dev_df$avg_temp[dev_df$rate_l2a > 0]
  Tmin_upper <- min(pos_temps) - 0.1

  tryCatch({
    fit <- nlsLM(
      rate_l2a ~ briere1(avg_temp, a, Tmin, Tmax),
      data    = dev_df,
      weights = n_initial,
      start   = list(a = 5e-5, Tmin = 10, Tmax = 38),
      lower   = c(a = 1e-8, Tmin = 0, Tmax = max(dev_df$avg_temp) + 0.5),
      upper   = c(a = 1e-2, Tmin = Tmin_upper, Tmax = 45),
      control = nls.lm.control(maxiter = 500)
    )
    as.list(coef(fit))
  }, error = function(e) {
    warning("Briere fit failed: ", conditionMessage(e))
    fallback %||% list(a = 5e-5, Tmin = 10, Tmax = 38)
  })
}

dev_data_pooled <- make_dev_data(constant_data)
briere_params   <- fit_briere_pooled(dev_data_pooled)

cat("\nPooled Briere parameters:\n"); print(briere_params)
cat("\nDevelopment data (with zero anchor):\n"); print(dev_data_pooled)

# Diagnostic plot (base R, quick check — not saved)
ts <- seq(8, 42, 0.2)
rs <- briere1(ts, briere_params$a, briere_params$Tmin, briere_params$Tmax)
plot(dev_data_pooled$avg_temp, dev_data_pooled$rate_l2a,
     pch = 16, col = "#6A00A8", xlab = "Temperature (°C)", ylab = "Dev rate (1/day)",
     main = "Brière-1 — pooled", xlim = range(ts),
     ylim = c(0, max(rs, dev_data_pooled$rate_l2a, na.rm = TRUE) * 1.15))
lines(ts, rs, col = "#b973e1", lwd = 2)

# ==============================================================================
# 4. LONGEVITY — Poisson GLM (pooled, no city)
# ==============================================================================

lon_data <- constant_data %>% filter(!is.na(adult_longevity), adult_longevity > 0)

m_longevity <- glm(adult_longevity ~ avg_temp, data = lon_data, family = poisson(),
                   control = glm.control(maxit = 100))
cat("\n--- Longevity model ---\n"); summary(m_longevity)

pred_longevity <- function(temp, mod = m_longevity) {
  nd <- data.frame(avg_temp = temp)
  pmax(1, predict(mod, newdata = nd, type = "response", allow.new.levels = TRUE))
}

# ==============================================================================
# 5. LARVAL SURVIVAL — binomial GLM (pooled, no city)
# ==============================================================================

surv_data <- constant_data %>%
  group_by(replic_uid, avg_temp) %>%
  summarise(n_total = n(), n_survived = sum(adult, na.rm = TRUE),
            n_failed = n_total - n_survived, .groups = "drop")

m_survival <- glm(cbind(n_survived, n_failed) ~ avg_temp + I(avg_temp^2),
                  data = surv_data, family = binomial(), control = glm.control(maxit = 100))
summary(m_survival)

pred_survival <- function(temp, mod = m_survival) {
  nd <- data.frame(avg_temp = temp)
  pmax(0, pmin(1, predict(mod, newdata = nd, type = "response", allow.new.levels = TRUE)))
}

# ==============================================================================
# 6. WING LENGTH — Gaussian LM, females only (pooled, no city)
# ==============================================================================

# wing_length is converted from um to mm by dividing by 1000; remove the
# division if the source data are already in mm.
wing_data <- constant_data %>%
  filter(!is.na(wing_length), sex == "F") %>%
  mutate(wing_mm = wing_length / 1000)

m_wing <- lm(wing_mm ~ avg_temp, data = wing_data)
cat("\n--- Wing length model ---\n"); summary(m_wing)

pred_wing <- function(temp, mod = m_wing) {
  nd <- data.frame(avg_temp = temp)
  pmax(0.5, predict(mod, newdata = nd, type = "response", allow.new.levels = TRUE))
}

# ==============================================================================
# 7. BIOLOGICAL REALISM CHECK
# ==============================================================================

check_temps <- seq(8, 42, by = 1)

bio_check <- tibble(
  temp      = check_temps,
  dev_rate  = briere1(check_temps, briere_params$a, briere_params$Tmin, briere_params$Tmax),
  dev_time  = dev_time_days(check_temps, briere_params$a, briere_params$Tmin, briere_params$Tmax),
  surv      = pred_survival(check_temps),
  longevity = pred_longevity(check_temps),
  wing_mm   = pred_wing(check_temps),
  fecundity = fecundity_from_wl(pred_wing(check_temps)),
  hatch     = hatch_rate_T(check_temps)
)

cat("\nBiological realism check (selected temperatures):\n")
print(bio_check %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        filter(temp %in% c(10, 13, 18, 23, 28, 33, 38)))

# Hard assertions — stop with an error if violated
stopifnot("Survival out of [0,1]"   = all(bio_check$surv >= 0 & bio_check$surv <= 1))
stopifnot("Longevity below 1 day"   = all(bio_check$longevity >= 1))
stopifnot("Wing length below 0.5mm" = all(bio_check$wing_mm >= 0.5))
stopifnot("Hatch rate out of [0,1]" = all(bio_check$hatch >= 0 & bio_check$hatch <= 1))
stopifnot("Negative dev rate"       = all(bio_check$dev_rate >= 0))
cat("All biological bounds passed.\n")

par(mfrow = c(2, 3))
plot(check_temps, bio_check$dev_rate,  type = "l", main = "Dev rate",       xlab = "°C", ylab = "1/day")
plot(check_temps, bio_check$surv,      type = "l", main = "L2A survival",   xlab = "°C", ylab = "Prob", ylim = c(0, 1))
plot(check_temps, bio_check$longevity, type = "l", main = "Adult longevity",xlab = "°C", ylab = "Days")
plot(check_temps, bio_check$wing_mm,   type = "l", main = "Wing length",    xlab = "°C", ylab = "mm")
plot(check_temps, bio_check$hatch,     type = "l", main = "Hatch rate",     xlab = "°C", ylab = "Prob", ylim = c(0, 1))
plot(check_temps, bio_check$fecundity, type = "l", main = "Fecundity",      xlab = "°C", ylab = "Eggs/cycle")
par(mfrow = c(1, 1))

# ==============================================================================
# 8. CORE METRIC FUNCTION
#    Computes R0 for one temperature.
#    Model objects are passed as arguments so the bootstrap can supply
#    refitted models.
# ==============================================================================

compute_metrics <- function(temp, bp, m_lon, m_wl, m_sv, max_age = 150, sex_ratio = 0.5) {
  na_out <- list(R0 = NA_real_)

  dev_t <- dev_time_days(temp, bp$a, bp$Tmin, bp$Tmax)
  gc    <- gonotrophic_cycle_T(temp)
  if (!is.finite(dev_t) || dev_t <= 0) return(na_out)
  if (!is.finite(gc)    || gc    <= 0) return(na_out)

  nd <- data.frame(avg_temp = temp, replic_uid = NA_character_)

  # Works for both glm (observed-data models) and glm/lm bootstrap refits
  lon <- tryCatch(pmax(1, predict(m_lon, newdata = nd, type = "response", allow.new.levels = TRUE)),
                  error = function(e) pmax(1, predict(m_lon, newdata = nd, type = "response")))
  wl  <- tryCatch(pmax(0.5, predict(m_wl, newdata = nd, type = "response", allow.new.levels = TRUE)),
                  error = function(e) pmax(0.5, predict(m_wl, newdata = nd, type = "response")))
  ls  <- tryCatch(pmax(0, pmin(1, predict(m_sv, newdata = nd, type = "response", allow.new.levels = TRUE))),
                  error = function(e) pmax(0, pmin(1, predict(m_sv, newdata = nd, type = "response"))))

  mu <- 1 / lon
  fe <- fecundity_from_wl(wl)
  h  <- hatch_rate_T(temp)

  n_cycles <- lon / gc
  R0       <- h * ls * sex_ratio * fe * n_cycles
  R0       <- if (is.finite(R0) && R0 >= 0) R0 else NA_real_

  list(R0 = R0)
}

# ==============================================================================
# 9. POINT ESTIMATES AT OBSERVED TEMPERATURES
# ==============================================================================

results <- tibble(temperature = temperatures_obs) %>%
  mutate(
    met  = map(temperature, ~ compute_metrics(.x, briere_params, m_longevity, m_wing, m_survival)),
    R0   = map_dbl(met, "R0")
  ) %>%
  select(-met)

cat("\n=== Point estimates at observed temperatures ===\n")
print(results)

# ==============================================================================
# 10. NON-PARAMETRIC BOOTSTRAP
# ==============================================================================

refit_survival <- function(dat) {
  sv <- dat %>%
    group_by(replic_uid, avg_temp) %>%
    summarise(n_survived = sum(adult, na.rm = TRUE), n_failed = n() - n_survived, .groups = "drop")
  tryCatch(
    suppressWarnings(glm(cbind(n_survived, n_failed) ~ avg_temp + I(avg_temp^2),
                         data = sv, family = binomial(), control = glm.control(maxit = 100))),
    error = function(e) m_survival
  )
}

refit_longevity <- function(dat) {
  ld <- dat %>% filter(!is.na(adult_longevity), adult_longevity > 0)
  tryCatch(
    suppressWarnings(glm(adult_longevity ~ avg_temp, data = ld, family = poisson(),
                         control = glm.control(maxit = 100))),
    error = function(e) m_longevity
  )
}

refit_wing <- function(dat) {
  wd <- dat %>% filter(!is.na(wing_length), sex == "F") %>% mutate(wing_mm = wing_length / 1000)
  tryCatch(lm(wing_mm ~ avg_temp, data = wd), error = function(e) m_wing)
}

n_boot    <- 500
temp_fine <- seq(13, 38, by = 0.1)
set.seed(42)
boot_all <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  if (b %% 50 == 0) cat("  Replicate", b, "/", n_boot, "\n")

  boot_dat <- constant_data %>% group_by(avg_temp) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()

  dev_b <- make_dev_data(boot_dat)
  bp_b  <- fit_briere_pooled(dev_b, fallback = briere_params)
  lon_b <- refit_longevity(boot_dat)
  wl_b  <- refit_wing(boot_dat)
  sv_b  <- refit_survival(boot_dat)

  boot_all[[b]] <- tibble(temperature = temp_fine) %>%
    mutate(
      met  = map(temperature, ~ compute_metrics(.x, bp_b, lon_b, wl_b, sv_b)),
      R0   = map_dbl(met, "R0"),
      boot = b
    ) %>%
    select(-met)
}

boot_df <- bind_rows(boot_all)

boot_ci <- boot_df %>%
  group_by(temperature) %>%
  summarise(
    R0_lo = quantile(R0,   0.025, na.rm = TRUE), R0_hi = quantile(R0,   0.975, na.rm = TRUE),
    .groups = "drop"
  )

curves <- tibble(temperature = temp_fine) %>%
  mutate(
    met  = map(temperature, ~ compute_metrics(.x, briere_params, m_longevity, m_wing, m_survival)),
    R0   = map_dbl(met, "R0")
  ) %>%
  select(-met) %>%
  left_join(boot_ci, by = "temperature")

# ==============================================================================
# 11. PLOTS
# ==============================================================================
T_opt_R0 <- curves$temperature[which.max(curves$R0)]
R0_opt   <- max(curves$R0,   na.rm = TRUE)

upper_R0_above_1 <- max(curves$temperature[curves$R0 > 1], na.rm = TRUE)
lower_R0_above_1 <- min(curves$temperature[curves$R0 > 1], na.rm = TRUE)
cat("R0 > 1 window: [", lower_R0_above_1, ",", upper_R0_above_1, "]°C\n")
cat("T_opt R0(T) =", round(T_opt_R0, 1), "°C  R0_opt =", round(R0_opt, 2), "\n")

# R0(T): net reproductive rate, with extrapolation zones shaded beyond the
# observed 13-33°C range
r0_plot <- ggplot(curves, aes(x = temperature)) +
  # annotate("rect", xmin = -Inf, xmax = 13, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5) +
  # annotate("rect", xmin = 33, xmax = Inf,  ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5) +
  geom_ribbon(aes(ymin = pmax(0, R0_lo), ymax = R0_hi), fill = "#6A00A8", alpha = 0.20) +
  geom_line(aes(y = R0), colour = "#6A00A8", linewidth = 1) +
  geom_point(data = results, aes(x = temperature, y = R0),
             shape = 21, size = 2.5, stroke = 0.6, fill = "white", colour = "#6A00A8") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  # annotate("text", x = 10.5, y = max(curves$R0, na.rm = TRUE) * 0.95,
  #          label = "extrapolated", size = 3, colour = "grey50", angle = 90) +
  # annotate("text", x = 35.5, y = max(curves$R0, na.rm = TRUE) * 0.95,
  #          label = "extrapolated", size = 3, colour = "grey50", angle = 90) +
  labs(x = "Temperature (°C)", y = expression(R[0](T)), tag = "E") +
  theme_classic(base_size = 13) + base_theme()

print(r0_plot)
ggsave("figures/thermal_dependence_of_r0.png", r0_plot, width = 9, height = 6, dpi = 300)

################################################################################
# References
# Blackmore & Lord (2000) J Vector Ecol 25:212-217
# Delatte et al. (2009) J Med Entomol 46:33-41
# Marini et al. (2020) Insects 11:808
################################################################################
