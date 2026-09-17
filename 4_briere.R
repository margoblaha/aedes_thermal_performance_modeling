################################################################################
# Three thermal performance curves (TPCs) are compared on a single panel:
#   1. CT          — Brière-1 fitted to constant-temperature data (pooled)
#   2. FT-empirical — Brière-1 fitted to fluctuating-temperature data (pooled)
#   3. RS          — predicted FT rate from the CT TPC via Rate Summation
#                    (mean of CT-Brière over a stepwise ±5°C diel cycle)
#
# Requires: _common.R; briere_no_city.R must be sourced first (uses its
#           `rate_data` and `briere1()`)
# Produces: p_briere_rs, curve_ci, ct_curve, ft_curve, rs_curve, pooled_rate
################################################################################

source("~/Dati/Experiments/0_common.R")

suppressPackageStartupMessages({
  library(tidyr)
  library(purrr)
  library(broom)
  library(minpack.lm)
  library(patchwork)
})

# ==============================================================================
# 1. POOL ACROSS CITIES: one mean rate per (avg_temp × method_temp)
# ==============================================================================
all_data <- readRDS("experimental_data.RDS") %>%
  mutate(
    avg_temp    = as.numeric(avg_temp),
    method_temp = factor(method_temp, levels = METHOD_LEVELS),
    replic      = factor(replic),
    sex         = factor(sex, levels = c("F", "M")),
    adult       = as.integer(!is.na(adult_emergence)),
    replic_uid  = interaction(city, method_temp, avg_temp, replic, ID_individual, drop = TRUE)
  )

cat("=== Main dataset ===\n")
cat("Total individuals:", nrow(all_data), "\n")
cat("City × method_temp × avg_temp counts:\n")
print(all_data %>% count(city, method_temp, avg_temp) %>% as.data.frame())

n_initial <- all_data %>%
  count(method_temp, avg_temp, name = "n_initial")

# ==============================================================================
# 1. PREPARE RATE DATA
#
# Development RATE = 1 / days (larva-to-adult here; swap for l2p if needed).
# Only individuals that completed development (non-NA days) are included —
# deaths before pupation have no rate and are correctly excluded. Zero-survival
# conditions (e.g. Frankfurt at 13°C) simply do not appear in the data; the
# curve is fit to where development occurred, and T_min is estimated from the
# curve shape rather than from an explicit zero-survival anchor.
# ==============================================================================

rate_data <- all_data %>%
  filter(!is.na(days_larva_adult), days_larva_adult > 0) %>%
  mutate(rate_l2a = 1 / days_larva_adult) %>%
  left_join(n_initial, by = c("method_temp", "avg_temp"))

pooled_rate <- rate_data %>%
  group_by(method_temp, avg_temp) %>%
  summarise(
    mean_rate = mean(rate_l2a, na.rm = TRUE),
    sd_rate   = sd(rate_l2a,   na.rm = TRUE),
    n         = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    se       = sd_rate / sqrt(n),
    rate_l2a = mean_rate    # rename for downstream code
  ) %>%
  select(method_temp, avg_temp, rate_l2a, n, se)

cat("=== Pooled rates by method × temperature ===\n")
print(pooled_rate)

# ==============================================================================
# 2. FIT BRIÈRE-1 TO POOLED DATA — CT AND FT SEPARATELY
# ==============================================================================

fit_briere_pooled <- function(df) {
  pos_temps  <- df$avg_temp[df$rate_l2a > 0]
  Tmin_upper <- if (length(pos_temps)) min(pos_temps) - 0.1 else 15

  tryCatch(
    nlsLM(
      rate_l2a ~ briere1(avg_temp, a, Tmin, Tmax),
      data    = df,
      weights = n,
      start   = list(a = 5e-5, Tmin = 10, Tmax = 38),
      lower   = c(a = 1e-8, Tmin = 0, Tmax = max(df$avg_temp) + 0.5),
      upper   = c(a = 1e-2, Tmin = Tmin_upper, Tmax = 45),
      control = nls.lm.control(maxiter = 500)
    ),
    error = function(e) {
      message("Pooled Brière fit failed: ", conditionMessage(e))
      NULL
    }
  )
}

ct_data <- pooled_rate %>% filter(method_temp == "constant")
ft_data <- pooled_rate %>% filter(method_temp == "fluctuating")

fit_ct <- fit_briere_pooled(ct_data)
fit_ft <- fit_briere_pooled(ft_data)

ct_params <- as.list(coef(fit_ct))
ft_params <- as.list(coef(fit_ft))

cat("\nCT Brière parameters:\n"); print(ct_params)
cat("\nFT Brière parameters:\n"); print(ft_params)

ct_params$Topt <- t_opt_briere1(ct_params$Tmin, ct_params$Tmax)
ft_params$Topt <- t_opt_briere1(ft_params$Tmin, ft_params$Tmax)

# ==============================================================================
# 3. RATE SUMMATION: project the CT TPC through the FT daily cycle
#
# Stepwise diel profile: 4h neutral / 8h warm (+amp) / 4h neutral / 8h cool
# (-amp); the two 4h transition blocks sit at the mean, so the weighted
# average below collapses to (8 cold + 8 mean + 8 warm) / 24.
# ==============================================================================

FT_AMPLITUDE <- 5   # ±5°C around the mean -> 10°C total daily range

rs_predict <- function(mean_temp, params, amp = FT_AMPLITUDE) {
  r_cold <- briere1(mean_temp - amp, params$a, params$Tmin, params$Tmax)
  r_mean <- briere1(mean_temp,       params$a, params$Tmin, params$Tmax)
  r_warm <- briere1(mean_temp + amp, params$a, params$Tmin, params$Tmax)
  (8 * r_cold + 8 * r_mean + 8 * r_warm) / 24
}

temp_fine <- seq(5, 40, by = 0.2)
rs_curve  <- tibble(
  avg_temp  = temp_fine,
  predicted = sapply(temp_fine, rs_predict, params = ct_params)
)

rs_Tmin <- temp_fine[which(rs_curve$predicted > 0)[1]] %||% NA_real_
rs_Tmax <- {
  pos <- which(rs_curve$predicted > 0)
  if (length(pos)) temp_fine[max(pos)] else NA_real_
}
rs_Topt <- temp_fine[which.max(rs_curve$predicted)]

cat("\nRS-curve parameters (from CT TPC + ±5°C cycle):\n")
cat("  Tmin =", round(rs_Tmin, 2), "  Topt =", round(rs_Topt, 2),
    "  Tmax =", round(rs_Tmax, 2), "\n")

# ==============================================================================
# 4. BOOTSTRAP — 1000 replicates, resampling individuals within (method × temp)
# ==============================================================================

N_BOOT <- 1000
set.seed(42)

# Refit pooled Brière on a resampled individual-level dataset
fit_pooled_safe <- function(indiv_df) {
  pooled <- indiv_df %>%
    group_by(avg_temp) %>%
    summarise(rate_l2a = mean(rate_l2a, na.rm = TRUE), n = n(), .groups = "drop")
  fit <- fit_briere_pooled(pooled)
  if (is.null(fit)) return(NULL)
  as.list(coef(fit))
}

cat("\nRunning bootstrap (n =", N_BOOT, ")...\n")
boot_records <- vector("list", N_BOOT)

for (b in seq_len(N_BOOT)) {
  if (b %% 100 == 0) cat("  Replicate", b, "/", N_BOOT, "\n")

  rs_boot <- rate_data %>%
    group_by(method_temp, avg_temp) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    ungroup()

  ct_b <- rs_boot %>% filter(method_temp == "constant")
  ft_b <- rs_boot %>% filter(method_temp == "fluctuating")

  ct_p <- fit_pooled_safe(ct_b)
  ft_p <- fit_pooled_safe(ft_b)
  if (is.null(ct_p) || is.null(ft_p)) next

  boot_records[[b]] <- tibble(
    boot     = b,
    avg_temp = temp_fine,
    ct_rate  = briere1(temp_fine, ct_p$a, ct_p$Tmin, ct_p$Tmax),
    ft_rate  = briere1(temp_fine, ft_p$a, ft_p$Tmin, ft_p$Tmax),
    rs_rate  = sapply(temp_fine, rs_predict, params = ct_p),
    ct_Tmin  = ct_p$Tmin, ct_Tmax = ct_p$Tmax,
    ct_Topt  = t_opt_briere1(ct_p$Tmin, ct_p$Tmax),
    ft_Tmin  = ft_p$Tmin, ft_Tmax = ft_p$Tmax,
    ft_Topt  = t_opt_briere1(ft_p$Tmin, ft_p$Tmax)
  )
}

boot_df <- bind_rows(boot_records)
cat("Bootstrap complete — successful replicates:",
    length(unique(boot_df$boot)), "/", N_BOOT, "\n")

# ==============================================================================
# 5. BUILD CURVE CIs AND PARAMETER CIs
# ==============================================================================

curve_ci <- boot_df %>%
  group_by(avg_temp) %>%
  summarise(
    ct_lo = quantile(ct_rate, 0.025, na.rm = TRUE),
    ct_hi = quantile(ct_rate, 0.975, na.rm = TRUE),
    ft_lo = quantile(ft_rate, 0.025, na.rm = TRUE),
    ft_hi = quantile(ft_rate, 0.975, na.rm = TRUE),
    rs_lo = quantile(rs_rate, 0.025, na.rm = TRUE),
    rs_hi = quantile(rs_rate, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ct_curve <- tibble(avg_temp = temp_fine,
                   predicted = briere1(temp_fine, ct_params$a, ct_params$Tmin, ct_params$Tmax))
ft_curve <- tibble(avg_temp = temp_fine,
                   predicted = briere1(temp_fine, ft_params$a, ft_params$Tmin, ft_params$Tmax))

param_ci <- boot_df %>%
  distinct(boot, ct_Tmin, ct_Tmax, ct_Topt, ft_Tmin, ft_Tmax, ft_Topt) %>%
  pivot_longer(-boot, names_to = c("source", "param"), names_sep = "_", values_to = "value") %>%
  group_by(source, param) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    lo     = quantile(value, 0.025, na.rm = TRUE),
    hi     = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# RS-curve parameters: Tmin/Topt/Tmax extracted from each bootstrap replicate's
# own RS curve, giving bootstrap CIs consistent with param_ci above.
rs_param_ci <- boot_df %>%
  group_by(boot) %>%
  summarise(
    Tmin = avg_temp[which(rs_rate > 0)[1]],
    Tmax = {
      pos <- which(rs_rate > 0)
      if (length(pos)) max(avg_temp[pos]) else NA_real_
    },
    Topt = avg_temp[which.max(rs_rate)],
    .groups = "drop"
  ) %>%
  pivot_longer(-boot, names_to = "param", values_to = "value") %>%
  group_by(param) %>%
  summarise(
    source = "rs",
    median = median(value, na.rm = TRUE),
    lo     = quantile(value, 0.025, na.rm = TRUE),
    hi     = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

all_params <- bind_rows(param_ci, rs_param_ci) %>%
  mutate(
    source = factor(source, levels = c("ct", "ft", "rs"),
                    labels = c("Constant", "Fluctuating", "Rate Summation")),
    param  = factor(param, levels = c("Tmin", "Topt", "Tmax"))
  )

cat("\nParameter point estimates with 95% bootstrap CIs:\n")
print(all_params)

# ==============================================================================
# 6. FIGURE (Shocket-style: TPC curves on top, parameter strip below)
# ==============================================================================

palette <- c(
  Constant         = unname(METHOD_PALETTE["constant"]),
  Fluctuating      = unname(METHOD_PALETTE["fluctuating"]),
  "Rate Summation" = "grey40"
)

base_theme_pooled <- function() {
  base_theme() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey60"))
}

# Top panel: TPC curves with a reduced legend (only "Rate Summation" shown —
# Constant/Fluctuating are self-evident from the coloured ribbons/points).
p_top <- ggplot() +
  geom_ribbon(data = curve_ci,
              aes(x = avg_temp, ymin = ct_lo, ymax = ct_hi, fill = "Constant"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = ct_curve,
            aes(x = avg_temp, y = predicted, colour = "Constant"),
            linewidth = 1, show.legend = FALSE) +
  geom_ribbon(data = curve_ci,
              aes(x = avg_temp, ymin = ft_lo, ymax = ft_hi, fill = "Fluctuating"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = ft_curve,
            aes(x = avg_temp, y = predicted, colour = "Fluctuating"),
            linewidth = 1, show.legend = FALSE) +
  geom_ribbon(data = curve_ci,
              aes(x = avg_temp, ymin = rs_lo, ymax = rs_hi, fill = "Rate Summation"),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(data = rs_curve,
            aes(x = avg_temp, y = predicted, linetype = "Rate Summation"),
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
  labs(x = "Temperature (°C)", y = "Development rate (day⁻¹)") +
  base_theme_pooled() +
  theme(legend.position = c(0.97, 0.97), legend.justification = c(1, 1),
        legend.background = element_blank())

# Bottom panel: thermal parameter strip (Shocket Fig. 3 D–F style)
p_bottom <- ggplot(all_params, aes(x = median, y = source, colour = source)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  facet_wrap(~ param, scales = "free_x", nrow = 1,
             labeller = labeller(param = c(Tmin = "T_min", Topt = "T_opt", Tmax = "T_max"))) +
  scale_colour_manual(values = palette, guide = "none") +
  scale_y_discrete(limits = rev(levels(all_params$source))) +
  labs(x = "Temperature (°C)", y = NULL) +
  base_theme_pooled() +
  theme(strip.text = element_text(face = "bold"))

p_briere_rs <- p_top / p_bottom + plot_layout(heights = c(2, 1))

print(p_briere_rs)
ggsave("figures/briere_l2a_pooled_RS.png", p_top, width = 9, height = 6, dpi = 300)

# ==============================================================================
# 7. RS BIAS QUANTIFICATION AT OBSERVED TEMPERATURES
# ==============================================================================

bias_summary <- ft_data %>%
  rename(ft_observed = rate_l2a) %>%
  mutate(
    rs_predicted = sapply(avg_temp, rs_predict, params = ct_params),
    bias_abs     = rs_predicted - ft_observed,
    bias_pct     = 100 * bias_abs / ft_observed
  )

cat("\n=== RS bias at observed FT temperatures (pooled) ===\n")
print(bias_summary %>% select(avg_temp, ft_observed, rs_predicted, bias_abs, bias_pct))

cat("\nFigure saved to figures/briere_l2a_pooled_RS.png\n")

# ── CT vs FT bias ────────────────────────────────────────────────────────────
# How much do observed constant-temperature rates exceed observed
# fluctuating-temperature rates at matched mean temperatures?
# (Distinct from RS bias, which compares the RS *prediction* to observed FT.)
ct_ft_bias_indiv <- rate_data %>%
  filter(avg_temp %in% ft_data$avg_temp, rate_l2a > 0) %>%
  group_by(avg_temp) %>%
  summarise(
    ct_n     = sum(method_temp == "constant"),
    ft_n     = sum(method_temp == "fluctuating"),
    ct_rate  = mean(rate_l2a[method_temp == "constant"]),
    ft_rate  = mean(rate_l2a[method_temp == "fluctuating"]),
    p_value  = tryCatch(
      t.test(rate_l2a[method_temp == "constant"],
             rate_l2a[method_temp == "fluctuating"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    bias_abs = ct_rate - ft_rate,
    bias_pct = 100 * bias_abs / ft_rate,
    sig = case_when(
      is.na(p_value)  ~ "n/a",
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

print(ct_ft_bias_indiv, n = Inf)
