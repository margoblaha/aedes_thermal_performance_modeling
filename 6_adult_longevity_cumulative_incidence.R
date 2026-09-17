################################################################################
# longevity_cumulative_incidence.R
#
# Adapted from person_period_competing_risks.R for ADULT LONGEVITY analysis.
#
# DIFFERENCES FROM THE LARVAL/PUPAL COMPETING-RISKS ANALYSIS:
#   1. Single terminal event: ADULT DEATH. No competing risks here — every
#      adult eventually dies, and the only outcome is "days alive as adult".
#      So this is a standard discrete-time hazard model, NOT competing risks.
#   2. Constant-temperature data only: longevity was not recorded under
#      fluctuating regime, so method_temp is dropped from all models.
#   3. Population structure: city retained because the trait-level Poisson
#      GLMM showed Trento ≠ Frankfurt (β = +0.628, p = 0.003) with a
#      significant city × temperature interaction.
#   4. Sex retained as a fixed effect (males live ~7.6% shorter, p = 0.015).
#
# MODEL STRUCTURE (discrete-time cloglog hazard, mirrors larval framework):
#   log(day)              → Weibull-shaped baseline hazard
#   population            → city main effect
#   avg_temp + quad       → thermal response
#   sex                   → male vs female
#   log(day):avg_temp     → non-proportional hazards term
#   population:avg_temp   → does the thermal response differ between cities?
#   (1 | replic_uid)      → random intercept for replicate cup
#
# INTERPRETATION:
#   m_adult_death: HR > 1 → covariate increases the daily rate of adult
#                            mortality (conditional on still being alive).
#   Higher HR = shorter lifespan; lower HR = longer lifespan.
################################################################################

library(tidyverse)
library(glmmTMB)
library(broom.mixed)

setwd("~/Dati/Experiments")

constant_data <- readRDS("experimental_data.RDS") %>%
  filter(method_temp == "constant")

constant_data <- constant_data %>% 
  mutate(longevity=(days_larva_adult+adult_longevity))

CITY_LEVELS <- c("Frankfurt", "Trento", "Palermo")
CITY_COLS   <- c(Frankfurt = "#5E4BB5", Trento = "#2E86AB", Palermo = "#E07A5F")

adult_lifetable <- constant_data %>%
  filter(!is.na(adult_longevity)) %>%
  mutate(
    population = factor(city, levels = CITY_LEVELS),
    sex        = factor(sex,  levels = c("F", "M")),
    # +1 ensures that adults dying on day 0 still contribute a day-at-risk;
    # without this, log(day) would be undefined for those individuals.
    longevity_days = as.integer(longevity) + 1L
  ) %>%
  dplyr::select(
    individual_id = ID_individual,
    replic_uid,
    population,
    avg_temp,
    sex,
    longevity_days
  ) %>%
  filter(longevity_days >= 1)

cat("=== adult_lifetable diagnostics ===\n")
cat("n adults                :", nrow(adult_lifetable), "\n")
cat("Longevity range (days)  :", range(adult_lifetable$longevity_days), "\n")
cat("Median longevity        :", median(adult_lifetable$longevity_days), "\n")
cat("\nBy city × temperature:\n")
print(adult_lifetable %>%
        group_by(population, avg_temp) %>%
        summarise(n = n(), median_days = median(longevity_days),
                  .groups = "drop"))

# ==============================================================================
# 2.  BUILD person_period FOR ADULT DEATH
#
# Each adult contributes one row per day they were alive AS AN ADULT.
# died_pp = 1 on the final day (the day they died); 0 on every other day.
# ==============================================================================

MAX_DAY <- max(adult_lifetable$longevity_days)

adult_pp <- adult_lifetable %>%
  rowwise() %>%
  mutate(day = list(seq(1L, longevity_days))) %>%
  unnest(day) %>%
  mutate(died_pp = as.integer(day == longevity_days)) %>%
  ungroup()

cat("\n=== adult person-period diagnostics ===\n")
cat("Total person-days  :", nrow(adult_pp), "\n")
cat("Death events       :", sum(adult_pp$died_pp), "\n")
cat("Expected deaths    :", nrow(adult_lifetable),
    "(= n individuals)\n")

# Sanity check: each individual fires exactly once
sanity <- adult_pp %>%
  group_by(replic_uid, individual_id) %>%
  summarise(n_died = sum(died_pp), .groups = "drop")
cat("Individuals with != 1 death event (should be 0):\n")
print(sanity %>% filter(n_died != 1))

# ==============================================================================
# 3.  FIT DISCRETE-TIME HAZARD GLMM FOR ADULT DEATH
# ==============================================================================

cat("\n=== Fitting adult death hazard model ===\n")

m_adult_death <- glmmTMB(
  died_pp ~ log(day)
    + population
    + avg_temp + I(avg_temp^2)
    + sex
    + log(day):avg_temp
    + population:avg_temp
    + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = adult_pp
)

summary(m_adult_death)
# 
# # ==============================================================================
# # 4.  MODEL SELECTION: LRTs for key terms
# # ==============================================================================
# 
# cat("\n=== LRT: population effect ===\n")
# m_nopop <- glmmTMB(
#   died_pp ~ log(day) + avg_temp + I(avg_temp^2) + sex
#     + log(day):avg_temp + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = adult_pp
# )
# print(anova(m_nopop, m_adult_death))
# 
# cat("\n=== LRT: non-proportional hazards (log(day):avg_temp) ===\n")
# m_prop <- glmmTMB(
#   died_pp ~ log(day) + population + avg_temp + I(avg_temp^2) + sex
#     + population:avg_temp + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = adult_pp
# )
# print(anova(m_prop, m_adult_death))
# 
# cat("\n=== LRT: sex effect ===\n")
# m_nosex <- glmmTMB(
#   died_pp ~ log(day) + population + avg_temp + I(avg_temp^2)
#     + log(day):avg_temp + population:avg_temp + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = adult_pp
# )
# print(anova(m_nosex, m_adult_death))
# 
# cat("\n=== LRT: population × temperature interaction ===\n")
# m_nopopint <- glmmTMB(
#   died_pp ~ log(day) + population + avg_temp + I(avg_temp^2) + sex
#     + log(day):avg_temp + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = adult_pp
# )
# print(anova(m_nopopint, m_adult_death))

# ==============================================================================
# 5.  HAZARD RATIO TABLE
# ==============================================================================

hr_table <- tidy(m_adult_death, effects = "fixed", conf.int = TRUE) %>%
  filter(component == "cond") %>%
  mutate(
    HR    = round(exp(estimate), 3),
    HR_lo = round(exp(conf.low),  3),
    HR_hi = round(exp(conf.high), 3),
    p_lab = case_when(
      p.value < 0.001 ~ "<0.001",
      p.value < 0.01  ~ paste0(round(p.value, 3), " **"),
      p.value < 0.05  ~ paste0(round(p.value, 3), " *"),
      TRUE            ~ as.character(round(p.value, 3))
    )
  ) %>%
  dplyr::select(term, HR, HR_lo, HR_hi, p_lab)

cat("\n=== Hazard ratios (adult death) ===\n")
print(hr_table, n = Inf)

# ==============================================================================
# 5b. PAIRWISE POPULATION COMPARISONS AT EACH TEMPERATURE
#     Reuses the already-fitted m_adult_death — no refitting needed.
# ==============================================================================

library(emmeans)

pop_pairs <- emmeans(m_adult_death, ~ population | avg_temp,
                     at = list(avg_temp = sort(unique(adult_lifetable$avg_temp))),
                     type = "link") %>%
  pairs(adjust = "tukey") %>%
  as_tibble()

cat("\n=== Pairwise population comparisons by temperature (log-hazard scale) ===\n")
print(pop_pairs, n = Inf)

# ==============================================================================
# 6.  PREDICT CUMULATIVE INCIDENCE CURVES
#
# For a single-event survival problem:
#   S(t)        = prod_{s=1}^{t} (1 - h(s))      ← probability still alive
#   CI_death(t) = 1 - S(t)                        ← probability already dead
#
# This is just standard discrete-time survival — no competing risks needed
# because every adult eventually dies.
# ==============================================================================

# Predict for each combination of population × temperature × sex
# Females shown by default; sex can be flipped or faceted in the plot.
pred_grid <- expand.grid(
  day        = seq(1, MAX_DAY, by = 1),
  population = CITY_LEVELS,
  avg_temp   = sort(unique(adult_lifetable$avg_temp)),
  sex        = c("F", "M"),
  replic_uid = NA_character_
) %>%
  as_tibble() %>%
  mutate(
    population = factor(population, levels = CITY_LEVELS),
    sex        = factor(sex,        levels = c("F", "M"))
  )

pred_grid$hazard <- predict(
  m_adult_death,
  newdata          = pred_grid,
  type             = "response",
  allow.new.levels = TRUE
)

pred_ci <- pred_grid %>%
  group_by(population, avg_temp, sex) %>%
  arrange(day) %>%
  mutate(
    S        = cumprod(1 - hazard),
    CI_death = 1 - S
  ) %>%
  ungroup()

library(MASS)   # mvrnorm

set.seed(123)
N_SIM <- 500

# Fixed-effects-only formulas (strip the random-intercept term)
form_adult_death  <- delete.response(terms(lme4::nobars(formula(m_adult_death))))

X_death  <- model.matrix(form_adult_death,  data = pred_grid)
beta_death  <- fixef(m_adult_death)$cond
vcov_death  <- vcov(m_adult_death)$cond

sim_beta_death  <- mvrnorm(N_SIM, mu = beta_death,  Sigma = vcov_death)

inv_cloglog <- function(eta) 1 - exp(-exp(eta))  # cloglog inverse link

sim_ci_list <- vector("list", N_SIM)

for (s in seq_len(N_SIM)) {
  
  h_death_s  <- inv_cloglog(as.numeric(X_death  %*% sim_beta_death[s, ]))
  
  sim_grid_s <- pred_grid %>%
    mutate(h_death = h_death_s) %>%
    group_by(population, avg_temp, sex) %>%
    arrange(day) %>%
    mutate(
      S_prev    = cumprod(lag(1 - h_death, default = 1)),
      CI_death  = cumsum(h_death  * S_prev)
    ) %>%
    ungroup() %>%
    dplyr::select(population, avg_temp, sex, day, CI_death)
  
  sim_grid_s$sim <- s
  sim_ci_list[[s]] <- sim_grid_s
}

sim_ci <- bind_rows(sim_ci_list)

ci_bands <- sim_ci %>%
  group_by(population, avg_temp, sex, day) %>%
  summarise(
    CI_death_lo  = quantile(CI_death,  0.025),
    CI_death_hi  = quantile(CI_death,  0.975),
    .groups = "drop"
  )

pred_ci <- pred_ci %>%
  left_join(ci_bands, by = c("population", "avg_temp", "sex", "day"))

# ==============================================================================
# 7.  PLOTS
# ==============================================================================

base_theme <- function() {
  theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold", size = 11)
    )
}

# ── Plot 1: Cumulative incidence of adult death — facet by avg_temp, ──────────
#            colour-coded by population, NOT split by sex.
#
# There is no method_temp facet here (longevity was only recorded under
# constant temperature — see header notes), so the equivalent of the
# person_period_competing_risks.R faceting (avg_temp x method_temp) collapses
# to avg_temp alone. Sex is retained in the fitted model (m_adult_death) and
# in pred_grid/pred_ci, but the curve shown here is evaluated at sex = "F"
# (females), consistent with the reference-level convention already used
# elsewhere in this script; males sit at a uniformly higher hazard per the
# fitted sex effect (see hr_table).

p_ci_long <- pred_ci %>%
  mutate(avg_temp = factor(avg_temp))

# p_ci_pop <- p_ci_long %>%
#   filter(sex == "F")

plot_ci <- ggplot(p_ci_long,
                  aes(x = day, y = CI_death,
                      colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_death_lo, ymax = CI_death_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    sex ~ avg_temp,
    labeller = labeller(avg_temp = function(x) paste0(x, " °C"))
  ) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(suffix = ""),
    limits = c(0, 1)
  ) +
  labs(
    x   = "Days from hatching",
    y = "Cumulative incidence \n of death across the life span (%)",
    tag = "C"
  ) +
  base_theme() +
  theme(legend.position = "none")

print(plot_ci)
ggsave("figures/longevity_cumulative_incidence_ci.png",
       plot_ci, width = 11, height = 3.7, dpi = 300)
saveRDS(plot_ci, "figures/plot_c.rds")

# ── Combined A + B + C for the manuscript ──────────────────────────────────
# Requires plot_1a.rds and plot_1b.rds already saved by
# person_period_competing_risks.R (run that script first, or in the same
# session). Stacked vertically with patchwork's "/" operator; each panel's
# tag ("A"/"B"/"C") was already set via labs(tag = ...) in its own plot.

if (file.exists("figures/plot_1a.rds") && file.exists("figures/plot_1b.rds")) {

  library(patchwork)

  plot_1a <- readRDS("figures/plot_1a.rds")
  plot_1b <- readRDS("figures/plot_1b.rds")

  combined_ABC <- plot_1a / plot_1b / plot_ci +
    plot_layout(guides = "collect") +
    plot_annotation(theme = theme(legend.position = "bottom"))

  ggsave("figures/CI_combined_ABC.png",
         combined_ABC, width = 11, height = 16, dpi = 300)

} else {
  cat("\nNOTE: figures/plot_1a.rds and/or plot_1b.rds not found — run",
      "person_period_competing_risks.R first to generate the A/B panels,",
      "then re-run this script to produce the combined A+B+C figure.\n")
}

# ── Plot 2: HR dot-whisker plot (for transparency on effect sizes) ──────────

hr_plot_data <- hr_table %>%
  filter(!grepl("Intercept|log.day.", term))

plot_hr <- ggplot(hr_plot_data,
                  aes(x = HR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = HR_lo, xmax = HR_hi),
                 height = 0.25, linewidth = 0.7, colour = "#5E4BB5") +
  geom_point(size = 3, shape = 21, fill = "#5E4BB5", colour = "#5E4BB5") +
  scale_x_log10() +
  labs(
    x = "Hazard ratio for adult death (log scale)",
    y = NULL
  ) +
  base_theme() +
  theme(axis.text.y = element_text(size = 9))

print(plot_hr)
ggsave("figures/longevity_HR.png",
       plot_hr, width = 9, height = 5, dpi = 300)

cat("\n=== Done. Figures saved to figures/ ===\n")

# ---- sensitivity analysis ----
# Identify the single Trento individual at 13°C
trento_13 <- adult_lifetable %>%
  filter(population == "Trento", avg_temp == 13)

cat("The n=1 observation:\n")
print(trento_13)

# Refit excluding that individual
adult_pp_sens <- adult_pp %>%
  filter(!(population == "Trento" & avg_temp == 13))

m_adult_death_sens <- glmmTMB(
  died_pp ~ log(day)
  + population
  + avg_temp + I(avg_temp^2)
  + sex
  + log(day):avg_temp
  + population:avg_temp
  + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = adult_pp_sens
)

summary(m_adult_death_sens)

# Compare coefficient estimates
comparison <- bind_rows(
  tidy(m_adult_death,      effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Full"),
  tidy(m_adult_death_sens, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "Excluding Trento 13°C")
) %>%
  filter(component == "cond") %>%
  select(model, term, estimate, std.error, p.value, conf.low, conf.high) %>%
  pivot_wider(names_from = model,
              values_from = c(estimate, std.error, p.value, conf.low, conf.high))

print(comparison, n = Inf, width = Inf)
