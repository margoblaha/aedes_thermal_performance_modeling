################################################################################
# person_period_competing_risks.R
#
# Extends person_period_analysis_fixed.R to model eclosion and larval/pupal
# death as COMPETING EVENTS using cause-specific discrete-time hazard GLMMs.
#
# WHY glmmTMB INSTEAD OF coxme?
#   coxme fits continuous-time Cox models with random effects but cannot
#   include the log(day) x temperature non-proportional hazards term that
#   was decisive in the original analysis (ΔAIC = -100). The person-period
#   glmmTMB framework already in this script handles:
#     - random intercepts for replicate cup
#     - non-proportional hazards via log(day):avg_temp
#     - cause-specific competing risks (two separate outcome columns)
#   Switching to coxme would lose all of these.
#
# COMPETING RISKS APPROACH: cause-specific hazards
#   - Two outcomes: death (larva/pupa dies) and eclosion (reaches adulthood)
#   - For each outcome, the other is treated as censored on its event day
#   - This correctly answers: "what drives mortality?" and "what drives
#     eclosion?" separately, while accounting for the competing nature of
#     the two events
#   - Alternative (Fine-Gray subdistribution) would be used if you wanted
#     cumulative incidence functions for absolute risk; cause-specific is
#     preferred here because the biological question is about mechanisms
#
# INTERPRETATION:
#   m_death:  HR > 1 means covariate increases the daily rate of dying
#             (conditional on not having died or eclosed yet)
#   m_eclose: HR > 1 means covariate increases the daily rate of eclosing
#             (conditional on not having died or eclosed yet)
#   Together: if temperature increases m_eclose hazard but not m_death,
#             warming accelerates development without increasing mortality.
#             If treatment increases m_death AND decreases m_eclose,
#             fluctuating temperatures are both more lethal AND slower.
################################################################################

library(tidyverse)
library(glmmTMB)
library(broom.mixed)
library(survival)    # for Surv() diagnostics only

setwd("~/Dati/Experiments")
all_data <- readRDS("experimental_data.RDS")

CITY_LEVELS <- c("Frankfurt", "Trento", "Palermo")
CITY_COLS   <- c(Frankfurt = "#5E4BB5", Trento = "#2E86AB", Palermo = "#E07A5F")

# ==============================================================================
# 1.  BUILD survival_raw WITH COMPETING EVENT INDICATOR
# ==============================================================================

survival_raw <- all_data %>%
  mutate(
    days_to_adult   = as.integer(adult_emergence - start),
    days_to_pupa_d  = as.integer(pupa_death      - start),
    days_to_larva_d = as.integer(larva_death      - start),

    death_day = case_when(
      !is.na(days_to_adult)   & days_to_adult   > 0 ~ days_to_adult,
      !is.na(days_to_pupa_d)  & days_to_pupa_d  > 0 ~ days_to_pupa_d,
      !is.na(days_to_larva_d) & days_to_larva_d > 0 ~ days_to_larva_d,
      TRUE ~ NA_integer_
    ),

    # Competing event indicator:
    #   0 = censored (not used here — all outcomes resolved)
    #   1 = died as larva or pupa
    #   2 = eclosed (reached adulthood)
    event = case_when(
      adult == 1 ~ 2L,   # eclosed
      TRUE       ~ 1L    # died as larva or pupa
    ),

    population = as.character(city)
  ) %>%
  filter(!is.na(death_day), death_day > 0) %>%
  dplyr::select(
    individual_id = ID_individual,
    population,
    avg_temp,
    method_temp,
    replic_uid,
    death_day,
    event,          # 1 = died, 2 = eclosed
    sex,
    adult,
    generation
  )

cat("=== survival_raw: competing events ===\n")
cat("n individuals :", nrow(survival_raw), "\n")
cat("Died (event=1):", sum(survival_raw$event == 1), "\n")
cat("Eclosed (event=2):", sum(survival_raw$event == 2), "\n\n")

# ==============================================================================
# 2.  BUILD person_period WITH TWO OUTCOME COLUMNS
#
#   - died_pp:    1 on final day if the individual died (event == 1), else 0
#   - eclosed_pp: 1 on final day if the individual eclosed (event == 2), else 0
#
#   Both are 0 on every non-final day.
#   On the final day, exactly one of the two fires (never both).
#   This is the standard cause-specific person-period encoding.
# ==============================================================================

MAX_DAY <- 150L

person_period <- survival_raw %>%
  mutate(death_day = pmin(death_day, MAX_DAY)) %>%
  rowwise() %>%
  mutate(day = list(seq(1L, death_day))) %>%
  unnest(day) %>%
  mutate(
    died_pp    = as.integer(day == death_day & event == 1L),
    eclosed_pp = as.integer(day == death_day & event == 2L)
  ) %>%
  ungroup()

cat("=== person_period diagnostics ===\n")
cat("Total rows         :", nrow(person_period), "\n")
cat("Death events       :", sum(person_period$died_pp), "\n")
cat("Eclosion events    :", sum(person_period$eclosed_pp), "\n")
cat("Expected deaths    :", sum(survival_raw$event == 1),
    "(= n individuals who died)\n")
cat("Expected eclosions :", sum(survival_raw$event == 2),
    "(= n individuals who eclosed)\n\n")

# Sanity: each individual fires exactly once across the two outcomes
sanity <- person_period %>%
  group_by(replic_uid, individual_id) %>%
  summarise(
    n_died    = sum(died_pp),
    n_eclosed = sum(eclosed_pp),
    .groups   = "drop"
  ) %>%
  mutate(total_events = n_died + n_eclosed)

cat("Individuals with != 1 total event (should be 0):\n")
print(sanity %>% filter(total_events != 1))

# ==============================================================================
# 3.  CAUSE-SPECIFIC DISCRETE-TIME HAZARD GLMMs
#
#   Both models share the same formula structure and the same person_period
#   data. The only difference is the outcome column.
#
#   For m_death:  died_pp    is the outcome; on eclosion days died_pp = 0
#                 (eclosure treated as censored in cause-specific framing)
#   For m_eclose: eclosed_pp is the outcome; on death days eclosed_pp = 0
#                 (death treated as censored in cause-specific framing)
#
#   MODEL STRUCTURE (same as original primary model m_nonprop / m_final):
#     log(day)              → Weibull baseline hazard
#     population            → city main effect
#     avg_temp + quad       → thermal response (quadratic)
#     method_temp           → treatment main effect
#     avg_temp:method_temp  → does treatment shift the thermal response?
#     log(day):avg_temp     → non-proportional: does temp effect change over time?
#     population*avg_temp   → does population shift the thermal response?
#     (1 | replic_uid)      → random intercept for replicate cup
# ==============================================================================

cat("\n=== Fitting cause-specific models ===\n")

# ── Model for DEATH ───────────────────────────────────────────────────────────
m_death <- glmmTMB(
  died_pp ~ log(day)
  + population
  + avg_temp + I(avg_temp^2)
  + method_temp
  + avg_temp:method_temp
  + I(avg_temp^2):method_temp
  + log(day):avg_temp
  + population:avg_temp
  + population:I(avg_temp^2)
  + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = person_period
)

# ── Model for ECLOSION ────────────────────────────────────────────────────────
m_eclose <- glmmTMB(
  eclosed_pp ~ log(day)
  + population
  + avg_temp + I(avg_temp^2)
  + method_temp
  + avg_temp:method_temp
  + I(avg_temp^2):method_temp
  + log(day):avg_temp
  + population:avg_temp
  + population:I(avg_temp^2)
  + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = person_period
)

# ==============================================================================
# 4.  MODEL SELECTION: test key terms by LRT for each cause
# ==============================================================================

m_death_nopop <- glmmTMB(died_pp
  ~ log(day) + avg_temp + I(avg_temp^2) + method_temp +
    avg_temp:method_temp + I(avg_temp^2):method_temp +
    log(day):avg_temp + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = person_period)

anova(m_death_nopop, m_death)

m_eclose_nopop <- glmmTMB(eclosed_pp
  ~ log(day) + avg_temp + I(avg_temp^2) + method_temp +
    avg_temp:method_temp + I(avg_temp^2):method_temp +
    log(day):avg_temp + (1 | replic_uid),
  family = binomial(link = "cloglog"),
  data   = person_period)

anova(m_eclose_nopop, m_eclose) # population city effect non significant!

#LRT: non-proportional hazards (log(day):avg_temp)
# m_death_prop <- glmmTMB(died_pp
#   ~ log(day) + population + avg_temp + I(avg_temp^2) + method_temp +
#     avg_temp:method_temp + I(avg_temp^2):method_temp +
#     population:avg_temp + population:I(avg_temp^2) + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = person_period)
# 
# anova(m_death_prop, m_death)

# m_eclose_prop <- glmmTMB(eclosed_pp
#   ~ log(day) + population + avg_temp + I(avg_temp^2) + method_temp +
#     avg_temp:method_temp + I(avg_temp^2):method_temp +
#     population:avg_temp + population:I(avg_temp^2) + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = person_period)
# 
# anova(m_eclose_prop, m_eclose)

# cat("\n=== LRT: treatment effect ===\n")
# m_death_notreat <- glmmTMB(died_pp
#    ~ log(day) + population + avg_temp + I(avg_temp^2) +
#     log(day):avg_temp +
#     population:avg_temp + population:I(avg_temp^2) + (1 | replic_uid),
#    family = binomial(link = "cloglog"),
#    data   = person_period)
# cat("Death — treatment LRT:\n"); print(anova(m_death_notreat, m_death))
# 
# m_eclose_notreat <- glmmTMB(eclosed_pp
#   ~ log(day) + population + avg_temp + I(avg_temp^2) +
#     log(day):avg_temp +
#     population:avg_temp + population:I(avg_temp^2) + (1 | replic_uid),
#   family = binomial(link = "cloglog"),
#   data   = person_period)
# cat("Eclosion — treatment LRT:\n"); print(anova(m_eclose_notreat, m_eclose))

# ==============================================================================
# 5.  HAZARD RATIO TABLE (both outcomes side by side)
# ==============================================================================

tidy_cause <- function(mod, cause_label) {
  tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    filter(component == "cond") %>%
    mutate(
      cause  = cause_label,
      HR     = round(exp(estimate), 3),
      HR_lo  = round(exp(conf.low),  3),
      HR_hi  = round(exp(conf.high), 3),
      p_lab  = case_when(
        p.value < 0.001 ~ "<0.001",
        p.value < 0.01  ~ paste0(round(p.value, 3), " **"),
        p.value < 0.05  ~ paste0(round(p.value, 3), " *"),
        TRUE            ~ as.character(round(p.value, 3))
      )
    ) %>%
    dplyr::select(cause, term, HR, HR_lo, HR_hi, p_lab)
}

hr_table <- bind_rows(
  tidy_cause(m_death,  "Death"),
  tidy_cause(m_eclose_nopop, "Eclosion")
)

cat("\n=== Hazard ratio table: Death vs Eclosion ===\n")
print(hr_table, n = Inf)

# ==============================================================================
# 6.  PREDICT CAUSE-SPECIFIC HAZARD CURVES
# ==============================================================================

pred_grid <- expand.grid(
  day         = seq(1, 100, by = 1),
  population  = CITY_LEVELS,
  avg_temp    = c(13, 18, 23, 28, 33),
  method_temp = c("constant", "fluctuating"),
  replic_uid  = NA_character_
) %>%
  as_tibble() %>%
  mutate(
    population  = factor(population,  levels = CITY_LEVELS),
    method_temp = factor(method_temp, levels = c("constant", "fluctuating"))
  )

pred_grid$h_death   <- predict(m_death,  newdata = pred_grid,
                                type = "response", allow.new.levels = TRUE)
pred_grid$h_eclose  <- predict(m_eclose_nopop, newdata = pred_grid,
                                type = "response", allow.new.levels = TRUE)

# Overall survival = probability of neither event having occurred
# S(t) = prod_{s=1}^{t} [1 - h_death(s) - h_eclose(s)]
# Cumulative incidence of each cause:
#   CI_death(t)   = sum_{s=1}^{t} h_death(s)  * S(s-1)
#   CI_eclose(t)  = sum_{s=1}^{t} h_eclose(s) * S(s-1)

pred_ci <- pred_grid %>%
  group_by(population, avg_temp, method_temp) %>%
  arrange(day) %>%
  mutate(
    # Probability of surviving both events to start of each day
    S_prev    = cumprod(lag(1 - h_death - h_eclose, default = 1)),
    # Overall survival to end of each day
    S         = cumprod(1 - h_death - h_eclose),
    # Cause-specific cumulative incidence
    CI_death  = cumsum(h_death  * S_prev),
    CI_eclose = cumsum(h_eclose * S_prev)
  ) %>%
  ungroup()

# ==============================================================================
# 6b. SIMULATION-BASED CONFIDENCE INTERVALS FOR CUMULATIVE INCIDENCE
#
# WHY NOT THE DELTA METHOD:
#   CI_death(t) and CI_eclose(t) are path-dependent transforms of the fitted
#   hazards (cumulative products/sums over every prior day). Propagating a
#   per-day SE analytically through that chain is intractable.
#
# APPROACH (parametric bootstrap of the fixed effects):
#   Draw N_SIM coefficient vectors from the asymptotic sampling distribution
#   of each cause-specific model — MVN(fixef, vcov) — and, for EACH draw,
#   rebuild the entire hazard -> S_prev -> CI curve. Because one draw is
#   reused across all days within a curve, the day-to-day correlation induced
#   by shared parameter uncertainty is preserved (unlike simulating each day
#   independently). The 2.5th/97.5th percentiles across draws, at each day,
#   give the uncertainty band. No refitting of the mixed model is required.
#
# ASSUMPTION: m_death and m_eclose_nopop are simulated independently. Both
#   are fit on the same person_period data, so their estimates are not
#   literally independent, but treating them as such is a standard
#   simplification for visualisation-level bands and errs conservative
#   (slightly wider bands than a fully joint simulation would give).
# ==============================================================================

library(MASS)   # mvrnorm

set.seed(123)
N_SIM <- 500

# Fixed-effects-only formulas (strip the random-intercept term)
form_death  <- delete.response(terms(lme4::nobars(formula(m_death))))
form_eclose <- delete.response(terms(lme4::nobars(formula(m_eclose_nopop))))

X_death  <- model.matrix(form_death,  data = pred_grid)
X_eclose <- model.matrix(form_eclose, data = pred_grid)

beta_death  <- fixef(m_death)$cond
vcov_death  <- vcov(m_death)$cond
beta_eclose <- fixef(m_eclose_nopop)$cond
vcov_eclose <- vcov(m_eclose_nopop)$cond

sim_beta_death  <- mvrnorm(N_SIM, mu = beta_death,  Sigma = vcov_death)
sim_beta_eclose <- mvrnorm(N_SIM, mu = beta_eclose, Sigma = vcov_eclose)

inv_cloglog <- function(eta) 1 - exp(-exp(eta))  # cloglog inverse link

sim_ci_list <- vector("list", N_SIM)

for (s in seq_len(N_SIM)) {

  h_death_s  <- inv_cloglog(as.numeric(X_death  %*% sim_beta_death[s, ]))
  h_eclose_s <- inv_cloglog(as.numeric(X_eclose %*% sim_beta_eclose[s, ]))

  sim_grid_s <- pred_grid %>%
    mutate(h_death = h_death_s, h_eclose = h_eclose_s) %>%
    group_by(population, avg_temp, method_temp) %>%
    arrange(day) %>%
    mutate(
      S_prev    = cumprod(lag(1 - h_death - h_eclose, default = 1)),
      CI_death  = cumsum(h_death  * S_prev),
      CI_eclose = cumsum(h_eclose * S_prev)
    ) %>%
    ungroup() %>%
    dplyr::select(population, avg_temp, method_temp, day, CI_death, CI_eclose)

  sim_grid_s$sim <- s
  sim_ci_list[[s]] <- sim_grid_s
}

sim_ci <- bind_rows(sim_ci_list)

ci_bands <- sim_ci %>%
  group_by(population, avg_temp, method_temp, day) %>%
  summarise(
    CI_death_lo  = quantile(CI_death,  0.025),
    CI_death_hi  = quantile(CI_death,  0.975),
    CI_eclose_lo = quantile(CI_eclose, 0.025),
    CI_eclose_hi = quantile(CI_eclose, 0.975),
    .groups = "drop"
  )

pred_ci <- pred_ci %>%
  left_join(ci_bands, by = c("population", "avg_temp", "method_temp", "day"))

cat("\n=== Simulation-based CIs added to pred_ci (N_SIM =", N_SIM, ") ===\n")

# ==============================================================================
# 7.  PLOTS
# ==============================================================================

# temp_palette <- setNames(
#   viridis::viridis(5, option = "D"),
#   as.character(c(13, 18, 23, 28, 33))
# )

base_theme <- function() {
  theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold", size = 11)
    )
}

# ── Plot 1a: Cumulative incidence of DEATH (larval/pupal) ─────────────────────
# Faceted by method_temp (rows) x avg_temp (cols); populations colour-coded.
# Ribbon = 95% simulation-based CI (see section 6b).

p_death_long <- pred_ci %>%
  dplyr::select(population, avg_temp, method_temp, day,
         CI_death, CI_death_lo, CI_death_hi) %>%
  mutate(
    avg_temp    = factor(avg_temp),
    method_temp = factor(method_temp,
                         levels = c("constant", "fluctuating"),
                         labels = c("Constant", "Fluctuating"))
  )

plot_1a <- ggplot(p_death_long,
                  aes(x = day, y = CI_death,
                      colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_death_lo, ymax = CI_death_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    method_temp ~ avg_temp,
    labeller = labeller(avg_temp = function(x) paste0(x, " °C"))
  ) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(
    labels = scales::label_percent(suffix = ""),
    limits = c(0, 1)
  ) +
  labs(
    x        = "Days from hatching",
    y        = "Cumulative incidence \n  of larval and pupal death (%)",
    tag = "B"
  ) +
  base_theme()+
  theme(legend.position = "none")

print(plot_1a)
ggsave("figures/lp death facet city.png", plot_1a, width = 11, height = 6, dpi = 300)
saveRDS(plot_1a, "figures/plot_1a.rds")


# ── Plot 1b: Cumulative incidence of ECLOSION ──────────────────────────────────
# Faceted by method_temp (rows) x avg_temp (cols); populations colour-coded.
# Ribbon = 95% simulation-based CI (see section 6b).

p_eclose_long <- pred_ci %>%
  dplyr::select(population, avg_temp, method_temp, day,
         CI_eclose, CI_eclose_lo, CI_eclose_hi) %>%
  mutate(
    avg_temp    = factor(avg_temp),
    method_temp = factor(method_temp,
                         levels = c("constant", "fluctuating"),
                         labels = c("Constant", "Fluctuating"))
  )

plot_1b <- ggplot(p_eclose_long,
                  aes(x = day, y = CI_eclose,
                      colour = population, fill = population)) +
  geom_ribbon(aes(ymin = CI_eclose_lo, ymax = CI_eclose_hi),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_grid(
    method_temp ~ avg_temp,
    labeller = labeller(avg_temp = function(x) paste0(x, " °C"))
  ) +
  scale_colour_manual(values = CITY_COLS, name = "Population") +
  scale_fill_manual(values = CITY_COLS, name = "Population") +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(
    labels = scales::label_percent(suffix = ""),
    limits = c(0, 1)
  ) +
  labs(
    x        = "Days from hatching",
    y        = "Cumulative incidence \n  of adult emergence (%)",
    tag = "A"
  ) +
  base_theme()+ 
  theme(legend.position = "none")

print(plot_1b)
ggsave("figures/adult emergence cumulative incidence.png", plot_1b, width = 11, height = 6, dpi = 300)
saveRDS(plot_1b, "figures/plot_1b.rds")

all_data %>%
  filter(avg_temp == 28) %>%
  group_by(city, method_temp) %>%
  summarise(
    n_total  = n(),
    n_adults = sum(adult, na.rm = TRUE),
    pct      = round(100 * n_adults / n_total, 1),
    .groups  = "drop"
  )

# ── Combined 1a + 1b with patchwork ──────────────────────────────────────────
# library(patchwork)
# 
# combined <- plot_1a / plot_1b +
#   plot_layout(guides = "collect") +
#   plot_annotation(
#     #title = "Cause-specific cumulative incidence by population and treatment",
#     theme = theme(plot.title = element_text(face = "bold", size = 13),
#                   legend.position = "bottom")
#   )
# 
# ggsave("figures/CI_combined.png",
#        combined, width = 11, height = 11, dpi = 300)

# ── Plot B: HR comparison — same predictor, both causes ──────────────────────
# Dot-whisker plot: if death HR and eclosion HR point in opposite directions
# for a covariate, that covariate has antagonistic effects on the two processes

hr_plot_data <- hr_table %>%
  filter(!grepl("Intercept|log.day.", term)) %>%
  mutate(
    term = factor(term),
    
    outcome = recode(
      cause,
      "Death"     = "Larval or Pupal death",
      "Eclosion" = "Adult emergence"
    ),
    
    outcome = factor(
      outcome,
      levels = c("Larval or Pupal death", "Adult emergence")
    )
  )


# plot_hr <- ggplot(hr_plot_data,
#                   aes(x = HR, y = term, colour = outcome)) +
#   geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
#   geom_errorbarh(
#     aes(xmin = HR_lo, xmax = HR_hi),
#     height = 0.25,
#     linewidth = 0.7,
#     position = position_dodge(width = 0.5)
#   ) +
#   
#   geom_point(
#     size = 3,
#     shape = 21,
#     aes(fill = outcome),
#     position = position_dodge(width = 0.5)
#   ) +
#   geom_point(size = 3, shape = 21, aes(fill = outcome),
#              position = position_dodge(width = 0.3)) +  # Dodging points horizontally
#   scale_colour_manual(
#     values = c("Larval or Pupal death" = "#C0392B",
#                "Adult emergence" = "#1B7837"),
#     name = "Outcome"
#   ) +
#   scale_fill_manual(
#     values = c("Larval or Pupal death" = "#C0392B",
#                "Adult emergence" = "#1B7837"),
#     name = "Outcome"
#   ) +
#   scale_x_log10() +
#   labs(
#     x = "Hazard ratio (log scale)",
#     y = NULL,
#     #title = "Cause-specific hazard ratios for mortality and adult emergence",
#     #subtitle = "HR > 1 indicates an increased event rate; HR < 1 indicates a decreased event rate"
#   ) +
#   base_theme() +
#   theme(axis.text.y = element_text(size = 9))
# 
# print(plot_hr)
# ggsave("figures/hr_death_vs_eclosion.png",
#        plot_hr, width = 10, height = 7, dpi = 300)

# cat("\n=== Done. Figures saved to figures/ ===\n")
