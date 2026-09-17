# Shared configuration, palettes, and helper functions for the Aedes albopictus
# thermal-biology analysis pipeline.
#
# Contents:
#   - working directory
#   - population / treatment factor levels and colour palettes
#   - shared ggplot theme
#   - model-diagnostic helpers (overdispersion, cloglog inverse link)
#   - Brière-1 / Brière-2 thermal-performance-curve functions

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

setwd("~/Dati/Experiments")

# Base-R does not define `%||%` before R 4.4; define it defensively so every
# script can rely on it regardless of R version.
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# POPULATION / TREATMENT LEVELS AND COLOURS
# ==============================================================================

CITY_LEVELS <- c("Frankfurt", "Trento", "Palermo")
CITY_COLS   <- c(Frankfurt = "#5E4BB5", Trento = "#2E86AB", Palermo = "#E07A5F")

METHOD_LEVELS  <- c("constant", "fluctuating")
METHOD_PALETTE <- viridis::viridis(2, option = "C", begin = 0.2, end = 0.8)
names(METHOD_PALETTE) <- METHOD_LEVELS
METHOD_LTY <- c(constant = "solid",  fluctuating = "dashed")
METHOD_SHP <- c(constant = 16,       fluctuating = 17)

TEMP_BREAKS <- c(13, 18, 23, 28, 33)   # tested rearing temperatures (°C)

# ==============================================================================
# PLOT THEME
# ==============================================================================

# Default theme for individual trait/model figures (survival, development,
# wing length, population growth).
base_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey60"),
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = 13)
    )
}

# Variant used for the discrete-time hazard / cumulative-incidence figures
# (person_period_competing_risks, adult_longevity_cumulative_incidence).
base_theme_hazard <- function() {
  theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      strip.text       = element_text(face = "bold", size = 11)
    )
}

# ==============================================================================
# MODEL DIAGNOSTIC HELPERS
# ==============================================================================

# Pearson chi-square overdispersion ratio for a fitted GL(M)M.
# ratio >> 1 indicates overdispersion relative to the assumed variance function.
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  chi <- sum(rp^2)
  data.frame(
    Pearson_ChiSq = round(chi, 2),
    df            = rdf,
    ratio         = round(chi / rdf, 3),
    p_value       = round(pchisq(chi, df = rdf, lower.tail = FALSE), 4)
  )
}

# Inverse of the complementary log-log link, used to turn linear predictors
# from cloglog-linked discrete-time hazard models back into hazard probabilities.
inv_cloglog <- function(eta) 1 - exp(-exp(eta))

# ==============================================================================
# BRIÈRE THERMAL-PERFORMANCE-CURVE FUNCTIONS
# ==============================================================================

# Brière-1: r(T) = a * T * (T - Tmin) * sqrt(Tmax - T), clipped to [0, Inf)
briere1 <- function(T, a, Tmin, Tmax) {
  r <- a * T * (T - Tmin) * sqrt(pmax(Tmax - T, 0))
  pmax(r, 0)
}

# Brière-2: r(T) = a * T * (T - Tmin) * (Tmax - T)^(1/m)
briere2 <- function(T, a, Tmin, Tmax, m) {
  r <- a * T * (T - Tmin) * pmax(Tmax - T, 0)^(1 / m)
  pmax(r, 0)
}

# Analytical T_opt for Brière-1 (vertex of the curve on the temperature axis).
t_opt_briere1 <- function(Tmin, Tmax) {
  b <- 4 * Tmax + 3 * Tmin
  (b + sqrt(b^2 - 40 * Tmin * Tmax)) / 10
}