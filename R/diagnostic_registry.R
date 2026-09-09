# R/diagnostic_registry.R
# DIAGNOSTIC_REGISTRY: named list of family-specific diagnostic functions.
# Each entry is a function(fit, wf) returning list(checks, warnings, plots).
# ADR-007: extensible registry architecture.

DIAGNOSTIC_REGISTRY <- list()

# family_key(): extract a canonical family name from a brmsfit object.
# Returns one of: "bernoulli", "poisson", "negbinomial",
#   "gaussian_hierarchical", "time_series", "unknown"
family_key <- function(fit) {
  stopifnot(inherits(fit, "brmsfit"))

  fam <- tryCatch(family(fit)$family, error = function(e) "unknown")

  # Check for time-series indicator: a time/date column in the data
  has_time_col <- any(c("time", "date", "year", "week", "month") %in%
                        tolower(names(fit$data)))

  # Check for random-effects (hierarchical) structure
  re_terms <- lme4::findbars(formula(fit))
  is_hierarchical <- !is.null(re_terms) && length(re_terms) > 0

  # Key resolution — order matters: time-series before gaussian_hierarchical
  key <- if (has_time_col && fam %in% c("gaussian", "student")) {
    "time_series"
  } else if (fam == "bernoulli") {
    "bernoulli"
  } else if (fam == "poisson") {
    "poisson"
  } else if (fam %in% c("negbinomial", "neg_binomial")) {
    "negbinomial"
  } else if (fam == "gaussian" && is_hierarchical) {
    "gaussian_hierarchical"
  } else if (fam %in% names(DIAGNOSTIC_REGISTRY)) {
    fam
  } else {
    warning("No family-specific diagnostic entry for family '", fam,
            "'. Using generic checks only.")
    "unknown"
  }

  key
}

# --- Bernoulli family entry (DESIGN.md §8.1) ---

DIAGNOSTIC_REGISTRY[["bernoulli"]] <- function(fit, wf) {
  # Extract outcome variable name from formula LHS
  outcome_var <- as.character(formula(fit)[[2]])
  y           <- fit$data[[outcome_var]]

  # Sample 100 rows from posterior predictive draws
  pp      <- brms::posterior_predict(fit, ndraws = 100)
  obs_p   <- mean(y)
  pred_p  <- mean(pp)
  ratio   <- if (pred_p > 0) obs_p / pred_p else NA_real_

  warn <- character()
  if (!is.na(obs_p) && obs_p < 0.05) {
    warn <- c(warn, sprintf(
      "Rare event detected (observed rate %.1f%%). Consider Firth penalized logistic or a penalized prior.",
      obs_p * 100
    ))
  }

  list(
    checks = list(
      event_rate        = obs_p,
      predicted_rate    = pred_p,
      calibration_ratio = ratio
    ),
    warnings = warn,
    plots    = list(
      calibration = tryCatch(
        bayesplot::ppc_bars(y, pp),
        error = function(e) NULL
      )
    )
  )
}
