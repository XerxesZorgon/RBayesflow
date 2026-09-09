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
  re_terms <- reformulas::findbars(formula(fit))
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
  } else if (fam != "unknown" && fam %in% names(DIAGNOSTIC_REGISTRY)) {
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

# --- Poisson family entry (DESIGN.md §8.2) ---

DIAGNOSTIC_REGISTRY[["poisson"]] <- function(fit, wf) {
  outcome_var <- as.character(formula(fit)[[2]])
  y           <- fit$data[[outcome_var]]
  pp          <- brms::posterior_predict(fit, ndraws = 100)

  # Variance-to-mean ratio from posterior predictive draws
  pp_means <- apply(pp, 1, mean)
  pp_vars  <- apply(pp, 1, var)
  ratio    <- mean(pp_vars / (pp_means + 1e-8))

  warn <- character()
  if (!is.na(ratio) && ratio > 2)
    warn <- c(warn, sprintf(
      "Overdispersion detected: posterior predictive variance-to-mean ratio = %.2f (> 2). Consider negative binomial family.",
      ratio
    ))
  if (!is.na(ratio) && ratio < 0.5)
    warn <- c(warn, sprintf(
      "Underdispersion detected: posterior predictive variance-to-mean ratio = %.2f (< 0.5). Check model specification.",
      ratio
    ))

  list(
    checks   = list(variance_to_mean_ratio = ratio),
    warnings = warn,
    plots    = list(
      ppc_mean = tryCatch(
        bayesplot::ppc_stat(y, pp, stat = "mean"),
        error = function(e) NULL
      )
    )
  )
}

# --- NegBinomial family entry (DESIGN.md §8.2) ---

DIAGNOSTIC_REGISTRY[["negbinomial"]] <- function(fit, wf) {
  outcome_var <- as.character(formula(fit)[[2]])
  y           <- fit$data[[outcome_var]]
  pp          <- brms::posterior_predict(fit, ndraws = 100)

  pp_means <- apply(pp, 1, mean)
  pp_vars  <- apply(pp, 1, var)
  ratio    <- mean(pp_vars / (pp_means + 1e-8))

  warn <- character()
  if (!is.na(ratio) && ratio > 2)
    warn <- c(warn, sprintf(
      "Overdispersion detected in NegBinomial fit: variance-to-mean ratio = %.2f (> 2). Check model specification.",
      ratio
    ))
  if (!is.na(ratio) && ratio < 0.5)
    warn <- c(warn, sprintf(
      "Underdispersion in NegBinomial fit: variance-to-mean ratio = %.2f. Consider Poisson or verifying the dispersion parameter.",
      ratio
    ))

  list(
    checks   = list(variance_to_mean_ratio = ratio),
    warnings = warn,
    plots    = list(
      ppc_mean = tryCatch(
        bayesplot::ppc_stat(y, pp, stat = "mean"),
        error = function(e) NULL
      )
    )
  )
}

# --- Gaussian hierarchical entry (DESIGN.md §8.3) ---

DIAGNOSTIC_REGISTRY[["gaussian_hierarchical"]] <- function(fit, wf) {
  # Identify the grouping variable from the first random-effects term
  re_terms   <- reformulas::findbars(formula(fit))
  group_var  <- if (length(re_terms) > 0) {
    as.character(re_terms[[1]][[3]])
  } else {
    NULL
  }

  n_groups <- if (!is.null(group_var) && group_var %in% names(fit$data)) {
    length(unique(fit$data[[group_var]]))
  } else {
    NA_integer_
  }

  warn <- character()
  if (!is.na(n_groups) && n_groups < 5) {
    warn <- c(warn, sprintf(
      "Few groups detected (n_groups = %d < 5). Centered parameterization may cause divergences. Consider calling refit_noncentered(wf).",
      n_groups
    ))
  }

  list(
    checks   = list(n_groups = n_groups, group_var = group_var),
    warnings = warn,
    plots    = list()
  )
}

# --- Time-series entry (DESIGN.md §8.4) ---

DIAGNOSTIC_REGISTRY[["time_series"]] <- function(fit, wf) {
  # Check ACF of residuals at lag 1
  res  <- tryCatch(residuals(fit)[, "Estimate"], error = function(e) NULL)
  acf1 <- if (!is.null(res) && length(res) > 5) {
    as.numeric(acf(res, plot = FALSE)$acf[2])
  } else {
    NA_real_
  }

  warn <- character()
  if (!is.na(acf1) && abs(acf1) > 0.3) {
    warn <- c(warn, sprintf(
      "Significant autocorrelation at lag 1 (ACF = %.3f, |ACF| > 0.3). Consider adding an AR(1) error structure in brms: autocor = ~ar(p = 1).",
      acf1
    ))
  }

  list(
    checks   = list(acf_lag1 = acf1),
    warnings = warn,
    plots    = list()
  )
}

# --- Unknown family fallback entry ---

DIAGNOSTIC_REGISTRY[["unknown"]] <- function(fit, wf) {
  list(
    checks   = list(),
    warnings = character(),
    plots    = list()
  )
}


