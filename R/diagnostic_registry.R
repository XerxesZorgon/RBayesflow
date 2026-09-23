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
  } else if (fam %in% c("beta_binomial", "betabinomial")) {
    "beta_binomial"
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
  outcome_var <- as.character(as.formula(formula(fit))[[2]])
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
  outcome_var <- as.character(as.formula(formula(fit))[[2]])
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
  outcome_var <- as.character(as.formula(formula(fit))[[2]])
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

# --- Beta-binomial family entry ---
# Triggered for bounded count outcomes (e.g. cannabis-use days in [0, 28])
# where binomial underdispersion has been identified and the user has moved
# to a beta-binomial family to model the extra variance. Added in v0.2.0
# after Chapter 18 (clinical trial case study) exposed the gap.
#
# Key checks:
#   1. Boundary-mass calibration: observed proportion at 0 and at `trials`
#      should fall inside the posterior predictive distribution.
#   2. Variance-to-mean ratio from posterior predictive draws, interpreted
#      on the same scale as the Poisson/NB checks but normalised by the
#      trials-adjusted maximum variance.
#
# NOTE: brms exposes the `phi` (precision) parameter for beta-binomial.
#   Larger phi => less overdispersion (pure binomial in the limit phi->inf).
#   We report the posterior median of phi as a diagnostic number.

DIAGNOSTIC_REGISTRY[["beta_binomial"]] <- function(fit, wf) {
  stopifnot(inherits(fit, "brmsfit"))

  # Outcome variable and trials constant
  outcome_var  <- as.character(as.formula(formula(fit))[[2]])
  y            <- fit$data[[outcome_var]]

  # Trials: the trials() addition sets `trials` as a column or a scalar.
  # brms stores it in fit$data under the column name used in trials(...).
  # Try to recover it; fall back to max(y) if unavailable.
  trials_val <- tryCatch({
    trials_col <- as.character(
      rlang::call_args(brms::brmsterms(formula(fit))$adforms$trials)[[1]]
    )
    if (trials_col %in% names(fit$data)) {
      fit$data[[trials_col]]
    } else {
      rep(as.integer(trials_col), length(y))   # scalar constant
    }
  }, error = function(e) rep(max(y), length(y)))

  n_trials <- if (length(unique(trials_val)) == 1L) trials_val[1L] else NA_integer_

  pp <- brms::posterior_predict(fit, ndraws = 100)

  # 1. Boundary mass: proportion of observations at 0 and at trials
  obs_zero_rate <- mean(y == 0)
  obs_max_rate  <- mean(y == n_trials, na.rm = TRUE)
  pred_zero_rates <- apply(pp, 1, function(r) mean(r == 0))
  pred_max_rates  <- apply(pp, 1, function(r) mean(r == n_trials, na.rm = TRUE))

  zero_pval <- mean(pred_zero_rates <= obs_zero_rate)
  max_pval  <- mean(pred_max_rates  <= obs_max_rate)

  warn <- character()
  if (zero_pval < 0.05 || zero_pval > 0.95) {
    warn <- c(warn, sprintf(
      "Boundary mass at 0 poorly calibrated: posterior predictive p-value = %.2f. ",
      zero_pval
    ))
  }
  if (!is.na(n_trials) && (max_pval < 0.05 || max_pval > 0.95)) {
    warn <- c(warn, sprintf(
      "Boundary mass at %d poorly calibrated: posterior predictive p-value = %.2f.",
      n_trials, max_pval
    ))
  }

  # 2. Posterior median of the phi (precision) parameter
  phi_draws <- tryCatch(
    as.numeric(posterior::as_draws_matrix(
      brms::as_draws(fit, variable = "phi")
    )),
    error = function(e) NA_real_
  )
  phi_median <- if (all(is.na(phi_draws))) NA_real_ else median(phi_draws, na.rm = TRUE)

  if (!is.na(phi_median) && phi_median > 50) {
    warn <- c(warn, sprintf(
      "phi median = %.1f (> 50): beta-binomial is approaching binomial. ",
      phi_median
    ))
  }

  # Plots
  stat_zero <- function(y_vec) mean(y_vec == 0)
  stat_max  <- if (!is.na(n_trials)) {
    function(y_vec) mean(y_vec == n_trials)
  } else {
    NULL
  }

  plots <- list(
    ppc_zero_mass = tryCatch(
      bayesplot::ppc_stat(y, pp, stat = "stat_zero"),
      error = function(e) NULL
    )
  )
  if (!is.null(stat_max)) {
    plots$ppc_max_mass <- tryCatch(
      bayesplot::ppc_stat(y, pp, stat = "stat_max"),
      error = function(e) NULL
    )
  }

  list(
    checks = list(
      n_trials         = n_trials,
      obs_zero_rate    = obs_zero_rate,
      obs_max_rate     = obs_max_rate,
      zero_pvalue      = zero_pval,
      max_pvalue       = max_pval,
      phi_median       = phi_median
    ),
    warnings = warn,
    plots    = plots
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


