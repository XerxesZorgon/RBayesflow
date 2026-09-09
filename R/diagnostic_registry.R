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
