# tests/scenarios/scenario4_offramps.R
# SCENARIO-4: Off-ramp equal weighting verification.
# TEST_PLAN.md §3 SCENARIO-4 (assertions 4.1–4.4).
# Run from project root: source("tests/scenarios/scenario4_offramps.R")

cat("=== SCENARIO-4: Off-Ramp Presentation ===\n\n")

source("R/wf_state.R")
source("R/context.R")
source("R/offramps.R")

set.seed(42)
dat_bad <- data.frame(
  y = c(rep(0L, 145L), rep(1L, 5L)),
  x = rnorm(150)
)

result <- assess_offramps(
  data         = dat_bad,
  outcome_var  = "y",
  outcome_type = "binary",
  goal         = "coefficient estimation",
  n            = 150
)

# 4.1: returns list with >= 2 alternatives
stopifnot("4.1: fewer than 2 alternatives" = length(result$alternatives) >= 2)
cat("[4.1] PASS: at least 2 alternatives returned\n")

# 4.2: no alternative is marked "recommended"
has_recommended <- any(sapply(result$alternatives, function(a) {
  "recommended" %in% names(a)
}))
stopifnot("4.2: at least one alternative marked 'recommended'" = !has_recommended)
cat("[4.2] PASS: no alternative marked 'recommended' (equal visual weight)\n")

# 4.3: event rate < 5% triggers Firth penalized logistic
methods <- sapply(result$alternatives, function(a) a$method)
stopifnot("4.3: firth_logistic not offered for rare event" =
  "firth_logistic" %in% methods)
cat("[4.3] PASS: firth_logistic offered for event rate < 5%\n")

# 4.4: full_stan always included
stopifnot("4.4: full_stan not included in alternatives" =
  "full_stan" %in% methods)
cat("[4.4] PASS: full_stan always included\n")

# Print the off-ramp output to confirm equal visual weight
cat("\n--- Off-ramp output (verify equal indentation) ---\n")
print_offramps(result)

cat("\nSCENARIO-4: ALL 4 ASSERTIONS PASSED\n")
