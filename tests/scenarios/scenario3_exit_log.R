# tests/scenarios/scenario3_exit_log.R
# SCENARIO-3: SC-3 — exit_workflow() writes a valid, machine-parsable YAML log.
# TEST_PLAN.md §3 SCENARIO-3 (assertions 3.1–3.7).
# Run from project root: source("tests/scenarios/scenario3_exit_log.R")
# Requires: all R source files loadable; no Stan fit needed.

cat("=== SCENARIO-3: Exit Workflow and YAML Log ===\n\n")

# Load all workflow source files
source("R/wf_state.R")
source("R/context.R")
source("R/diagnostic_registry.R")
source("R/diagnostics.R")
source("R/offramps.R")
source("R/exit_workflow.R")

# --- Setup: known-bad dataset (binary, n=150, ~3% event rate) ---
set.seed(42)
dat_bad <- data.frame(
  y = c(rep(0L, 145L), rep(1L, 5L)),
  x = rnorm(150)
)

wf <- new_wf_state(mode = "practice", stage = "explore")
wf$declared_goal  <- "coefficient estimation"
wf$n_observations <- 150L
wf$event_rate     <- mean(dat_bad$y)

# --- Mock readline to auto-confirm ---
local_readline_yes <- function(code) {
  old <- options(rbayesflow_test_readline = "yes")
  on.exit(options(old))
  eval(code)
}

# Patch readline in base to return "yes" for duration of scenario
unlockBinding("readline", baseenv())
original_readline <- base::readline
assign("readline", function(prompt = "") "yes", envir = baseenv())

# =========================================================
# Assertion 3.1: assess_offramps returns non-empty list
# including "Firth penalized logistic" for binary rare event
# =========================================================
offramp_result <- assess_offramps(
  data         = dat_bad,
  outcome_var  = "y",
  outcome_type = "binary",
  goal         = "coefficient estimation",
  n            = 150
)
methods_3_1 <- sapply(offramp_result$alternatives, function(a) a$method)
stopifnot("3.1: off-ramp list empty" = length(offramp_result$alternatives) > 0)
stopifnot("3.1: firth_logistic not offered" = "firth_logistic" %in% methods_3_1)
cat("[3.1] PASS: off-ramp list non-empty; firth_logistic offered\n")

# =========================================================
# Assertion 3.2: exit_workflow() returns without error;
# exit.yaml written to tempdir
# =========================================================
yaml_path <- file.path(tempdir(), "exit.yaml")
wf2 <- tryCatch(
  exit_workflow(
    wf           = wf,
    method       = "firth_logistic",
    alternatives = c("logistic_bootstrap", "full_stan"),
    path         = yaml_path
  ),
  error = function(e) stop("3.2: exit_workflow() errored: ", e$message)
)
stopifnot("3.2: exit.yaml not written" = file.exists(yaml_path))
cat("[3.2] PASS: exit_workflow() returned without error; exit.yaml written\n")

# =========================================================
# Assertion 3.3: YAML parses; all required fields present
# =========================================================
log <- tryCatch(
  yaml::read_yaml(yaml_path),
  error = function(e) stop("3.3: yaml::read_yaml() failed: ", e$message)
)
required_fields <- c("method", "justification", "evidence_inspected",
                     "alternatives_considered", "user_confirmed",
                     "timestamp", "rbayesflow_version")
missing_fields <- setdiff(required_fields, names(log))
stopifnot("3.3: missing YAML fields" = length(missing_fields) == 0)
cat("[3.3] PASS: YAML parses; all required fields present\n")

# =========================================================
# Assertion 3.4: justification drawn from wf_state objects
# (contains sample_size and event_rate references)
# =========================================================
stopifnot("3.4: justification null or empty" =
  !is.null(log$justification) && nchar(log$justification) > 0)
stopifnot("3.4: justification contains n=" =
  grepl("n = 150", log$justification))
cat("[3.4] PASS: justification drawn from wf_state (n=150 referenced)\n")

# =========================================================
# Assertion 3.5: user_confirmed is TRUE
# =========================================================
stopifnot("3.5: user_confirmed not TRUE" = isTRUE(log$user_confirmed))
cat("[3.5] PASS: user_confirmed is TRUE\n")

# =========================================================
# Assertion 3.6: wf$stage set to "exit" after exit_workflow()
# =========================================================
stopifnot("3.6: stage not set to 'exit'" = identical(wf2$stage, "exit"))
cat("[3.6] PASS: wf$stage == 'exit'\n")

# =========================================================
# Assertion 3.7: exit at Phase 4 (after diagnostics run)
# evidence_inspected includes diagnostics_run: true
# =========================================================
wf_with_diag <- new_wf_state(mode = "practice", stage = "explore")
wf_with_diag$declared_goal  <- "coefficient estimation"
wf_with_diag$n_observations <- 150L
wf_with_diag$event_rate     <- mean(dat_bad$y)
# Simulate post-diagnostic state
wf_with_diag$diagnostics$passed         <- FALSE
wf_with_diag$diagnostics$failed_criteria <- c("rhat_max")
wf_with_diag$diagnostics$rhat_max       <- 1.045

yaml_path2 <- file.path(tempdir(), "exit_phase4.yaml")
wf3 <- exit_workflow(
  wf           = wf_with_diag,
  method       = "firth_logistic",
  alternatives = c("logistic_bootstrap", "full_stan"),
  path         = yaml_path2
)
log2 <- yaml::read_yaml(yaml_path2)
stopifnot("3.7: diagnostics_run not TRUE in evidence_inspected" =
  isTRUE(log2$evidence_inspected$diagnostics_run))
cat("[3.7] PASS: evidence_inspected$diagnostics_run == TRUE\n")

# --- Restore readline ---
assign("readline", original_readline, envir = baseenv())
lockBinding("readline", baseenv())

cat("\nSCENARIO-3: ALL 7 ASSERTIONS PASSED\n")
