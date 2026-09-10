# RBayesflow — Test Plan (TEST_PLAN.md)

**Status:** Active
**Last updated:** 2026-09-10
**Owner:** John Peach
**Spec Kit role:** Defines the acceptance test strategy for all success criteria (SC-1 through SC-5) and the supporting unit tests. No new statistical code means no statistical unit tests. The primary strategy is scenario-based: run known-good and known-bad models through each mode/stage combination and verify the display contract.

---

## 1. Test Strategy

### 1.1 Primary: Scenario-Based Acceptance Tests

Each success criterion maps to one or more complete workflow scenarios executed end-to-end in an R session. A scenario passes when the observable outputs (console output, files written, plots displayed, YAML content) match the specification.

Scenarios are defined as R scripts in `tests/scenarios/`. They are run manually (not via `testthat` CI) because they require an interactive R session with Stan installed and operational. Each scenario produces a log file that records pass/fail per assertion.

### 1.2 Secondary: Unit Tests on `wf_state` and Related Functions

Functions that do not require fitting (constructors, display logic, exit log generation, registry dispatch) are tested with `testthat`. These tests run without Stan.

### 1.3 Not In Scope

- Statistical correctness of brms, bayesplot, loo, or any other dependency
- Posit Assistant behavior (LLM output is non-deterministic and not testable programmatically)
- Performance benchmarking
- Cross-platform CI (tests require Stan; CI setup is a future task)

---

## 2. Test Models

### 2.1 Known-Good Model — `examples/glmm_gaussian/`

**Dataset:** `sleepstudy` from the `lme4` package (180 observations; reaction time by days of sleep deprivation, grouped by subject). Chosen because: small, fast to fit, hierarchical, Gaussian family, well-understood posterior.

**Model:**
```r
fit_good <- brm(
  Reaction ~ Days + (Days | Subject),
  data   = lme4::sleepstudy,
  family = gaussian(),
  prior  = c(prior(normal(250, 50), class = Intercept),
             prior(normal(10, 5), class = b, coef = Days),
             prior(exponential(1), class = sigma)),
  backend = "cmdstanr",
  seed   = 42
)
```

Expected diagnostic outcome: Rhat < 1.01 for all parameters; bulk ESS > 400; tail ESS > 400; 0 divergences; BFMI > 0.3 on all chains; no max treedepth hit.

### 2.2 Known-Bad Model — `examples/bernoulli_rare_event/`

**Dataset:** Simulated data: n = 150; binary outcome with true event rate 3%; one continuous predictor. Chosen because: triggers the rare-event off-ramp warning, and with a highly misspecified prior will produce diagnostic failures.

**Model (intentionally misspecified prior to force diagnostic failure):**
```r
set.seed(99)
n <- 150
x <- rnorm(n)
p <- plogis(-4 + 0.5 * x)  # ~3% baseline event rate
y <- rbinom(n, 1, p)
dat_bad <- data.frame(y = y, x = x)

fit_bad <- brm(
  y ~ x,
  data    = dat_bad,
  family  = bernoulli(),
  prior   = prior(normal(0, 10), class = Intercept),  # diffuse, may cause sampling issues
  backend = "cmdstanr",
  seed    = 42
)
```

This model should trigger the Bernoulli rare-event warning from the diagnostic registry.

---

## 3. Acceptance Test Scenarios

### SCENARIO-1: SC-1 — Learn Mode Full Loop (Known-Good Model)

**Objective:** Verify SC-1. A user in `mode = "learn"` completes the full pipeline without writing code outside the templates.

**Setup:** Fresh R session; `source("R/init.R")`; `init_workflow(mode = "learn", stage = "explore")`

**Steps and assertions:**

| Step | Action | Assertion |
|---|---|---|
| 1.1 | Run Phase 1 template on `sleepstudy` | Off-ramp assessment runs; at least one alternative is offered; choice logged to audit trail |
| 1.2 | User selects "Full Stan (Bayesian)" | `wf_state$audit_trail` contains Phase 1 entry with `"bayesian_selected"` |
| 1.3 | Run Phase 2 template | Prior objects stored in `wf_state$priors_objects`; prior predictive draws stored in `wf_state$prior_pred_draws` |
| 1.4 | Run Phase 3 template (fit known-good model) | `fit_good` produced; `wf_state$fit_timestamp` set; `wf_state$parameterization` logged |
| 1.5 | `print(wf)` before Phase 4 | Output shows mode/stage/formula header only (no coefficients, no overlay — fit complete but diagnostics not yet run) |
| 1.6 | Run Phase 4 template | `wf_state$diagnostics$passed == TRUE`; `wf_state$diagnostics$rhat_max < 1.01` |
| 1.7 | `print(wf)` after clean diagnostics | Output shows **prior-vs-posterior overlay plot** (not coefficient table) + one-line health summary |
| 1.8 | Run Phase 5 template | `wf_state$ppc_complete == TRUE` |
| 1.9 | Run Phase 6 template | `wf_state$loo_complete == TRUE`; `wf_state$loo_table` is a data.frame |
| 1.10 | Render Phase 7 Quarto report | Report renders without error; all sections present; all claims traceable to `wf_state` |

**Pass criterion:** All 10 assertions pass.

---

### SCENARIO-2: SC-2 — Diagnostic Gate (Any Mode, Known-Bad Model)

**Objective:** Verify SC-2. A user cannot read coefficient output from a failed fit.

**Substeps run in both `mode = "learn"` and `mode = "practice"`.**

| Step | Action | Assertion |
|---|---|---|
| 2.1 | Fit known-bad model; run Phase 4 | `wf_state$diagnostics$passed == FALSE`; `wf_state$diagnostics$failed_criteria` is non-empty |
| 2.2 | `print(wf)` [mode = "learn"] | Output contains failure message + `wf$diagnose()` prompt; **no overlay plot shown** |
| 2.3 | `print(wf)` [mode = "practice"] | Output contains failure message + `wf$diagnose()` prompt; **no coefficient table shown** |
| 2.4 | Inspect `wf$diagnostics$acknowledged` | Value is `FALSE` |
| 2.5 | Call `wf$diagnose()`; enter "yes" at prompt | `wf$diagnostics$acknowledged` becomes `TRUE`; audit trail updated |
| 2.6 | `print(wf)` [mode = "practice"] after acknowledgment | Coefficient table shown with caveat warning |
| 2.7 | `print(wf)` [mode = "learn"] after acknowledgment | Failure message shown with acknowledgment note; overlay still suppressed (failed fit) |

**Pass criterion:** All 7 assertions pass for both modes.

---

### SCENARIO-3: SC-3 — Exit Workflow and YAML Log

**Objective:** Verify SC-3. `exit_workflow()` writes a valid, machine-parsable YAML log.

| Step | Action | Assertion |
|---|---|---|
| 3.1 | Run Phase 1 on known-bad dataset; off-ramps offered | Off-ramp list non-empty; includes "Firth penalized logistic" for binary rare event |
| 3.2 | Call `exit_workflow(wf, method = "firth_logistic")` | Returns without error; `exit.yaml` written to project root |
| 3.3 | Read `exit.yaml` | File parses with `yaml::read_yaml()`; all required fields present (method, justification, evidence_inspected, alternatives_considered, user_confirmed, timestamp, rbayesflow_version) |
| 3.4 | Check `justification` field | Field content is drawn from `wf_state` objects (sample_size, event_rate, declared_goal), not a user-typed string |
| 3.5 | Check `user_confirmed` field | Value is `true` |
| 3.6 | Check `wf_state$stage` after exit | Value is `"exit"` |
| 3.7 | Exit at Phase 4 (after fit, after diagnostics) | Same assertions 3.2–3.6 pass; `evidence_inspected` includes `diagnostics_run: true` |

**Pass criterion:** All 7 assertions pass.

---

### SCENARIO-4: Off-Ramp Presentation (Phase 1)

**Objective:** Verify that off-ramps are offered with equal visual weight before any fitting occurs.

| Step | Action | Assertion |
|---|---|---|
| 4.1 | `assess_offramps(data = dat_bad, outcome_type = "binary", goal = "coefficient estimation", n = 150)` | Returns list with ≥ 2 alternatives |
| 4.2 | Inspect printed output | All alternatives displayed at same indentation/weight; none marked "recommended" |
| 4.3 | Event rate < 5% | "Firth penalized logistic" or equivalent included in alternatives |
| 4.4 | "Full Stan" always included | `alternatives$method` vector includes `"full_stan"` |

---

### SCENARIO-5: Parameterization Transparency (Hierarchical Model)

**Objective:** Verify that brms's parameterization choice is detected and logged.

| Step | Action | Assertion |
|---|---|---|
| 5.1 | Fit `sleepstudy` model (centered by default in brms) | `wf_state$parameterization == "centered"` |
| 5.2 | Inspect audit trail | Entry with `action = "parameterization_logged"` present |
| 5.3 | Simulate divergences by fitting a pathological hierarchical model | `wf_state$diagnostics$n_divergences > 0`; workflow message includes recommendation to call `refit_noncentered(wf)` |

---

## 4. Unit Tests (`testthat`)

Located in `tests/testthat/`. Run with `testthat::test_dir("tests/testthat/")`. Do not require Stan.

### UT-1: `new_wf_state()` constructor

```r
test_that("new_wf_state validates mode and stage", {
  expect_error(new_wf_state(mode = "invalid", stage = "explore"))
  expect_error(new_wf_state(mode = "learn",   stage = "invalid"))
  wf <- new_wf_state(mode = "learn", stage = "explore")
  expect_s3_class(wf, "wf_state")
  expect_identical(wf$mode, "learn")
  expect_false(wf$diagnostics$acknowledged)
})
```

### UT-2: Display contract logic (no fit required)

```r
test_that("print.wf_state withholds coefficients on failed diagnostics", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$diagnostics$passed       <- FALSE
  wf$diagnostics$acknowledged <- FALSE
  out <- capture.output(print(wf))
  expect_false(any(grepl("Estimate", out)))   # no coefficient table
  expect_true(any(grepl("wf\\$diagnose", out))) # prompt present
})

test_that("print.wf_state shows coefficients after acknowledgment in practice mode", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$diagnostics$passed       <- FALSE
  wf$diagnostics$acknowledged <- TRUE
  wf$fit <- mock_brmsfit()   # minimal brmsfit stub with fixef() method
  out <- capture.output(print(wf))
  expect_true(any(grepl("Estimate", out)))
})
```

### UT-3: Exit log schema

```r
test_that("exit_workflow produces a valid YAML log with all required fields", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$declared_goal  <- "coefficient estimation"
  wf$n_observations <- 150
  wf$event_rate     <- 0.03

  withr::with_tempdir({
    log_path <- exit_workflow(wf, method = "firth_logistic",
                              alternatives = c("logistic_bootstrap", "full_stan"))
    log <- yaml::read_yaml(log_path)
    expect_true(log$user_confirmed)
    expect_equal(log$method, "firth_logistic")
    expect_true(!is.null(log$justification))
    expect_true(!is.null(log$evidence_inspected$sample_size))
    expect_true("full_stan" %in% log$alternatives_considered)
  })
})
```

### UT-4: Diagnostic registry dispatch

```r
test_that("run_diagnostics dispatches to the bernoulli registry entry", {
  registry_entry <- DIAGNOSTIC_REGISTRY[["bernoulli"]]
  expect_is(registry_entry, "function")
  # Smoke test with a minimal mock fit (stub posterior_predict)
  # Full test requires fitted model — deferred to scenario tests
})

test_that("unknown family falls back to unknown entry with warning", {
  expect_warning(
    family_key_result <- family_key(mock_fit(family = "custom_family")),
    regexp = "No family-specific"
  )
})
```

### UT-5: `export_context()` JSON schema

```r
test_that("export_context produces valid JSON with required fields", {
  wf <- new_wf_state(mode = "learn", stage = "explore")
  withr::with_tempdir({
    path <- export_context(wf, path = "wf_context.json")
    json <- jsonlite::fromJSON(path)
    expect_equal(json$mode, "learn")
    expect_equal(json$stage, "explore")
    expect_false(is.null(json$rbayesflow_version))
  })
})
```

**[v0.2.0]** UT-5 is extended to test the subfolder default path. Added case:

```r
test_that("export_context writes to data/<name>/ when called from an analysis subfolder [v0.2.0]", {
  withr::with_tempdir({
    # Simulate an RBayesflow project root with a DESCRIPTION file
    file.create("DESCRIPTION")
    dir.create(file.path("data", "my_analysis"), recursive = TRUE)
    file.create(file.path("data", "my_analysis", ".Rprofile"))
    withr::with_dir(file.path("data", "my_analysis"), {
      wf <- new_wf_state(mode = "learn", stage = "explore")
      path <- export_context(wf)   # no explicit path — use new default
      expect_true(file.exists("wf_context.json"))
      expect_equal(normalizePath(path), normalizePath("wf_context.json"))
    })
  })
})
```

### UT-6: `assess_offramps()` output contract

```r
test_that("assess_offramps returns at least 2 alternatives for binary rare event", {
  result <- assess_offramps(outcome_type = "binary", n = 150, event_rate = 0.03,
                            goal = "coefficient estimation")
  expect_gte(length(result$alternatives), 2)
  expect_true("full_stan" %in% sapply(result$alternatives, `[[`, "method"))
})
```

### UT-7: `rbf_install()` step vector **[v0.2.0]**

```r
test_that("rbf_install returns a named logical vector with 7 entries", {
  # Mock all system checks to return TRUE
  # (actual Stan install not triggered in unit test)
  with_mocked_bindings(
    check_cmdstan_toolchain = function(...) invisible(NULL),
    cmdstan_version         = function(...) "2.35.0",
    .package = "cmdstanr",
    {
      result <- rbf_install(dry_run = TRUE)   # dry_run skips install steps
      expect_type(result, "logical")
      expect_length(result, 7)
      expect_named(result)
    }
  )
})
```

### UT-8: `rbf_new()` folder structure **[v0.2.0]**

```r
test_that("rbf_new creates correct subfolder structure", {
  withr::with_tempdir({
    dir.create("data")
    file.create("DESCRIPTION")   # needed by rbf_analysis_path()
    rbf_new("my_analysis")
    expect_true(dir.exists(file.path("data", "my_analysis")))
    expect_true(file.exists(file.path("data", "my_analysis", ".Rprofile")))
    expect_true(file.exists(file.path("data", "my_analysis", "wf_context.json")))
    expect_true(file.exists(file.path("data", "my_analysis", "README.md")))
  })
})

test_that("rbf_new errors if analysis folder already exists", {
  withr::with_tempdir({
    file.create("DESCRIPTION")
    dir.create(file.path("data", "existing"), recursive = TRUE)
    expect_error(rbf_new("existing"), regexp = "already exists")
  })
})
```

### UT-9: `guide()` phase detection **[v0.2.0]**

```r
test_that("guide detects Phase 1 for a fresh wf_state", {
  wf <- new_wf_state(mode = "learn", stage = "explore")
  out <- capture.output(guide(wf))
  expect_true(any(grepl("Phase 1", out)))
  expect_true(any(grepl("phase1_exploration", out)))
})

test_that("guide detects diagnostic block when diagnostics failed and unacknowledged", {
  wf <- new_wf_state(mode = "practice", stage = "explore")
  wf$fit_timestamp             <- Sys.time()
  wf$diagnostics$passed        <- FALSE
  wf$diagnostics$acknowledged  <- FALSE
  out <- capture.output(guide(wf))
  expect_true(any(grepl("diagnose", out)))
})

test_that("guide detects exit stage", {
  wf <- new_wf_state(mode = "practice", stage = "exit")
  out <- capture.output(guide(wf))
  expect_true(any(grepl("exit", tolower(out))))
})
```

---

## 5. Test Coverage Targets

| Component | Coverage target | Test type |
|---|---|---|
| `new_wf_state()` | 100% of validation branches | Unit |
| `print.wf_state()` | All 4 display-contract branches | Unit |
| `diagnose.wf_state()` | Acknowledged / not-acknowledged paths | Unit |
| `exit_workflow()` | YAML schema compliance | Unit |
| `assess_offramps()` | All decision-matrix rows | Unit |
| `export_context()` | JSON schema; subfolder path default **[v0.2.0]** | Unit |
| `DIAGNOSTIC_REGISTRY` entries | Bernoulli, Poisson, Hierarchical, Time-series | Unit (smoke) + Scenario |
| Full pipeline (SC-1) | `mode = "learn"`, `stage = "explore"` | Scenario |
| Diagnostic gate (SC-2) | `mode = "learn"` and `mode = "practice"` | Scenario |
| Exit log (SC-3) | Phase 1 exit and Phase 4 exit | Scenario |
| **[v0.2.0]** `rbf_install()` | 7-step vector; dry_run path | Unit |
| **[v0.2.0]** `rbf_new()` | Folder structure; duplicate error | Unit |
| **[v0.2.0]** `guide()` | All 10 phase-detection branches | Unit |
| **[v0.2.0]** `rbf_analysis_path()` | Inside subfolder; project root fallback | Unit |
| **[v0.2.0]** User-guide Markdown | All 5 files render without error | Manual |

---

*End of TEST_PLAN.md*
