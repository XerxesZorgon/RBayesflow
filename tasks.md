# RBayesflow — tasks.md

**Generated:** 2026-09-08 (v0.1.0) | **Updated:** 2026-09-10 (v0.2.0)
**Source documents:** SDD.md, DESIGN.md, PLAN.md, TEST_PLAN.md, ADR-001–012
**Execution model:** Antigravity Ask mode, one task at a time, gated on confirmed green.

---

## Milestone Map

| Milestone | Phase | Tasks | Gate |
|---|---|---|---|
| M-0 | Environment Setup | 001–007 | Stan smoke test passes |
| M-1 | wf_state + Display Contract | 008–017 | UT-1, UT-2, UT-5 all pass |
| M-2 | Diagnostic Infrastructure | 018–027 | UT-3 (smoke), UT-4 pass |
| M-3 | Phase Scripts + Off-Ramps + Exit | 028–038 | UT-6, UT-3 (exit log) pass; SCENARIO-3 assertions 3.1–3.6 pass |
| M-4 | Quarto Phase Templates | 039–050 | All six templates run interactively without error in learn + practice modes |
| M-5 | Report Template + Integration | 051–056 | SCENARIO-1 all 10 assertions pass |
| M-6 | Acceptance Testing + Polish | 057–064 | SC-1, SC-2, SC-3 verified; v0.1.0 tag applied |
| M-7a **[v0.2.0]** | install.R + export_context patch | 058–062 | All four functions parse; export_context writes to analysis subfolder |
| M-7b **[v0.2.0]** | Unit tests UT-5 ext + UT-7–9 | 063–066 | testthat reports 10 tests, 0 failures |
| M-7c **[v0.2.0]** | User-guide documents | 067–071 | All five docs render as valid Markdown |
| M-7d **[v0.2.0]** | Verification + v0.2.0 tag | 072–074 | SC-4, SC-5 verified; v0.2.0 tag applied |

---

## Milestone M-0: Environment Setup

---

## Task 001: Create project folder structure
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** —

### What to do
Create the canonical RBayesflow project directory tree as specified in DESIGN.md §1. All directories must exist; placeholder `.gitkeep` files go in empty leaf directories. Do not create any R source files yet.

### Files touched
- `R/` — create empty directory
- `templates/` — create empty directory
- `examples/glmm_gaussian/` — create empty directory with `.gitkeep`
- `examples/bernoulli_rare_event/` — create empty directory with `.gitkeep`
- `tests/testthat/` — create empty directory with `.gitkeep`
- `tests/scenarios/` — create empty directory with `.gitkeep`

### Acceptance Criterion
Running `list.dirs(recursive = TRUE)` from the project root returns all six directories listed above. No R source files exist yet.

### On Failure
`TASK 001 FAILED — directory missing: [name of missing directory]`

---

## Task 002: Write DESCRIPTION dependency manifest
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 001

### What to do
Create the `DESCRIPTION` file in the project root. This is not a package DESCRIPTION — it is a dependency manifest used by `renv` to discover required packages. List every required package from SDD.md §9 plus testing utilities. Use the `Imports:` field only; do not include `Package:`, `Version:`, or other package-only fields. Format is one package per line.

### Files touched
- `DESCRIPTION` — create with Imports field listing: brms, cmdstanr, bayesplot, posterior, tidybayes, loo, ggplot2, esquisse, jsonlite, yaml, digest, testthat, withr, lme4

### Acceptance Criterion
`readLines("DESCRIPTION")` contains an `Imports:` line and every package listed above appears in the file. `desc::desc_get_deps()` parses it without error (install `desc` temporarily if needed to verify, then uninstall).

### On Failure
`TASK 002 FAILED — DESCRIPTION missing or malformed: [parse error or missing package name]`

---

## Task 003: Initialise renv and snapshot dependencies
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 002

### What to do
Run `renv::init()` to initialise the renv lockfile, then install all packages listed in `DESCRIPTION`, then run `renv::snapshot()` to write `renv.lock`. If any package fails to install, report the error verbatim and stop — do not proceed.

### Files touched
- `renv.lock` — created by `renv::snapshot()`
- `renv/` — created by `renv::init()`
- `.Rprofile` — created or modified by `renv::init()`

### Acceptance Criterion
`renv::status()` reports "No issues found." and `renv.lock` exists in the project root and is valid JSON parseable by `jsonlite::fromJSON("renv.lock")`.

### On Failure
`TASK 003 FAILED — renv::status() issue or renv.lock missing: [error text]`

---

## Task 004: Install and verify CmdStan
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 003

### What to do
Run `cmdstanr::install_cmdstan()` if CmdStan is not already installed at the version required (≥ 2.33). Then run `cmdstanr::check_cmdstan_toolchain()` to verify the C++ toolchain is operational. Report the installed CmdStan version.

### Files touched
- No project files changed. CmdStan installs to `cmdstanr::cmdstan_path()`.

### Acceptance Criterion
`cmdstanr::check_cmdstan_toolchain(fix = FALSE)` runs without error and `cmdstanr::cmdstan_version()` returns a version string ≥ "2.33.0".

### On Failure
`TASK 004 FAILED — toolchain check failed or version too low: [error text] [version found]`

---

## Task 005: Write and run Stan smoke test
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 004

### What to do
Create `examples/stan_demo.R`. The script must: (1) load `cmdstanr`; (2) write an inline Stan model string for a simple normal-normal conjugate model (8 observations, known sigma, estimate mu); (3) compile it with `cmdstan_model()`; (4) sample with `model$sample(data = list(...), chains = 2, iter_warmup = 200, iter_sampling = 200, refresh = 0)`; (5) print `fit$summary()`. The script must run to completion without error.

### Files touched
- `examples/stan_demo.R` — create

### Acceptance Criterion
`source("examples/stan_demo.R")` completes without error and prints a summary table containing a row for `mu` with `rhat` < 1.05.

### On Failure
`TASK 005 FAILED — Stan smoke test error or rhat ≥ 1.05: [error text or rhat value]`

---

## Task 006: Write known-good example model script
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 005

### What to do
Create `examples/glmm_gaussian/fit_sleepstudy.R`. The script must fit the sleepstudy GLMM defined in TEST_PLAN.md §2.1 using `brm()` with `backend = "cmdstanr"` and `seed = 42`. Save the fitted object to `examples/glmm_gaussian/fit_good.rds` using `saveRDS()`. Do not run the script yet — it will be run manually as part of M-2 testing.

### Files touched
- `examples/glmm_gaussian/fit_sleepstudy.R` — create

### Acceptance Criterion
`readLines("examples/glmm_gaussian/fit_sleepstudy.R")` contains the exact `brm()` call from TEST_PLAN.md §2.1 (formula, data, family, prior, backend, seed fields all present). File is valid R syntax: `parse(file = "examples/glmm_gaussian/fit_sleepstudy.R")` returns without error.

### On Failure
`TASK 006 FAILED — script missing brm() fields or parse error: [missing field or parse error]`

---

## Task 007: Write known-bad example model script
**Status:** [x] Done
**Milestone:** M-0
**Depends on:** Task 006

### What to do
Create `examples/bernoulli_rare_event/fit_rare_event.R`. The script must contain the data-simulation code and `brm()` call defined in TEST_PLAN.md §2.2 (`set.seed(99)`, n=150, 3% event rate, diffuse prior). Save the fitted object to `examples/bernoulli_rare_event/fit_bad.rds`. Do not run the script yet.

### Files touched
- `examples/bernoulli_rare_event/fit_rare_event.R` — create

### Acceptance Criterion
`parse(file = "examples/bernoulli_rare_event/fit_rare_event.R")` returns without error and `readLines(...)` contains `set.seed(99)`, `n <- 150`, `family = bernoulli()`, and `prior(normal(0, 10), class = Intercept)`.

### On Failure
`TASK 007 FAILED — script missing required elements or parse error: [element or error]`

---

## Milestone M-1: wf_state S3 Object and Display Contract

---

## Task 008: Write new_wf_state() constructor
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 007

### What to do
Create `R/wf_state.R`. Implement `new_wf_state(mode, stage)` returning an S3 object of class `c("wf_state", "list")` with the full schema from DESIGN.md §2. Validate that `mode` is one of `"learn"`, `"practice"`, `"expert"` and `stage` is one of `"explore"`, `"confirm"`, `"exit"` — stop with an informative error if either is invalid. Set all other fields to their `NULL` / `NA` / `FALSE` defaults as specified in the schema.

### Files touched
- `R/wf_state.R` — create with `new_wf_state()` only (other methods follow in later tasks)

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_wf_state.R")` passes UT-1 (all four `expect_*` calls green): invalid mode errors, invalid stage errors, valid construction returns `wf_state` class, `acknowledged` defaults to `FALSE`.

### On Failure
`TASK 008 FAILED — UT-1 failure: [failed expect_ call and message]`

---

## Task 009: Write UT-1 test file
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 008

### What to do
Create `tests/testthat/test_wf_state.R` containing exactly the UT-1 test block from TEST_PLAN.md §4. Do not add any other tests yet.

### Files touched
- `tests/testthat/test_wf_state.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_wf_state.R")` reports 1 test, 4 expectations, 0 failures.

### On Failure
`TASK 009 FAILED — UT-1 failure count [N]: [failed expectation and message]`

---

## Task 010: Write print.wf_state() display contract
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 009

### What to do
Add `print.wf_state()` to `R/wf_state.R`. Implement all four display-contract branches from DESIGN.md §3 exactly: (1) diagnostics not yet run — show header only; (2) diagnostics passed + learn mode — call `plot_prior_posterior_overlay(wf)` (stub returning invisible NULL for now) + one-line health summary; (3) diagnostics passed + practice/expert — coefficient table + health summary; (4) diagnostics failed + not acknowledged — failure message + `wf$diagnose()` prompt, no coefficients; (5) diagnostics failed + acknowledged + learn — failure message + acknowledgment note; (6) diagnostics failed + acknowledged + practice/expert — coefficient table + caveat warning. Add stub helpers `cat_wf_header()`, `cat_health_summary_one_line()`, `cat_coefficient_table()`, `cat_diagnostic_failure_message()`, and `plot_prior_posterior_overlay()` as no-ops that print a one-line placeholder — they will be fully implemented in later tasks.

### Files touched
- `R/wf_state.R` — add `print.wf_state()` and stub helpers

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_display_contract.R")` passes UT-2 (both `expect_*` calls green): failed+unacknowledged suppresses "Estimate"; acknowledged+practice shows "Estimate".

### On Failure
`TASK 010 FAILED — UT-2 failure: [failed expect_ call and output captured]`

---

## Task 011: Write UT-2 test file
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 010

### What to do
Create `tests/testthat/test_display_contract.R` containing the UT-2 test block from TEST_PLAN.md §4. The `mock_brmsfit()` helper must return a minimal named list with class `"brmsfit"` and a `fixef` method that returns a one-row matrix with column name "Estimate".

### Files touched
- `tests/testthat/test_display_contract.R` — create
- `tests/testthat/helpers.R` — create with `mock_brmsfit()` helper

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_display_contract.R")` reports 2 tests, 2 expectations, 0 failures.

### On Failure
`TASK 011 FAILED — UT-2 failure count [N]: [failed expectation and message]`

---

## Task 012: Write summary.wf_state() and format.wf_state()
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 011

### What to do
Add `summary.wf_state()` and `format.wf_state()` to `R/wf_state.R`. `summary.wf_state()` must print: mode, stage, fit timestamp (or "not yet fitted"), diagnostic pass/fail status, and audit trail length. `format.wf_state()` must return a single character string suitable for pasting into the Posit Assistant context (one-line summary: mode, stage, diagnostic status).

### Files touched
- `R/wf_state.R` — add `summary.wf_state()` and `format.wf_state()`

### Acceptance Criterion
`wf <- new_wf_state("learn", "explore"); out <- capture.output(summary(wf)); any(grepl("mode", out, ignore.case = TRUE))` returns `TRUE` and `is.character(format(wf))` returns `TRUE`.

### On Failure
`TASK 012 FAILED — summary or format method missing output field: [missing field]`

---

## Task 013: Write export_context() and wf_context.json schema
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 012

### What to do
Create `R/context.R`. Implement `export_context(wf, path = "wf_context.json")` that serialises `wf` to JSON using `jsonlite::toJSON(wf, auto_unbox = TRUE, pretty = TRUE)` and writes to `path`. The JSON must include at minimum: `mode`, `stage`, `rbayesflow_version`, `diagnostics$passed`, `diagnostics$acknowledged`, `audit_trail` length as an integer field `audit_trail_length`.

### Files touched
- `R/context.R` — create with `export_context()`

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_context.R")` passes UT-5: JSON file written, `fromJSON()` parses it, `mode == "learn"`, `stage == "explore"`, `rbayesflow_version` not null.

### On Failure
`TASK 013 FAILED — UT-5 failure: [failed expect_ call and observed value]`

---

## Task 014: Write UT-5 test file
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 013

### What to do
Create `tests/testthat/test_context.R` containing the UT-5 test block from TEST_PLAN.md §4.

### Files touched
- `tests/testthat/test_context.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_context.R")` reports 1 test, 3 expectations, 0 failures.

### On Failure
`TASK 014 FAILED — UT-5 failure count [N]: [failed expectation and message]`

---

## Task 015: Write init_workflow()
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 014

### What to do
Create `R/init.R`. Implement `init_workflow(mode = "learn", stage = "explore")` that: (1) calls `new_wf_state(mode, stage)`; (2) calls `export_context(wf)` to write the initial `wf_context.json`; (3) prints a one-line startup message showing mode, stage, and the path of the written context file; (4) returns `wf` invisibly. Per ADR-011, `mode = "expert"` must be accepted without error.

### Files touched
- `R/init.R` — create with `init_workflow()`

### Acceptance Criterion
`wf <- init_workflow(mode = "learn", stage = "explore")` returns an object of class `"wf_state"` and `file.exists("wf_context.json")` is `TRUE` after the call.

### On Failure
`TASK 015 FAILED — init_workflow() error or wf_context.json not written: [error text]`

---

## Task 016: Write record_fit() hash linkage helper
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 015

### What to do
Add `record_fit(wf, fit)` to `R/wf_state.R`. Per ADR-009, this function: (1) computes a SHA-256 hash of `list(formula = as.character(formula(fit)), data_hash = digest::digest(fit$data), timestamp = Sys.time())` using `digest::digest(..., algo = "sha256")`; (2) stores the hash in `wf$fit_hash`; (3) stores `Sys.time()` in `wf$fit_timestamp`; (4) calls `export_context(wf)` to update `wf_context.json`; (5) returns the updated `wf` object. Also add `check_fit_hash(wf, fit)` that recomputes the hash and stops with an informative error if it does not match `wf$fit_hash`.

### Files touched
- `R/wf_state.R` — add `record_fit()` and `check_fit_hash()`

### Acceptance Criterion
With a mock brmsfit (use `mock_brmsfit()` from `tests/testthat/helpers.R`): `wf2 <- record_fit(wf, mock_brmsfit()); !is.null(wf2$fit_hash) && nchar(wf2$fit_hash) == 64` returns `TRUE`.

### On Failure
`TASK 016 FAILED — record_fit() did not set fit_hash or hash wrong length: [observed value]`

---

## Task 017: Run full M-1 test suite
**Status:** [x] Done
**Milestone:** M-1
**Depends on:** Task 016

### What to do
Run the complete testthat suite for M-1: `testthat::test_dir("tests/testthat/")`. All tests written in Tasks 009, 011, 014 must pass. Fix any failures before marking this task done. Do not add new tests or new source code — this is a verification-only task.

### Files touched
- No files changed. Tests run against existing code.

### Acceptance Criterion
`testthat::test_dir("tests/testthat/")` reports 4 tests, 9 expectations, 0 failures, 0 errors, 0 warnings.

### On Failure
`TASK 017 FAILED — test suite failure: [N failures, failed test names and messages]`

---

## Milestone M-2: Diagnostic Infrastructure

---

## Task 018: Write DIAGNOSTIC_REGISTRY global and family_key() dispatcher
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 017

### What to do
Create `R/diagnostic_registry.R`. Define `DIAGNOSTIC_REGISTRY` as an empty named R list. Implement `family_key(fit)` that extracts the family name from a `brmsfit` object (`family(fit)$family`) and returns the matching key: "bernoulli", "poisson", "negbinomial", "gaussian_hierarchical" (if formula contains a random-effects term), "time_series" (if data has a time index column — check for columns named `time`, `date`, `year`), or "unknown". When key is "unknown", emit a warning: `"No family-specific diagnostic entry for family [name]. Using generic checks only."`.

### Files touched
- `R/diagnostic_registry.R` — create with `DIAGNOSTIC_REGISTRY` and `family_key()`

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_registry.R")` passes UT-4 (both tests green): `DIAGNOSTIC_REGISTRY[["bernoulli"]]` is a function after the bernoulli entry is added in Task 019; unknown family emits warning matching "No family-specific".

### On Failure
`TASK 018 FAILED — family_key() or warning test failure: [failed expect_ and message]`

---

## Task 019: Write Bernoulli registry entry
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 018

### What to do
Add the Bernoulli entry to `DIAGNOSTIC_REGISTRY` in `R/diagnostic_registry.R`, implementing the exact function body from DESIGN.md §8.1. The function takes `(fit, wf_state)` and returns a named list with `checks`, `warnings`, and `plots` fields. The rare-event threshold is 5% (obs_p < 0.05).

### Files touched
- `R/diagnostic_registry.R` — add `DIAGNOSTIC_REGISTRY[["bernoulli"]]`

### Acceptance Criterion
`is.function(DIAGNOSTIC_REGISTRY[["bernoulli"]])` returns `TRUE` after sourcing the file. The function body contains the string `"Rare event detected"`.

### On Failure
`TASK 019 FAILED — bernoulli entry missing or body wrong: [observed]`

---

## Task 020: Write Poisson/NegBinomial registry entries
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 019

### What to do
Add `DIAGNOSTIC_REGISTRY[["poisson"]]` and `DIAGNOSTIC_REGISTRY[["negbinomial"]]` to `R/diagnostic_registry.R`. Both implement the overdispersion check from DESIGN.md §8.2: compute variance-to-mean ratio from 100 rows of posterior predictive draws. If ratio > 2: warn overdispersion. If ratio < 0.5: warn underdispersion. Return the same `list(checks, warnings, plots)` structure.

### Files touched
- `R/diagnostic_registry.R` — add poisson and negbinomial entries

### Acceptance Criterion
`is.function(DIAGNOSTIC_REGISTRY[["poisson"]]) && is.function(DIAGNOSTIC_REGISTRY[["negbinomial"]])` returns `TRUE`. Both function bodies contain the string `"overdispersion"`.

### On Failure
`TASK 020 FAILED — poisson or negbinomial entry missing or body wrong: [observed]`

---

## Task 021: Write hierarchical and time-series registry entries
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 020

### What to do
Add `DIAGNOSTIC_REGISTRY[["gaussian_hierarchical"]]` (DESIGN.md §8.3: warns when n_groups < 5, recommends `refit_noncentered(wf)`) and `DIAGNOSTIC_REGISTRY[["time_series"]]` (DESIGN.md §8.4: checks `acf(residuals(fit))[2]` for `|value| > 0.3`, warns and recommends AR(1) structure) and `DIAGNOSTIC_REGISTRY[["unknown"]]` (returns empty checks/warnings/plots with no-op behaviour).

### Files touched
- `R/diagnostic_registry.R` — add gaussian_hierarchical, time_series, unknown entries

### Acceptance Criterion
All five keys exist as functions: `all(c("bernoulli","poisson","negbinomial","gaussian_hierarchical","time_series","unknown") %in% names(DIAGNOSTIC_REGISTRY))` returns `TRUE`.

### On Failure
`TASK 021 FAILED — missing registry keys: [list of missing keys]`

---

## Task 022: Write UT-4 test file
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 021

### What to do
Create `tests/testthat/test_registry.R` containing the UT-4 test block from TEST_PLAN.md §4. The smoke test for bernoulli must source the registry and check `is.function(DIAGNOSTIC_REGISTRY[["bernoulli"]])`. The unknown-family test must call `family_key(mock_fit(family = "custom_family"))` — add `mock_fit(family)` to `tests/testthat/helpers.R` returning a minimal list with `family(fit)$family` equal to the supplied string.

### Files touched
- `tests/testthat/test_registry.R` — create
- `tests/testthat/helpers.R` — add `mock_fit()` helper

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_registry.R")` reports 2 tests, 2 expectations, 0 failures.

### On Failure
`TASK 022 FAILED — UT-4 failure: [failed expect_ and message]`

---

## Task 023: Write run_diagnostics() generic checks
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 022

### What to do
Create `R/diagnostics.R`. Implement `run_diagnostics(fit, wf)` that: (1) calls `check_fit_hash(wf, fit)` to verify provenance; (2) uses `posterior::summarise_draws(fit)` to extract `rhat` and `ess_bulk`/`ess_tail` per parameter; (3) checks: `rhat_max < 1.01`, `bulk_ess_min > 400`, `tail_ess_min > 400`, `n_divergences == 0` (from `fit$diagnostics()$num_divergent`), `bfmi > 0.3` on all chains, `max_treedepth_hit == FALSE`; (4) stores all raw values and a `failed_criteria` character vector in `wf$diagnostics`; (5) sets `wf$diagnostics$passed` to `TRUE` iff all criteria pass; (6) dispatches to `DIAGNOSTIC_REGISTRY[[family_key(fit)]](fit, wf)` and appends family-specific warnings to `wf$diagnostics$family_checks`; (7) appends to `wf$audit_trail`; (8) calls `export_context(wf)`; (9) returns updated `wf`.

### Files touched
- `R/diagnostics.R` — create with `run_diagnostics()`

### Acceptance Criterion
`is.function(run_diagnostics)` is `TRUE` and `parse(file = "R/diagnostics.R")` returns without error.

### On Failure
`TASK 023 FAILED — parse error or function not defined: [error text]`

---

## Task 024: Write detect_parameterization() helper
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 023

### What to do
Add `detect_parameterization(fit)` to `R/diagnostics.R`. Inspect `brms::stancode(fit)` for the presence of `z_` prefix variables (non-centered) or `r_` prefix variables (centered) as described in DESIGN.md §6. Return `"centered"` or `"non-centered"`. If neither pattern is found (non-hierarchical model), return `NA_character_`.

### Files touched
- `R/diagnostics.R` — add `detect_parameterization()`

### Acceptance Criterion
`is.function(detect_parameterization)` is `TRUE`. Function body contains both `"z_"` and `"r_"` as pattern strings.

### On Failure
`TASK 024 FAILED — function missing or pattern strings absent: [observed]`

---

## Task 025: Write refit_noncentered() stub
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 024

### What to do
Add `refit_noncentered(wf, fit)` to `R/diagnostics.R`. For v1.0 this function: (1) checks that `wf$parameterization == "centered"` and stops with an informative message if not; (2) rebuilds the brms formula replacing `(1 | group)` terms with `(0 + Intercept | group)` syntax using string substitution on `deparse(formula(fit))`; (3) refits using `brm()` with the modified formula, same data, family, prior, and backend; (4) returns a list `list(fit_new = ..., wf_new = record_fit(wf, fit_new))`. Document with a comment that the formula substitution handles only the simple `(1 | group)` case; complex random-effects structures require manual formula revision.

### Files touched
- `R/diagnostics.R` — add `refit_noncentered()`

### Acceptance Criterion
`is.function(refit_noncentered)` is `TRUE` and the function body contains `"0 + Intercept"`.

### On Failure
`TASK 025 FAILED — function missing or non-centered syntax absent: [observed]`

---

## Task 026: Write diagnose.wf_state() method
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 025

### What to do
Add `diagnose.wf_state(wf, ...)` to `R/wf_state.R`. Implement exactly the method body from DESIGN.md §3.1: (1) add rendering guard per ADR-008 — `if (isTRUE(getOption("knitr.in.progress"))) stop("wf$diagnose() must be called interactively...")`; (2) call `show_diagnostic_plots(wf)` (stub: `invisible(NULL)` for now); (3) call `cat_failed_criteria_detail(wf)` (stub: print `wf$diagnostics$failed_criteria`); (4) for `mode != "expert"`: call `readline()` for acknowledgment and set `wf$diagnostics$acknowledged <- TRUE` on "yes"; (5) for `mode == "expert"`: skip readline, set `acknowledged <- TRUE` automatically per ADR-011; (6) append to `wf$audit_trail`; (7) return `wf` invisibly.

### Files touched
- `R/wf_state.R` — add `diagnose.wf_state()`

### Acceptance Criterion
`is.function(diagnose.wf_state)` is `TRUE`. Function body contains `"knitr.in.progress"` and `"readline"` and `"audit_trail"`.

### On Failure
`TASK 026 FAILED — function missing or required strings absent: [missing string]`

---

## Task 027: Run full M-2 test suite
**Status:** [x] Done
**Milestone:** M-2
**Depends on:** Task 026

### What to do
Run `testthat::test_dir("tests/testthat/")`. All tests from M-1 (Tasks 009, 011, 014) and M-2 (Task 022) must pass. Fix any regressions before proceeding. Do not add new tests or source code.

### Files touched
- No files changed.

### Acceptance Criterion
`testthat::test_dir("tests/testthat/")` reports 5 tests (UT-1 through UT-5 minus UT-3 and UT-6 which arrive in M-3), 11 expectations total, 0 failures.

### On Failure
`TASK 027 FAILED — regression in test suite: [N failures, test names and messages]`

---

## Milestone M-3: Phase Scripts, Off-Ramps, and Exit Workflow

---

## Task 028: Write assess_offramps() decision matrix
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 027

### What to do
Create `R/offramps.R`. Implement `assess_offramps(data, outcome_var, outcome_type, goal, n = nrow(data))` per DESIGN.md §4 and ADR-010. Compute `event_rate <- mean(as.integer(factor(data[[outcome_var]])) - 1L)` for binary outcomes with a warning if the variable is not already 0/1-coded. Implement the full decision matrix from DESIGN.md §4.1. Return a named list: `list(n = n, event_rate = event_rate, outcome_type = outcome_type, goal = goal, alternatives = list(...))` where each alternative is `list(method = "...", rationale = "...")`. Present all alternatives with equal weight — no `recommended` field.

### Files touched
- `R/offramps.R` — create with `assess_offramps()`

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_offramps.R")` passes UT-6: `length(result$alternatives) >= 2` and `"full_stan"` in the method vector.

### On Failure
`TASK 028 FAILED — UT-6 failure: [failed expect_ and observed]`

---

## Task 029: Write UT-6 test file
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 028

### What to do
Create `tests/testthat/test_offramps.R` containing the UT-6 test block from TEST_PLAN.md §4. Use a minimal synthetic data frame: `data.frame(y = c(rep(0, 145), rep(1, 5)), x = rnorm(150))`.

### Files touched
- `tests/testthat/test_offramps.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_offramps.R")` reports 1 test, 2 expectations, 0 failures.

### On Failure
`TASK 029 FAILED — UT-6 failure: [failed expect_ and message]`

---

## Task 030: Write exit_workflow() and YAML exit log
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 029

### What to do
Create `R/exit_workflow.R`. Implement `exit_workflow(wf, method, alternatives, path = "exit.yaml")` per DESIGN.md §7. The function must: (1) prompt the user with `readline("Confirm exit to [method]? (yes/no): ")`; (2) stop with a message if not confirmed; (3) construct the YAML exit log from `wf_state` fields — `sample_size` from `wf$n_observations`, `event_rate` from `wf$event_rate`, `declared_goal` from `wf$declared_goal`, `diagnostics_run` from `!is.na(wf$diagnostics$passed)` — not from user-typed strings; (4) write with `yaml::write_yaml()`; (5) set `wf$stage <- "exit"`; (6) set `wf$exit_log` to the constructed list; (7) append to audit trail; (8) return `wf` invisibly.

### Files touched
- `R/exit_workflow.R` — create with `exit_workflow()`

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_exit_workflow.R")` passes UT-3: YAML file written, all required fields present, `user_confirmed == TRUE`, `"full_stan"` in alternatives.

### On Failure
`TASK 030 FAILED — UT-3 failure: [failed expect_ and observed]`

---

## Task 031: Write UT-3 test file
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 030

### What to do
Create `tests/testthat/test_exit_workflow.R` containing the UT-3 test block from TEST_PLAN.md §4. Use `withr::with_tempdir()` to isolate file writes. Mock `readline()` to return `"yes"` using `mockery::mock("yes")` or `withr::local_mocked_bindings(readline = function(...) "yes")`.

### Files touched
- `tests/testthat/test_exit_workflow.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_exit_workflow.R")` reports 1 test, 5 expectations, 0 failures.

### On Failure
`TASK 031 FAILED — UT-3 failure: [N failures, failed expect_ and message]`

---

## Task 032: Write Phase 1 R script — data inspection and off-ramps
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 031

### What to do
Create `R/phase1_exploration.R`. The script must define a function `run_phase1(wf, data, outcome_var, outcome_type, goal)` that: (1) produces an exploratory `ggplot2` plot appropriate to `outcome_type` (histogram for continuous, bar for binary/count); (2) calls `assess_offramps()`; (3) presents the alternatives with equal visual weight using `cat()` formatted output, numbered list, no "recommended" label; (4) prompts user for choice via `readline()`; (5) logs the choice to `wf$audit_trail` with `action = "offramp_choice"` or `"bayesian_selected"`; (6) stores `n_observations`, `event_rate`, and `declared_goal` on `wf`; (7) calls `export_context(wf)`; (8) returns updated `wf`.

### Files touched
- `R/phase1_exploration.R` — create

### Acceptance Criterion
`parse(file = "R/phase1_exploration.R")` returns without error. `is.function(run_phase1)` is `TRUE` after sourcing.

### On Failure
`TASK 032 FAILED — parse error or function not defined: [error text]`

---

## Task 033: Write Phase 2 R script — prior specification and prior predictive
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 032

### What to do
Create `R/phase2_priors.R`. Define `run_phase2(wf, formula, family, priors, data)` that: (1) validates `priors` as a list of `brms::prior()` objects; (2) stores `priors` in `wf$priors_objects` and a human-readable description in `wf$priors_text` via `deparse(priors)`; (3) runs prior predictive simulation: `brm(formula, data, family, prior = priors, sample_prior = "only", backend = "cmdstanr", seed = 42, refresh = 0)`; (4) stores draws in `wf$prior_pred_draws`; (5) plots prior predictive distribution with `bayesplot::ppc_dens_overlay()`; (6) appends to audit trail; (7) calls `export_context(wf)`; (8) returns updated `wf`.

### Files touched
- `R/phase2_priors.R` — create

### Acceptance Criterion
`parse(file = "R/phase2_priors.R")` returns without error. `is.function(run_phase2)` is `TRUE` after sourcing.

### On Failure
`TASK 033 FAILED — parse error or function not defined: [error text]`

---

## Task 034: Write Phase 3 R script — model fitting
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 033

### What to do
Create `R/phase3_fit.R`. Define `run_phase3(wf, formula, data, family, priors, seed = 42, chains = 4, iter_warmup = 1000, iter_sampling = 1000)` that: (1) calls `brm()` with `backend = "cmdstanr"` and the supplied arguments; (2) calls `record_fit(wf, fit)` to store hash and timestamp; (3) calls `detect_parameterization(fit)` and stores result in `wf$parameterization`; (4) if parameterization is `"centered"` and the formula contains a random-effects term, logs a note to the audit trail; (5) calls `export_context(wf)`; (6) returns `list(fit = fit, wf = wf)`.

### Files touched
- `R/phase3_fit.R` — create

### Acceptance Criterion
`parse(file = "R/phase3_fit.R")` returns without error. `is.function(run_phase3)` is `TRUE` after sourcing.

### On Failure
`TASK 034 FAILED — parse error or function not defined: [error text]`

---

## Task 035: Write display.R — mode-aware plot selection
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 034

### What to do
Create `R/display.R`. Implement the full helper functions referenced in `print.wf_state()`: `cat_wf_header(wf)`, `cat_health_summary_one_line(wf)`, `cat_coefficient_table(wf, fit)`, `cat_diagnostic_failure_message(wf)`, `cat_failed_criteria_detail(wf)`, `plot_prior_posterior_overlay(wf, fit)` (calls `bayesplot::mcmc_areas()` overlaid with prior predictive draws using `ggplot2`), and `show_diagnostic_plots(wf, fit)` (calls `bayesplot::mcmc_trace()`, `bayesplot::mcmc_rhat()`, `bayesplot::mcmc_neff()`). Replace all stubs in `R/wf_state.R` with calls to these functions (import from `display.R` via `source()`).

### Files touched
- `R/display.R` — create with all display helpers
- `R/wf_state.R` — replace stub calls with `source("R/display.R")` at top and real function calls

### Acceptance Criterion
`parse(file = "R/display.R")` returns without error. After sourcing both files: `is.function(plot_prior_posterior_overlay) && is.function(show_diagnostic_plots)` returns `TRUE`.

### On Failure
`TASK 035 FAILED — parse error or function not defined: [error text]`

---

## Task 036: Write scenario test script for SCENARIO-3 (exit log)
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 035

### What to do
Create `tests/scenarios/scenario3_exit_log.R`. This script sources all R files, constructs a minimal `wf` object with `n_observations = 150`, `event_rate = 0.03`, `declared_goal = "coefficient estimation"`, mocks `readline()` to return `"yes"`, runs `exit_workflow(wf, method = "firth_logistic", alternatives = c("logistic_bootstrap", "full_stan"))`, reads the YAML file, and asserts all seven assertions from SCENARIO-3 in TEST_PLAN.md §3 using `stopifnot()`. Each failed assertion must print a message identifying which assertion failed.

### Files touched
- `tests/scenarios/scenario3_exit_log.R` — create

### Acceptance Criterion
`source("tests/scenarios/scenario3_exit_log.R")` completes without error and prints "SCENARIO-3: ALL ASSERTIONS PASSED" to console.

### On Failure
`TASK 036 FAILED — SCENARIO-3 assertion failed: [assertion number and observed value]`

---

## Task 037: Write master source script (source_all.R)
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 036

### What to do
Create `R/source_all.R`. This script sources all R files in dependency order: `wf_state.R`, `context.R`, `init.R`, `diagnostic_registry.R`, `diagnostics.R`, `offramps.R`, `exit_workflow.R`, `display.R`, `phase1_exploration.R`, `phase2_priors.R`, `phase3_fit.R`. Users run `source("R/source_all.R")` at the top of each phase template to load the full workflow. Note: `install.R` will be appended to this list in Task 059 (v0.2.0).

### Files touched
- `R/source_all.R` — create

### Acceptance Criterion
`source("R/source_all.R")` completes without error and `exists("init_workflow") && exists("run_diagnostics") && exists("assess_offramps") && exists("exit_workflow")` returns `TRUE`.

### On Failure
`TASK 037 FAILED — source_all.R error or function missing: [error text or missing function name]`

---

## Task 038: Run full M-3 test suite
**Status:** [x] Done
**Milestone:** M-3
**Depends on:** Task 037

### What to do
Run `testthat::test_dir("tests/testthat/")` and then `source("tests/scenarios/scenario3_exit_log.R")`. All six unit tests (UT-1 through UT-6) and SCENARIO-3 must pass.

### Files touched
- No files changed.

### Acceptance Criterion
`testthat::test_dir("tests/testthat/")` reports 6 tests, 0 failures AND `source("tests/scenarios/scenario3_exit_log.R")` prints "SCENARIO-3: ALL ASSERTIONS PASSED".

### On Failure
`TASK 038 FAILED — [unit test or scenario failure]: [N failures, names and messages]`

---

## Milestone M-4: Quarto Phase Templates

---

## Task 039: Write phase1_exploration.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 038

### What to do
Create `templates/phase1_exploration.qmd`. The template must: (1) begin with a Posit Assistant context block as specified in DESIGN.md §5.2 (phase 1); (2) have a setup chunk that sources `R/source_all.R` and calls `init_workflow(mode = params$mode, stage = params$stage)`; (3) include a data-loading chunk where the user supplies their data; (4) include an exploratory plot chunk calling `run_phase1()`; (5) display the off-ramp alternatives as a formatted table using `knitr::kable()`; (6) include a choice-logging chunk with a clear `eval=FALSE`-labelled console instruction for `readline()` use. YAML frontmatter must include `params: {mode: "learn", stage: "explore"}`.

### Files touched
- `templates/phase1_exploration.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase1_exploration.qmd")` returns without error (validates YAML and chunk structure).

### On Failure
`TASK 039 FAILED — quarto_inspect error: [error text]`

---

## Task 040: Write phase2_priors.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 039

### What to do
Create `templates/phase2_priors.qmd`. The template must: (1) Posit Assistant context block for phase 2; (2) setup chunk sourcing `R/source_all.R` and loading `wf` from the previous phase (user supplies `wf` object); (3) prior specification chunk showing example `brms::prior()` usage with comments; (4) prior predictive simulation chunk calling `run_phase2()`; (5) prior predictive plot display chunk.

### Files touched
- `templates/phase2_priors.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase2_priors.qmd")` returns without error.

### On Failure
`TASK 040 FAILED — quarto_inspect error: [error text]`

---

## Task 041: Write phase3_fit.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 040

### What to do
Create `templates/phase3_fit.qmd`. The template must: (1) Posit Assistant context block for phase 3; (2) setup and `wf` loading; (3) `brm()` call chunk via `run_phase3()`; (4) parameterization note chunk displaying `wf$parameterization` with a conditional message if centered and hierarchical; (5) `saveRDS(fit, ...)` chunk for the user to save the fit object.

### Files touched
- `templates/phase3_fit.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase3_fit.qmd")` returns without error.

### On Failure
`TASK 041 FAILED — quarto_inspect error: [error text]`

---

## Task 042: Write phase4_diagnostics.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 041

### What to do
Create `templates/phase4_diagnostics.qmd`. The template must: (1) Posit Assistant context block for phase 4 with explicit instruction not to display coefficients if `wf$diagnostics$passed != TRUE`; (2) `run_diagnostics(fit, wf)` chunk; (3) `print(wf)` chunk showing mode-appropriate output; (4) a chunk marked `eval=FALSE` with explicit user instruction: "Run `wf <- wf$diagnose()` in your R console if diagnostics failed, then continue"; (5) rendering guard note in a callout block explaining why `wf$diagnose()` must be run at the console.

### Files touched
- `templates/phase4_diagnostics.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase4_diagnostics.qmd")` returns without error and `readLines(...)` contains the string `"knitr.in.progress"` (the guard is referenced or documented in the template).

### On Failure
`TASK 042 FAILED — quarto_inspect error or guard string missing: [observed]`

---

## Task 043: Write phase5_ppc.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 042

### What to do
Create `templates/phase5_ppc.qmd`. The template must: (1) Posit Assistant context block for phase 5; (2) `pp_check()` chunk using `brms::pp_check(fit)`; (3) mode-aware display — in `mode = "learn"`: print a callout with a Socratic question ("What systematic differences do you see between the blue and dark lines?"); in `mode = "practice"`: standard PPC table; (4) update `wf$ppc_complete <- TRUE` and `export_context(wf)`.

### Files touched
- `templates/phase5_ppc.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase5_ppc.qmd")` returns without error.

### On Failure
`TASK 043 FAILED — quarto_inspect error: [error text]`

---

## Task 044: Write phase6_loo.qmd template
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 043

### What to do
Create `templates/phase6_loo.qmd`. The template must: (1) Posit Assistant context block for phase 6; (2) `loo()` call chunk — `loo_fit <- loo(fit)`; (3) if multiple models, `loo_compare()` chunk; (4) LOO table display with `knitr::kable()`; (5) update `wf$loo_complete <- TRUE`, `wf$loo_table <- as.data.frame(loo_fit$estimates)`, and `export_context(wf)`.

### Files touched
- `templates/phase6_loo.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/phase6_loo.qmd")` returns without error.

### On Failure
`TASK 044 FAILED — quarto_inspect error: [error text]`

---

## Task 045: Run SCENARIO-4 (off-ramp equal weighting)
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 044

### What to do
Create `tests/scenarios/scenario4_offramps.R`. Run `assess_offramps()` on `dat_bad` (binary, n=150, event_rate=0.03, goal="coefficient estimation") and verify the four assertions from TEST_PLAN.md §3 SCENARIO-4 using `stopifnot()`.

### Files touched
- `tests/scenarios/scenario4_offramps.R` — create

### Acceptance Criterion
`source("tests/scenarios/scenario4_offramps.R")` completes without error and prints "SCENARIO-4: ALL ASSERTIONS PASSED".

### On Failure
`TASK 045 FAILED — SCENARIO-4 assertion failed: [assertion number and observed]`

---

## Task 046: Run M-4 template inspection check
**Status:** [x] Done
**Milestone:** M-4
**Depends on:** Task 045

### What to do
Create `tests/scenarios/check_templates.R`. The script calls `quarto::quarto_inspect()` on all six phase templates and the report template path (even though the report template doesn't exist yet — skip it with a `tryCatch`). Print PASS/FAIL for each.

### Files touched
- `tests/scenarios/check_templates.R` — create

### Acceptance Criterion
`source("tests/scenarios/check_templates.R")` prints "PASS" for all six phase templates (phase1 through phase6) with zero errors.

### On Failure
`TASK 046 FAILED — template inspect failure: [template name and error]`

---

## Milestone M-5: Report Template and End-to-End Integration

---

## Task 047: Fit known-good model for integration testing
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 046

### What to do
Run `source("examples/glmm_gaussian/fit_sleepstudy.R")` to produce `fit_good.rds`. This task requires Stan to be operational (verified in M-0). The fit must complete successfully with `seed = 42`. This task may take 2–5 minutes.

### Files touched
- `examples/glmm_gaussian/fit_good.rds` — created by running the script

### Acceptance Criterion
`file.exists("examples/glmm_gaussian/fit_good.rds")` is `TRUE` and `fit_good <- readRDS("examples/glmm_gaussian/fit_good.rds"); inherits(fit_good, "brmsfit")` returns `TRUE`.

### On Failure
`TASK 047 FAILED — fit_good.rds missing or not brmsfit class: [error text]`

---

## Task 048: Fit known-bad model for integration testing
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 047

### What to do
Run `source("examples/bernoulli_rare_event/fit_rare_event.R")` to produce `fit_bad.rds`.

### Files touched
- `examples/bernoulli_rare_event/fit_bad.rds` — created by running the script

### Acceptance Criterion
`file.exists("examples/bernoulli_rare_event/fit_bad.rds")` is `TRUE` and `inherits(readRDS("examples/bernoulli_rare_event/fit_bad.rds"), "brmsfit")` returns `TRUE`.

### On Failure
`TASK 048 FAILED — fit_bad.rds missing or not brmsfit class: [error text]`

---

## Task 049: Write bayesflow_report.qmd template
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 048

### What to do
Create `templates/bayesflow_report.qmd`. This is the only template intended to be rendered. It must: (1) use YAML frontmatter `params: {mode: "learn", stage: "explore", wf_path: "wf.rds", fit_path: "fit.rds"}`; (2) load `wf <- readRDS(params$wf_path)` and `fit <- readRDS(params$fit_path)` in the setup chunk — no refitting; (3) render sections for each completed phase by reading from `wf_state` fields only; (4) include diagnostic summary from `wf$diagnostics`; (5) include coefficient table only if `wf$diagnostics$passed == TRUE` OR `wf$diagnostics$acknowledged == TRUE`; (6) include PPC plots from `wf$ppc_summary` if `wf$ppc_complete`; (7) include LOO table from `wf$loo_table` if `wf$loo_complete`; (8) include full audit trail as a collapsible section; (9) include exit log if `wf$stage == "exit"`.

### Files touched
- `templates/bayesflow_report.qmd` — create

### Acceptance Criterion
`quarto::quarto_inspect("templates/bayesflow_report.qmd")` returns without error.

### On Failure
`TASK 049 FAILED — quarto_inspect error: [error text]`

---

## Task 050: Write SCENARIO-1 integration test script
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 049

### What to do
Create `tests/scenarios/scenario1_learn_mode.R`. The script must execute all ten steps and assert all ten assertions from TEST_PLAN.md §3 SCENARIO-1 using `stopifnot()`. Steps 1.1–1.9 run programmatically (no interactive readline); mock readline to always return `"yes"` for choice prompts. Step 1.10 (render Phase 7 report) calls `quarto::quarto_render("templates/bayesflow_report.qmd", output_file = "test_report.html", execute_params = list(wf_path = "wf.rds", fit_path = "examples/glmm_gaussian/fit_good.rds", mode = "learn", stage = "explore"))` and asserts the output file exists.

### Files touched
- `tests/scenarios/scenario1_learn_mode.R` — create

### Acceptance Criterion
`source("tests/scenarios/scenario1_learn_mode.R")` completes without error and prints "SCENARIO-1: ALL 10 ASSERTIONS PASSED".

### On Failure
`TASK 050 FAILED — SCENARIO-1 assertion failed: [assertion number, step description, observed value]`

---

## Task 051: Write SCENARIO-2 integration test script
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 050

### What to do
Create `tests/scenarios/scenario2_diagnostic_gate.R`. The script executes all seven assertions from TEST_PLAN.md §3 SCENARIO-2 for both `mode = "learn"` and `mode = "practice"`, using `fit_bad.rds` and mocked readline.

### Files touched
- `tests/scenarios/scenario2_diagnostic_gate.R` — create

### Acceptance Criterion
`source("tests/scenarios/scenario2_diagnostic_gate.R")` completes without error and prints "SCENARIO-2: ALL 7 ASSERTIONS PASSED (learn + practice)".

### On Failure
`TASK 051 FAILED — SCENARIO-2 assertion failed: [assertion number, mode, observed value]`

---

## Task 052: Write SCENARIO-5 parameterization test script
**Status:** [x] Done
**Milestone:** M-5
**Depends on:** Task 051

### What to do
Create `tests/scenarios/scenario5_parameterization.R`. Execute the three assertions from TEST_PLAN.md §3 SCENARIO-5 using `fit_good.rds` (centered by default in brms for the sleepstudy random-effects structure).

### Files touched
- `tests/scenarios/scenario5_parameterization.R` — create

### Acceptance Criterion
`source("tests/scenarios/scenario5_parameterization.R")` completes without error and prints "SCENARIO-5: ALL 3 ASSERTIONS PASSED".

### On Failure
`TASK 052 FAILED — SCENARIO-5 assertion failed: [assertion number and observed]`

---

## Milestone M-6: Acceptance Testing and Documentation Polish

---

## Task 053: Run all unit tests (final)
**Status:** [x] Done
**Milestone:** M-6
**Depends on:** Task 052

### What to do
Run `testthat::test_dir("tests/testthat/")`. All six unit tests (UT-1 through UT-6) must pass with zero failures, zero errors, zero warnings.

### Files touched
- No files changed.

### Acceptance Criterion
`testthat::test_dir("tests/testthat/")` reports 6 tests, 0 failures, 0 errors, 0 warnings.

### On Failure
`TASK 053 FAILED — unit test failure: [N failures, test names and messages]`

---

## Task 054: Run all scenario tests (final)
**Status:** [x] Done
**Milestone:** M-6
**Depends on:** Task 053

### What to do
Source all five scenario scripts in order: scenario1 through scenario5. All must print their "ALL ASSERTIONS PASSED" message.

### Files touched
- No files changed.

### Acceptance Criterion
All five scenario scripts complete without error and each prints its "ALL ASSERTIONS PASSED" line.

### On Failure
`TASK 054 FAILED — scenario failure: [scenario name, assertion number, observed]`

---

## Task 055: Verify README Getting Started section
**Status:** [x] Done
**Milestone:** M-6
**Depends on:** Task 054

### What to do
Read `README.md` and verify it contains: (1) installation instructions for all dependencies including `cmdstanr::install_cmdstan()`; (2) instructions to clone/copy the project folder; (3) instructions to run `source("R/source_all.R")`; (4) instructions to call `init_workflow(mode = "learn")`; (5) a link to the Gelman et al. (2020) reference; (6) a description of each mode ("learn", "practice", "expert") in plain language. Update any sections that are missing or incomplete.

### Files touched
- `README.md` — update if missing any of the six required sections

### Acceptance Criterion
`readLines("README.md")` contains all six strings: "install_cmdstan", "source_all.R", "init_workflow", "mode = \"learn\"", "Gelman", "expert".

### On Failure
`TASK 055 FAILED — README missing required strings: [list of missing strings]`

---

## Task 056: Update CHANGELOG.md for v0.1.0
**Status:** [x] Done
**Milestone:** M-6
**Depends on:** Task 055

### What to do
Add a `## [0.1.0] — 2026-09-08` section to `CHANGELOG.md` listing: all seven phases implemented, all three success criteria verified, all eleven ADRs accepted.

### Files touched
- `CHANGELOG.md` — add v0.1.0 entry

### Acceptance Criterion
`readLines("CHANGELOG.md")[1:5]` contains `"0.1.0"` and `grepl("0\\.1\\.0", readLines("CHANGELOG.md"))` has at least one `TRUE`.

### On Failure
`TASK 056 FAILED — CHANGELOG missing v0.1.0 entry: [observed first 5 lines]`

---

## Task 057: Apply v0.1.0 git tag
**Status:** [x] Done
**Milestone:** M-6
**Depends on:** Task 056

### What to do
Run `git add -A && git commit -m "feat: RBayesflow v0.1.0 — all SC-1, SC-2, SC-3 verified" && git tag -a v0.1.0 -m "RBayesflow v0.1.0"`.

### Files touched
- Git history and tags only.

### Acceptance Criterion
`system("git tag")` output contains `"v0.1.0"` and `system("git log --oneline -1")` contains `"v0.1.0"`.

### On Failure
`TASK 057 FAILED — git tag missing or commit failed: [error text]`

---

## Milestone M-7a: install.R Functions + export_context() Patch **[v0.2.0]**

---

## Task 058: Create data/ directory and rbf_analysis_path() helper
**Status:** [ ] Pending
**Milestone:** M-7a
**Depends on:** Task 057

### What to do
Create the `data/` directory at the project root with a `.gitkeep` placeholder (analyses that go in it are user-generated and will not be committed by default — add `data/*/` to `.gitignore` but keep `data/.gitkeep`). Then create `R/install.R` and implement `rbf_analysis_path()` as the first and only function in the file at this stage. Per DESIGN.md §10.4: return `getwd()` when called from inside a `data/<name>/` folder (detect by: `.Rprofile` exists in `getwd()` AND `getwd()` is an immediate child of a `data/` folder under the project root as found by `rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))`). Otherwise return the project root path string.

### Files touched
- `data/.gitkeep` — create
- `.gitignore` — add `data/*/` entry (keep `!data/.gitkeep`)
- `R/install.R` — create with `rbf_analysis_path()` only

### Acceptance Criterion
`parse(file = "R/install.R")` returns without error. `is.function(rbf_analysis_path)` is `TRUE` after `source("R/install.R")`. Function body contains `"rprojroot"` and `"DESCRIPTION"`.

### On Failure
`TASK 058 FAILED — parse error, function missing, or required strings absent: [observed]`

---

## Task 059: Add install.R to source_all.R load order
**Status:** [ ] Pending
**Milestone:** M-7a
**Depends on:** Task 058

### What to do
Append `source(file.path(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")), "R", "install.R"))` as the last line of `R/source_all.R`. This ensures `guide()`, `rbf_new()`, and `rbf_analysis_path()` are available in every phase template session. Do not move this line above the core file sources — `install.R` depends on `wf_state.R` being loaded first (because `guide()` calls `new_wf_state()` indirectly via phase-detection field access).

### Files touched
- `R/source_all.R` — append one `source()` line at the end

### Acceptance Criterion
`tail(readLines("R/source_all.R"), 3)` contains `"install.R"`. `source("R/source_all.R")` completes without error and `exists("rbf_analysis_path")` returns `TRUE`.

### On Failure
`TASK 059 FAILED — install.R not in source_all.R or source error: [error text]`

---

## Task 060: Implement rbf_new(name)
**Status:** [ ] Pending
**Milestone:** M-7a
**Depends on:** Task 059

### What to do
Add `rbf_new(name)` to `R/install.R` per DESIGN.md §10.2. The function must: (1) validate `name` — stop with a clear message if `name` is empty, contains whitespace, or contains path separators (`/` or `\`); (2) determine the project root via `rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))`; (3) construct `target <- file.path(root, "data", name)`; (4) stop with `"Analysis '<name>' already exists at <target>"` if the directory exists; (5) create the directory with `dir.create(target, recursive = TRUE)`; (6) write `.Rprofile` containing exactly `source(file.path("..", "..", "R", "source_all.R"))`; (7) write `wf_context.json` containing `{}`; (8) write `README.md` containing `# Analysis: <name>\nCreated: <Sys.Date()>\n`; (9) print `"Analysis '<name>' created at data/<name>/. Open that folder as your working directory, then call init_workflow()."`.

### Files touched
- `R/install.R` — add `rbf_new()`

### Acceptance Criterion
`parse(file = "R/install.R")` returns without error. `is.function(rbf_new)` is `TRUE` after sourcing. The function body contains `"already exists"`, `".Rprofile"`, and `"source_all.R"`.

### On Failure
`TASK 060 FAILED — parse error, function missing, or required strings absent: [observed]`

---

## Task 061: Implement guide(wf)
**Status:** [ ] Pending
**Milestone:** M-7a
**Depends on:** Task 060

### What to do
Add `guide(wf)` to `R/install.R` per DESIGN.md §10.3. Implement all ten phase-detection branches in the exact precedence order specified. Each branch must print exactly three lines: `"Current phase : <name>"`, `"What to expect: <one sentence>"`, `"Next step     : <exact call or filename>"`. Return `wf` invisibly. Do not call `cat()` for anything else — no blank lines, no headers, no colour codes.

### Files touched
- `R/install.R` — add `guide()`

### Acceptance Criterion
`parse(file = "R/install.R")` returns without error. `is.function(guide)` is `TRUE` after sourcing. For a fresh `wf <- new_wf_state("learn", "explore")`: `out <- capture.output(guide(wf)); length(out) == 3 && grepl("Phase 1", out[1]) && grepl("phase1_exploration", out[3])` returns `TRUE`.

### On Failure
`TASK 061 FAILED — parse error, wrong output lines, or wrong phase detected: [observed output]`

---

## Task 062: Implement rbf_install(dry_run = FALSE) and patch export_context()
**Status:** [ ] Pending
**Milestone:** M-7a
**Depends on:** Task 061

### What to do
This task has two parts — implement them in this order:

**Part A — rbf_install():** Add `rbf_install(dry_run = FALSE)` to `R/install.R` per DESIGN.md §10.1. Implement all seven steps. Each step prints `"✓ Step N: <description>"` on pass or `"✗ Step N: <description>\n  → <inline remediation>"` on fail. Stop after Step 2 if the toolchain is broken. When `dry_run = TRUE`, perform checks but skip all `install*()` calls. Return invisibly a named logical vector of length 7. Print the final summary line `"RBayesflow environment: N/7 checks passed."`.

**Part B — export_context() path patch:** Edit `R/context.R` to change the default value of the `path` argument from `"wf_context.json"` to `file.path(rbf_analysis_path(), "wf_context.json")`. This is a one-line change. Do not alter any other logic in the function.

### Files touched
- `R/install.R` — add `rbf_install()`
- `R/context.R` — change `path` argument default value

### Acceptance Criterion
`parse(file = "R/install.R")` and `parse(file = "R/context.R")` both return without error. `is.function(rbf_install)` is `TRUE`. `rbf_install` function body contains `"dry_run"` and `"7 checks passed"`. `readLines("R/context.R")` contains `"rbf_analysis_path()"` in the function signature line.

### On Failure
`TASK 062 FAILED — parse error or required strings absent: [file and missing string]`

---

## Milestone M-7b: Unit Tests UT-5 Extension + UT-7, UT-8, UT-9 **[v0.2.0]**

---

## Task 063: Extend UT-5 for analysis subfolder path
**Status:** [ ] Pending
**Milestone:** M-7b
**Depends on:** Task 062

### What to do
Add the subfolder-path test case from TEST_PLAN.md §4 UT-5 extension to `tests/testthat/test_context.R`. The new test must use `withr::with_tempdir()` to create a mock project structure with a `DESCRIPTION` file at root and a `data/my_analysis/.Rprofile` file, then call `withr::with_dir(file.path("data", "my_analysis"), { export_context(wf) })` with no explicit `path` argument, and assert that `"wf_context.json"` exists in `data/my_analysis/`.

### Files touched
- `tests/testthat/test_context.R` — add one new `test_that()` block

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_context.R")` reports 2 tests, 0 failures. The new test name contains "subfolder".

### On Failure
`TASK 063 FAILED — subfolder path test failure: [observed path or error]`

---

## Task 064: Write UT-7 test file — rbf_install() step vector
**Status:** [ ] Pending
**Milestone:** M-7b
**Depends on:** Task 063

### What to do
Create `tests/testthat/test_install_rbf_install.R` containing the UT-7 test block from TEST_PLAN.md §4. Use `with_mocked_bindings()` to stub `cmdstanr::check_cmdstan_toolchain` and `cmdstanr::cmdstan_version`. Call `rbf_install(dry_run = TRUE)` — the `dry_run` flag prevents any actual install calls, so mocking is only needed for the check calls in Steps 1, 2, 5.

### Files touched
- `tests/testthat/test_install_rbf_install.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_install_rbf_install.R")` reports 1 test, 3 expectations (`type == "logical"`, `length == 7`, `!is.null(names(result))`), 0 failures.

### On Failure
`TASK 064 FAILED — UT-7 failure: [failed expect_ and observed]`

---

## Task 065: Write UT-8 test file — rbf_new() folder structure
**Status:** [ ] Pending
**Milestone:** M-7b
**Depends on:** Task 064

### What to do
Create `tests/testthat/test_install_rbf_new.R` containing both UT-8 test blocks from TEST_PLAN.md §4: (1) correct subfolder structure test — runs `rbf_new("my_analysis")` inside a tempdir that has a `data/` folder and a `DESCRIPTION` file, checks for `.Rprofile`, `wf_context.json`, and `README.md`; (2) duplicate error test — creates `data/existing/` first, then calls `rbf_new("existing")` and expects an error matching `"already exists"`.

### Files touched
- `tests/testthat/test_install_rbf_new.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_install_rbf_new.R")` reports 2 tests, 5 expectations, 0 failures.

### On Failure
`TASK 065 FAILED — UT-8 failure: [test name, failed expect_, observed]`

---

## Task 066: Write UT-9 test file — guide() phase detection
**Status:** [ ] Pending
**Milestone:** M-7b
**Depends on:** Task 065

### What to do
Create `tests/testthat/test_install_guide.R` containing the three UT-9 test blocks from TEST_PLAN.md §4: (1) fresh `wf_state` → Phase 1 detected, `phase1_exploration` in output; (2) failed+unacknowledged diagnostics → `"diagnose"` in output; (3) `stage == "exit"` → `"exit"` in output (case-insensitive).

### Files touched
- `tests/testthat/test_install_guide.R` — create

### Acceptance Criterion
`testthat::test_file("tests/testthat/test_install_guide.R")` reports 3 tests, 5 expectations, 0 failures.

### On Failure
`TASK 066 FAILED — UT-9 failure: [test name, failed expect_, observed output]`

---

## Milestone M-7c: User-Guide Documents **[v0.2.0]**

---

## Task 067: Write docs/user-guide/00-overview.md
**Status:** [ ] Pending
**Milestone:** M-7c
**Depends on:** Task 066

### What to do
Create `docs/user-guide/00-overview.md`. The document must cover: (1) what RBayesflow is and is not (sequencer, not a package; no new statistical methods); (2) the three-tier help system (written guide → `guide()` function → Posit Assistant); (3) the three user modes in plain language (a paragraph each — no table); (4) the seven workflow phases as a numbered list with a one-sentence description of each; (5) a "Quick start" section listing the three commands a working scientist runs to start a new analysis: `source("R/install.R"); rbf_install(); rbf_new("my_analysis")`. Write in plain English; no R jargon without inline definition.

### Files touched
- `docs/user-guide/00-overview.md` — create

### Acceptance Criterion
`file.exists("docs/user-guide/00-overview.md")` is `TRUE`. `readLines(...)` contains all five strings: `"rbf_install"`, `"rbf_new"`, `"guide"`, `"learn"`, `"practice"`. File is valid Markdown: no unclosed fenced code blocks, no broken header syntax.

### On Failure
`TASK 067 FAILED — file missing or required strings absent: [missing strings]`

---

## Task 068: Write docs/user-guide/01-installation.md
**Status:** [ ] Pending
**Milestone:** M-7c
**Depends on:** Task 067

### What to do
Create `docs/user-guide/01-installation.md`. The document must cover these steps in order, each as a numbered section: (1) Prerequisites — R ≥ 4.3 (link to r-project.org) and RTools 4.5 for Windows (link to cran.r-project.org/bin/windows/Rtools/); (2) Get RBayesflow — `git clone` command or download-as-zip instructions; (3) The one manual step — install cmdstanr from r-universe: exact two-line R code block (`install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))`); (4) Run `rbf_install()` — exact code block, explanation of each of the 7 printed steps, and what to do if a step fails (pointer back to inline remediation message); (5) Verify — `source("R/source_all.R"); wf <- init_workflow()` and expected output. Each code block must be fenced with ` ```r `.

### Files touched
- `docs/user-guide/01-installation.md` — create

### Acceptance Criterion
`file.exists("docs/user-guide/01-installation.md")` is `TRUE`. `readLines(...)` contains all four strings: `"r-universe"`, `"rbf_install"`, `"RTools"`, `"init_workflow"`. File contains at least three fenced ` ```r ` code blocks.

### On Failure
`TASK 068 FAILED — file missing or required strings/code blocks absent: [missing item]`

---

## Task 069: Write docs/user-guide/02-posit-assistant-setup.md
**Status:** [ ] Pending
**Milestone:** M-7c
**Depends on:** Task 068

### What to do
Create `docs/user-guide/02-posit-assistant-setup.md`. Structure: (1) introduction — what Posit Assistant does in RBayesflow (reads `wf_context.json`, does not run code, works in both RStudio and Positron); (2) "Option A — RStudio" section: numbered steps: install/enable Posit Assistant if not present → Gear icon → Providers → Add Provider → OpenRouter → paste API key → select model → Done; (3) "Option B — Positron" section: numbered steps: Command Palette → "Positron Assistant: Configure Language Model Providers" → Add → Custom Provider → Base URL `https://openrouter.ai/api/v1` → API key → model → Done; (4) "Getting an OpenRouter API key" sub-section: openrouter.ai → Sign up → Dashboard → Keys → Create Key → copy; (5) "Which model to use" sub-section: recommend `meta-llama/llama-3.3-70b-instruct:free` or `google/gemini-2.0-flash-exp:free` as free-tier options with adequate context windows; note that the user can substitute any model with context window ≥ 8 K tokens; (6) "Using Posit Assistant with RBayesflow" section: run `export_context(wf)` to refresh context → paste the standard prompt into the Assistant chat; (7) "What Posit Assistant will and won't do" — bulleted will/won't list as specified in the Cowork brief.

### Files touched
- `docs/user-guide/02-posit-assistant-setup.md` — create

### Acceptance Criterion
`file.exists("docs/user-guide/02-posit-assistant-setup.md")` is `TRUE`. `readLines(...)` contains all five strings: `"openrouter.ai"`, `"RStudio"`, `"Positron"`, `"export_context"`, `"wf_context.json"`.

### On Failure
`TASK 069 FAILED — file missing or required strings absent: [missing strings]`

---

## Task 070: Write docs/user-guide/03-starting-an-analysis.md and 04-workflow-phases.md
**Status:** [ ] Pending
**Milestone:** M-7c
**Depends on:** Task 069

### What to do
Create two documents in one task (both are short):

**03-starting-an-analysis.md** — cover: `rbf_new("name")` usage, opening the analysis subfolder as the working directory in RStudio/Positron, calling `init_workflow(mode, stage)` with a table of mode choices and when to use each, calling `guide(wf)` at any time to find out where you are, the `.Rprofile` note (what it does; what to do if you have a conflicting global `.Rprofile`).

**04-workflow-phases.md** — one section per phase (seven sections). Each section: phase name and number as heading, one sentence on purpose, the exact function or template to run, what to expect to see on success, when to move to the next phase. For Phase 4 only: add a callout block explaining that `wf$diagnose()` must be run at the console, not rendered.

### Files touched
- `docs/user-guide/03-starting-an-analysis.md` — create
- `docs/user-guide/04-workflow-phases.md` — create

### Acceptance Criterion
Both files exist. `readLines("docs/user-guide/03-starting-an-analysis.md")` contains `"rbf_new"`, `"init_workflow"`, `"guide"`, `".Rprofile"`. `readLines("docs/user-guide/04-workflow-phases.md")` contains `"Phase 1"`, `"Phase 4"`, `"diagnose"`, `"Phase 7"`.

### On Failure
`TASK 070 FAILED — file missing or required strings absent: [file name and missing strings]`

---

## Task 071: Write docs/user-guide/05-plotting-reference.md
**Status:** [ ] Pending
**Milestone:** M-7c
**Depends on:** Task 070

### What to do
Create `docs/user-guide/05-plotting-reference.md`. Cover every plot listed in the Cowork brief's plotting reference specification, in the order specified. For each plot: (1) a level-3 heading with the plot name; (2) a **Purpose** line (one sentence); (3) a **Call** fenced code block with the exact R call; (4) a **What to look for** paragraph: good result, bad result, action on bad result; (5) a **Learn-mode note** sentence where behaviour differs between modes. Sections:

- Phase 1: `esquisse::esquisser(data)` — drag-and-drop EDA; note that the generated ggplot2 code can be pasted into the template
- Phase 2: `bayesplot::ppc_dens_overlay(y, prior_pred_draws)` — prior predictive density
- Phase 4: `bayesplot::mcmc_trace(fit)`, `bayesplot::mcmc_rhat(brms::rhat(fit))`, `bayesplot::mcmc_neff(brms::neff_ratio(fit))`, `bayesplot::mcmc_pairs(fit, ...)`
- Phase 5: `bayesplot::ppc_dens_overlay(y, posterior_predict(fit))`, `bayesplot::ppc_stat(y, posterior_predict(fit), stat = "mean")`, `tidybayes::add_epred_draws(data, fit) |> ggplot(...)`
- Phase 6: `loo::loo_compare(loo1, loo2)` — reading ELPD difference
- Phase 7: note only (plots reproduced from saved objects; no new plots generated)

### Files touched
- `docs/user-guide/05-plotting-reference.md` — create

### Acceptance Criterion
`file.exists("docs/user-guide/05-plotting-reference.md")` is `TRUE`. `readLines(...)` contains all six strings: `"esquisse"`, `"mcmc_trace"`, `"ppc_dens_overlay"`, `"add_epred_draws"`, `"loo_compare"`, `"Learn-mode"`. File contains at least nine fenced ` ```r ` code blocks (one per distinct plot call).

### On Failure
`TASK 071 FAILED — file missing, required strings absent, or code block count wrong: [missing item or count found]`

---

## Milestone M-7d: Verification and v0.2.0 Tag **[v0.2.0]**

---

## Task 072: Run full test suite including UT-7 through UT-9
**Status:** [ ] Pending
**Milestone:** M-7d
**Depends on:** Task 071

### What to do
Run `testthat::test_dir("tests/testthat/")`. All tests UT-1 through UT-9 (10 test files after the UT-5 extension adds one test) must pass with zero failures, zero errors, zero warnings.

### Files touched
- No files changed.

### Acceptance Criterion
`testthat::test_dir("tests/testthat/")` reports 10 tests, 0 failures, 0 errors, 0 warnings.

### On Failure
`TASK 072 FAILED — test suite failure: [N failures, test file names and messages]`

---

## Task 073: Manual SC-4 and SC-5 verification
**Status:** [ ] Pending
**Milestone:** M-7d
**Depends on:** Task 072

### What to do
Manually verify SC-4 and SC-5 from SDD.md §6. SC-4: confirm all five user-guide files exist and each renders as valid Markdown (open each in RStudio's Preview or run `rmarkdown::render()` with `output_format = "md_document"`). SC-5: read `05-plotting-reference.md` and confirm it has a section for `esquisse::esquisser` and separate sections for at least three `bayesplot` / `tidybayes` calls, each with a "What to look for" paragraph. Record the result as a one-line note appended to `CHANGELOG.md`.

### Files touched
- `CHANGELOG.md` — append SC-4 and SC-5 verification note

### Acceptance Criterion
`readLines("CHANGELOG.md")` contains `"SC-4"` and `"SC-5"` and `"verified"` within the last 10 lines of the file.

### On Failure
`TASK 073 FAILED — CHANGELOG missing SC-4/SC-5 verification lines, or a user-guide file failed to render: [file name and error]`

---

## Task 074: Update CHANGELOG.md and apply v0.2.0 git tag
**Status:** [ ] Pending
**Milestone:** M-7d
**Depends on:** Task 073

### What to do
Add a `## [0.2.0] — 2026-09-10` section to `CHANGELOG.md` listing: (1) new file `R/install.R` with `rbf_install()`, `rbf_new()`, `guide()`, `rbf_analysis_path()`; (2) `export_context()` path default changed to `rbf_analysis_path()`; (3) `docs/user-guide/` with five documents; (4) ADR-012 accepted; (5) SDD, DESIGN, PLAN, TEST_PLAN updated with `[v0.2.0]` sections; (6) SC-4 and SC-5 verified. Then run `git add -A && git commit -m "feat: RBayesflow v0.2.0 — user guidance layer, rbf_install, rbf_new, guide, user-guide docs" && git tag -a v0.2.0 -m "RBayesflow v0.2.0"`.

### Files touched
- `CHANGELOG.md` — add v0.2.0 entry
- Git history and tags only

### Acceptance Criterion
`grepl("0\\.2\\.0", readLines("CHANGELOG.md"))` has at least one `TRUE` AND `system("git tag")` output contains `"v0.2.0"`.

### On Failure
`TASK 074 FAILED — CHANGELOG missing v0.2.0 or git tag failed: [observed]`

---

*End of tasks.md*
*Total tasks: 74 | Milestones: 11 (M-0 through M-7d)*
*v0.1.0 tasks: 001–057 (all complete) | v0.2.0 tasks: 058–074 (pending)*
