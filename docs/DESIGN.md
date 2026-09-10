# RBayesflow — Software Design Document (DESIGN.md)

**Status:** Active
**Last updated:** 2026-09-10
**Owner:** John Peach
**Spec Kit role:** The **HOW** document — implementation architecture, data schemas, function contracts, and behavioral logic. Derived from the SDD and the accepted ADRs. Feeds `PLAN.md` (WHEN) and `tasks.md` (atomic work items).

---

## 1. Architecture Overview

RBayesflow is a **project folder** (not a package) containing:

```
RBayesflow/
├── R/
│   ├── wf_state.R              # S3 constructor, print, summary, format, to_json
│   ├── diagnostic_registry.R   # DIAGNOSTIC_REGISTRY global; family-specific checks
│   ├── diagnostics.R           # run_diagnostics(), generic checks (Rhat, ESS, div, BFMI)
│   ├── offramps.R              # assess_offramps(); non-Bayesian alternative assessment
│   ├── exit_workflow.R         # exit_workflow(); YAML exit log generation
│   ├── display.R               # mode-aware plot selection and coefficient gating
│   ├── context.R               # export_context(); wf_state → JSON for Posit Assistant
│   ├── init.R                  # init_workflow(); sets mode, stage, loads registry
│   └── install.R               # [v0.2.0] rbf_install(), rbf_new(), guide(), rbf_analysis_path()
├── data/                       # [v0.2.0] one subfolder per analysis
│   └── <analysis_name>/
│       ├── wf_context.json
│       └── .Rprofile           # auto-sources source_all.R from project root
├── templates/
│   ├── bayesflow_report.qmd    # Phase 7 Quarto report template (parameterized)
│   ├── phase1_exploration.qmd  # Phase 1 goal declaration + exploratory graphics
│   ├── phase2_priors.qmd       # Phase 2 prior specification + prior predictive
│   ├── phase3_fit.qmd          # Phase 3 model fitting
│   ├── phase4_diagnostics.qmd  # Phase 4 MCMC diagnostics
│   ├── phase5_ppc.qmd          # Phase 5 posterior predictive checks
│   └── phase6_loo.qmd          # Phase 6 model comparison
├── examples/
│   ├── glmm_gaussian/          # Known-good model (used in SC-1 acceptance test)
│   └── bernoulli_rare_event/   # Known-bad diagnostic model (used in SC-2 test)
├── docs/
│   ├── SDD.md
│   ├── DESIGN.md
│   ├── PLAN.md
│   ├── TEST_PLAN.md
│   ├── adr/                    # ADR-001 … ADR-012
│   └── user-guide/             # [v0.2.0] written user guide (Markdown)
│       ├── 00-overview.md
│       ├── 01-installation.md
│       ├── 02-posit-assistant-setup.md
│       ├── 03-starting-an-analysis.md
│       ├── 04-workflow-phases.md
│       └── 05-plotting-reference.md
├── wf_context.json             # Legacy top-level context (created only when no analysis subfolder)
├── DESCRIPTION                 # Dependency manifest (not a package; used for renv)
├── renv.lock                   # Pinned dependency versions
├── README.md
└── CHANGELOG.md
```

The workflow has **no central database**, no background process, and no network calls. Everything lives in R objects in the user's session and in files in the project folder.

---

## 2. `wf_state` Schema

Constructed by `new_wf_state(mode, stage, formula, data, family)`. Returns an S3 object of class `c("wf_state", "list")`.

```r
structure(
  list(
    # --- Identity ---
    rbayesflow_version = "0.1.0",
    mode  = "learn",        # "learn" | "practice" | "expert"
    stage = "explore",      # "explore" | "confirm" | "exit"

    # --- Model specification ---
    formula  = NULL,        # brms formula object
    family   = NULL,        # brms family object
    data_hash = NULL,       # SHA-256 of data frame (reproducibility anchor)

    # --- Phase 2: Priors ---
    priors_text    = NULL,  # character vector: human-readable prior descriptions
    priors_objects = NULL,  # list of brms prior objects (brmsprior class)
    prior_pred_draws = NULL, # draws_array from prior predictive simulation

    # --- Phase 3: Fit ---
    fit_timestamp  = NULL,  # POSIXct
    stan_backend   = NULL,  # "cmdstanr" | "rstan"
    parameterization = NULL, # "centered" | "non-centered" (hierarchical only)

    # --- Phase 4: Diagnostics ---
    diagnostics = list(
      passed        = NA,   # logical: all criteria passed
      acknowledged  = FALSE, # logical: user has called wf$diagnose()
      rhat_max      = NA,   # numeric: max Rhat across all parameters
      bulk_ess_min  = NA,   # integer: min bulk ESS
      tail_ess_min  = NA,   # integer: min tail ESS
      n_divergences = NA,   # integer
      bfmi          = NA,   # numeric vector (one per chain)
      max_treedepth_hit = NA, # logical
      family_checks = list(), # output from registry entry (named list)
      failed_criteria = character() # names of failed checks
    ),

    # --- Phase 5: PPC ---
    ppc_complete = FALSE,
    ppc_summary  = NULL,    # list: test statistic results

    # --- Phase 6: LOO ---
    loo_complete = FALSE,
    loo_table    = NULL,    # data.frame from loo_compare()

    # --- Audit trail ---
    audit_trail = list(),   # list of named lists: {phase, action, timestamp, notes}

    # --- Exit log ---
    exit_log = NULL         # named list matching YAML exit log schema (or NULL)
  ),
  class = c("wf_state", "list")
)
```

### 2.1 Validation Rules (enforced by `new_wf_state()`)

- `mode` must be one of `"learn"`, `"practice"`, `"expert"`
- `stage` must be one of `"explore"`, `"confirm"`, `"exit"`
- `diagnostics$acknowledged` may only be set to `TRUE` by `wf$diagnose()`; never set directly

---

## 3. Display Contract

`print.wf_state(wf, ...)` implements the full display contract. Logic:

```
if (is.na(wf$diagnostics$passed)) {
  # Fit not yet run — show mode/stage/formula summary only
  cat_wf_header(wf)

} else if (wf$diagnostics$passed) {
  # Clean diagnostics
  if (wf$mode == "learn") {
    plot_prior_posterior_overlay(wf)
    cat_health_summary_one_line(wf)
  } else {  # practice or expert
    cat_coefficient_table(wf)
    cat_health_summary_one_line(wf)
  }

} else {
  # Failed diagnostics
  if (!wf$diagnostics$acknowledged) {
    cat_diagnostic_failure_message(wf)
    cat("→ Call wf$diagnose() to review the failing checks.\n")
  } else {
    # Acknowledged but failed — show what is safe to show
    if (wf$mode == "learn") {
      cat_diagnostic_failure_message(wf)
      cat("(Diagnostics reviewed and acknowledged.)\n")
    } else {
      cat_coefficient_table(wf)
      cat("⚠ Diagnostics failed (acknowledged). Coefficients shown with caveat.\n")
    }
  }
}
```

### 3.1 `wf$diagnose()` Method

```r
diagnose.wf_state <- function(wf, ...) {
  stopifnot(inherits(wf, "wf_state"))

  # 1. Display bayesplot diagnostic panels (family-appropriate)
  show_diagnostic_plots(wf)

  # 2. Print human-readable failure summary
  cat_failed_criteria_detail(wf)

  # 3. Prompt for acknowledgment
  ack <- readline("Have you reviewed the diagnostic plots? (yes/no): ")
  if (tolower(ack) == "yes") {
    wf$diagnostics$acknowledged <- TRUE
    wf$audit_trail <- append(wf$audit_trail, list(list(
      phase    = 4,
      action   = "diagnostic_acknowledged",
      timestamp = Sys.time(),
      failed_criteria = wf$diagnostics$failed_criteria
    )))
    message("Acknowledged. You may now access coefficient output via print(wf).")
  } else {
    message("Not acknowledged. Call wf$diagnose() again when ready.")
  }

  invisible(wf)
}
```

---

## 4. Non-Bayesian Off-Ramp Logic

Called at the end of Phase 1 (goal declaration) via `assess_offramps(data, outcome_type, goal, n)`.

### 4.1 Decision Matrix

| Outcome type | n | Hierarchical | Goal | Offered alternatives |
|---|---|---|---|---|
| Binary | < 200 | No | Coefficient estimation | Logistic + bootstrap CI; Firth penalized logistic; Full Stan |
| Binary | ≥ 200 | No | Prediction | Logistic + bootstrap CI; Full Stan |
| Count | Any | No | Overdispersion test | Poisson + quasi-likelihood; NegBinomial GLM; Full Stan |
| Continuous | Any | No | Simple comparison | t-test + Cohen's d; LM + bootstrap; Full Stan |
| Any | Any | Yes | Any | Warn: hierarchical + small n → low power; GLMM with REML; Full Stan |

### 4.2 Presentation Rule

Alternatives are presented with **equal visual weight**. No alternative is marked as "recommended" by default. Each alternative receives a one-sentence stakes rationale drawn from the data (sample size, event rate, declared goal). The user chooses; the choice is logged to `wf_state$audit_trail`.

---

## 5. Posit Assistant Integration

### 5.1 `export_context(wf)` function

Serializes `wf_state` to `wf_context.json` using `jsonlite::toJSON(wf_state, auto_unbox = TRUE, pretty = TRUE)`. Called:
- As the last step of each phase script
- Explicitly by the user via `export_context(wf)` at session start

**[v0.2.0]** The `path` default is now `file.path(rbf_analysis_path(), "wf_context.json")` rather than the project root — see §10.

### 5.2 Posit Assistant System Prompt Convention

A block comment at the top of each phase `.qmd` template instructs Posit Assistant:

```
# POSIT ASSISTANT CONTEXT
# This is phase {N} of the RBayesflow Bayesian workflow.
# Current mode: {mode} | Current stage: {stage}
# Read wf_context.json in the current analysis folder (data/<name>/)
# for full workflow state. Do NOT suggest running any code directly.
# Always present suggestions as non-executable previews. Do not display
# coefficient summaries unless wf_state$diagnostics$passed == TRUE.
```

### 5.3 Mode-Specific LLM Behavior

| Mode | LLM behavior |
|---|---|
| `"learn"` | After each diagnostic display, Posit Assistant generates a Socratic question ("What does the Rhat value of 1.89 tell you about this chain?") and a pointer to the relevant Gelman et al. (2020) section |
| `"practice"` | LLM interprets diagnostics on explicit request only (`wf$ask("Why did this diverge?")`) |
| `"expert"` | LLM silent unless directly queried |

### 5.4 Provider setup **[v0.2.0]**

Posit Assistant supports OpenRouter as a named provider in both RStudio (≥ 0.7.7) and Positron. Configuration is done through the IDE settings UI; no R code is involved. Full step-by-step instructions are in `docs/user-guide/02-posit-assistant-setup.md`. The `wf_context.json`-based integration (ADR-006) is unchanged.

---

## 6. Parameterization Transparency

When `brm()` is called with a random-effects formula, `detect_parameterization(fit)` inspects the generated Stan code (`stancode(fit)`) for the presence of `z_` prefix variables (non-centered) or `r_` prefix variables (centered):

- Result is stored in `wf$parameterization` (`"centered"` or `"non-centered"`)
- Always logged to `audit_trail` with a note
- If `n_divergences > 0` AND `parameterization == "centered"` AND the model is hierarchical: the workflow generates a one-call switch suggestion: `refit_noncentered(wf)`, which rebuilds the brms formula with `(0 + Intercept | group)` syntax and refits

---

## 7. Exit Workflow

`exit_workflow(wf, method, justification)` writes the YAML exit log and closes the stage:

```yaml
exit:
  method: <chosen method>           # e.g., "logistic regression + bootstrap CI"
  justification: <drawn from wf_state fields, not user self-report>
  evidence_inspected:
    - sample_size: <n>
    - event_rate: <p>               # only for binary outcomes
    - declared_goal: <string>
    - diagnostics_run: <true/false>
  alternatives_considered:
    - <alternative 1>
    - <alternative 2>
  user_confirmed: true
  timestamp: <ISO 8601>
  rbayesflow_version: "0.1.0"
```

The justification field is populated from `wf_state` objects (sample size, event rate extracted from the data, declared goal from Phase 1 input). It is not a user-typed string.

---

## 8. Diagnostic Registry — Built-in Entries

See ADR-007 for the registry architecture. The `DIAGNOSTIC_REGISTRY` lives in `R/diagnostic_registry.R` as an R list object. Built-in entries for v1.0:

### 8.1 Bernoulli

```r
DIAGNOSTIC_REGISTRY[["bernoulli"]] <- function(fit, wf_state) {
  y      <- fit$data[[as.character(formula(fit)[[2]])]]
  pp     <- posterior_predict(fit)
  obs_p  <- mean(y)
  pred_p <- mean(pp)
  ratio  <- obs_p / pred_p

  warn <- character()
  if (obs_p < 0.05) warn <- c(warn,
    sprintf("Rare event detected (observed rate %.1f%%). Consider Firth penalized logistic or a penalized prior.", obs_p * 100))

  list(
    checks  = list(event_rate = obs_p, predicted_rate = pred_p, calibration_ratio = ratio),
    warnings = warn,
    plots   = list(calibration = bayesplot::ppc_bars(y, pp[sample(nrow(pp), 100), ]))
  )
}
```

### 8.2 Poisson / Negative Binomial

Checks variance-to-mean ratio (from posterior predictive draws). If ratio > 2: warn overdispersion. If ratio < 0.5: warn underdispersion.

### 8.3 Hierarchical (any family, few groups)

Triggered when the model includes a random-effects term AND the number of unique groups is < 5. Warns about centered parameterization and recommends `refit_noncentered(wf)`.

### 8.4 Time-Series

Checks ACF of residuals (`acf(residuals(fit))`) for significant autocorrelation at lag 1. If `|acf[2]| > 0.3`: warn and recommend AR(1) error structure in brms.

---

## 9. Open Design Questions

All ODQ-1 through ODQ-4 raised in the v0.1.0 draft were resolved through ADRs before v0.1.0 tagging. No open design questions remain for v0.2.0; the sole new decision (analysis subfolder layout and `export_context()` path default) is resolved in ADR-012.

---

## 10. `R/install.R` — function contracts **[v0.2.0]**

`R/install.R` is a new file added in v0.2.0. It contains the four functions below and no others. It is sourced by `source_all.R` after the core files so that `guide()` has access to `wf_state` helpers, but `rbf_install()` and `rbf_new()` are designed to be callable from a fresh R session where `source_all.R` has not yet been sourced. Concretely, the file is written so that its top-level definitions have no side effects and reference only base R and the packages named in each function's body (loaded lazily inside the function).

### 10.1 `rbf_install()`

Standalone function; no `source_all.R` dependency. Checks and installs in order, printing `✓` or `✗` for each step:

```
Step 1  R version ≥ 4.3
Step 2  RTools / C++ toolchain  (cmdstanr::check_cmdstan_toolchain())
Step 3  renv installed           (install if missing)
Step 4  renv::restore()          (installs all pinned packages)
Step 5  CmdStan present          (cmdstanr::cmdstan_version())
Step 6  CmdStan install          (cmdstanr::install_cmdstan() if absent)
Step 7  Stan smoke test          (source("R/stan_demo.R"))
```

Signature: `rbf_install(dry_run = FALSE)`. When `dry_run = TRUE`, the function performs each check but skips any installation call; it is used by UT-7.

Returns invisibly: a named logical vector, one entry per step (`TRUE` = passed, `FALSE` = failed). Stops at Step 2 if the toolchain is broken, because subsequent steps cannot succeed. Prints a final one-line summary: `"RBayesflow environment: N/7 checks passed."` For any failing step, full remediation instructions are printed inline, with a pointer to `docs/user-guide/01-installation.md`.

The one prerequisite `rbf_install()` cannot bootstrap on its own is `cmdstanr` itself, which must be installed from r-universe before `rbf_install()` can call `cmdstanr::check_cmdstan_toolchain()`. This one manual step is documented in `docs/user-guide/01-installation.md` and printed as a hint if `requireNamespace("cmdstanr", quietly = TRUE)` returns `FALSE` at entry to Step 2.

### 10.2 `rbf_new(name)`

`name` is a character string (no spaces; use underscores). Creates `data/<name>/` with:
- `.Rprofile` containing `source(file.path("..", "..", "R", "source_all.R"))` so the workflow loads automatically when RStudio opens the analysis subfolder as a working directory
- `wf_context.json` initialised as an empty `{}` placeholder
- `README.md` containing the analysis name and creation date

Prints confirmation: `"Analysis 'name' created at data/name/. Open that folder as your working directory, then call init_workflow()."`

Errors with a clear message if `data/<name>/` already exists. Errors if `name` is empty, contains whitespace, or contains path separators.

### 10.3 `guide(wf)`

Takes a `wf_state` object. Reads `wf$audit_trail` and the diagnostic fields to determine the current phase. Prints three lines:

```
Current phase : <phase name and number>
What to expect: <one sentence describing normal output at this phase>
Next step     : <exact function call or template to open>
```

Phase detection logic (in order of precedence — the first branch that matches wins):

1. `wf$stage == "exit"` → "You have exited the Bayesian path. See `exit.yaml` for the log."
2. `wf$loo_complete` → Phase 7: open `templates/bayesflow_report.qmd`
3. `wf$ppc_complete` → Phase 6: run `wf <- run_phase6(wf, fit)`
4. `isTRUE(wf$diagnostics$passed)` → Phase 5: open `templates/phase5_ppc.qmd`
5. `isFALSE(wf$diagnostics$passed) && !wf$diagnostics$acknowledged` → Phase 4 (blocked): call `wf <- wf$diagnose()`
6. `isFALSE(wf$diagnostics$passed) && wf$diagnostics$acknowledged` → Phase 4 (unblocked): open `templates/phase5_ppc.qmd` with caveat
7. `!is.null(wf$fit_timestamp) && !is.na(wf$fit_timestamp)` → Phase 4: open `templates/phase4_diagnostics.qmd`
8. `!is.null(wf$priors_objects)` → Phase 3: open `templates/phase3_fit.qmd`
9. `!is.null(wf$declared_goal)` → Phase 2: open `templates/phase2_priors.qmd`
10. default → Phase 1: open `templates/phase1_exploration.qmd`

Returns the `wf` object invisibly so `guide()` can be chained without side effects.

### 10.4 `rbf_analysis_path()`

One-line helper. Returns `getwd()` when called from inside a `data/<name>/` folder (detected by the presence of `.Rprofile` in the working directory and the working directory being an immediate child of a `data/` folder under the RBayesflow project root, as determined by `rprojroot::find_root(rprojroot::has_file("DESCRIPTION"))`). Otherwise returns the project root.

### 10.5 `export_context()` path change

The existing `path` argument defaults to `"wf_context.json"` in v0.1.0. In v0.2.0, the default becomes `file.path(rbf_analysis_path(), "wf_context.json")`. This preserves full backward compatibility: callers that pass an explicit `path` are unaffected. `init_workflow()` no longer needs to pass a special path — the new default resolves correctly whether the session is running inside an analysis subfolder or at the project root.

---

*End of DESIGN.md*
