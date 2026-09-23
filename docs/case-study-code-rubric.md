# RBayesflow Case Study Code Rubric

*For R scripts accompanying the Bayesian Workflow case study series.*
*Derived from movies_ch16_analysis.R, sleep_ch17_analysis.R, and*
*clinical_trial_ch18.R, plus the correction history from each session.*

---

Text and data for each chapter analysis should be read from the book [website](https://avehtari.github.io/Bayesian-Workflow/).

## 1. File header

Every script begins with a header comment block. No code appears before it.
The block covers, in this order:

```r
# data/<chapter_folder>/<chapter_slug>_analysis.R
#
# Ch N — "<exact chapter title from the book>"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# <One paragraph: what this script does and what it deliberately omits.>
# Be explicit about scope so a student is not surprised by missing figures.
#
# Data:
#   <source, how acquired, N, key variables>
#
# Models covered:
#   <one line per model: variable name, formula, family, purpose>
#
# Figures produced:
#   <Fig N.M — short description>
#   ...
#
# Book-target posterior summaries (used for success criteria):
#   <parameter: book value, 95% CI>
#   ...
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/<chapter_folder>/.
```

The chapter title must match the book's chapter title exactly, including
capitalisation and punctuation. Model variable names in the header must match
the variable names used in the script body.

---

## 2. Downloading required files

Before writing any analysis code, Claude must download all required files for
the chapter from the book website and save them to the chapter folder at
`C:\Users\johnx\Documents\WildPeaches\Projects\RBayesflow\data\<chapter_folder>\`.

### What to download

Fetch the chapter HTML page from
`https://avehtari.github.io/Bayesian-Workflow/<chapter_slug>/<chapter_slug>.html`
and inspect it for the following asset types:

| Asset type | Where to look | Save as |
|---|---|---|
| Stan model files (`.stan`) | `print_stan_file()` calls in the chapter code, or linked directly | `<model_name>.stan` |
| CSV / RDS data files | `read_csv()`, `read_rds()`, or download links in the chapter | original filename |
| Any other data file | Data-loading calls in the chapter code | original filename |

Download every asset to the chapter folder before the analysis script is
written. Do not reference any file in the script that has not already been
saved locally.

### How to download

Use the Filesystem MCP tool to write Stan files and small data files directly.
For larger data files, include a download step in the `## Data acquisition`
section of the script (see §3 below) that caches the file locally on first run.

### Stan files

For each Stan file, reconstruct the exact model code from the `print_stan_file()`
output shown on the chapter page and write it verbatim to the chapter folder.
Do not paraphrase or simplify Stan code. If the chapter page does not show a
Stan file's full contents, fetch the raw file from the book's GitHub repository
at `https://github.com/avehtari/Bayesian-Workflow/tree/master/<chapter_slug>/`.

### Verification

After downloading, list the chapter folder contents and confirm that every
file referenced by the analysis script is present before delivery.

---

## 3. Environment setup

The first executable lines after the header, in this exact order:

```r
source("../../R/source_all.R")
```

Then package loads for any packages not covered by `source_all.R`, using
`library()` calls. Keep this list as short as possible; prefer `pkg::fn()`
over adding a library call for a package used only once.

Then global options:

```r
options(brms.backend = "cmdstanr", mc.cores = 4)
```

Then a seed constant if the chapter uses simulated data or if any downstream
call needs reproducibility:

```r
SEED <- <value matching the book, or 42 if the book does not specify>
```

Then create output directories and open the results log (see §17):

```r
dir.create("figs", showWarnings = FALSE)
results_con <- file("results.txt", open = "wt")
log_result  <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = results_con, append = FALSE)
}
log_result("=== Ch N results log — ", format(Sys.time()), " ===")
```

The `log_result()` helper writes every message to both the console and
`results.txt` in the analysis subfolder. All book comparison blocks,
diagnostic summaries, LOO tables, and the wrap-up block use `log_result()`
instead of bare `cat()`. See §14, §15, and §17 for usage.

Do not set a theme globally. Apply `ggplot2::theme_minimal()` per figure.
Do not call `theme_set()`.

**No `pause()` function.** The script is run line by line interactively.
Replace any `pause()` calls with a comment: `# <what to look for before
continuing>`.

---

## 4. Namespace discipline

`source_all.R` loads the RBayesflow machinery. It does not guarantee a safe
search-path order for CRAN packages. To avoid masking errors:

- **brms functions:** always qualify: `brms::prior()`, `brms::fixef()`,
  `brms::ranef()`, `brms::coef()`, `brms::VarCorr()`,
  `brms::posterior_predict()`, `brms::posterior_summary()`.
- **bayesplot functions:** always qualify: `bayesplot::ppc_dens_overlay()`,
  `bayesplot::ppc_stat()`, `bayesplot::ppc_loo_pit_overlay()`.
- **loo functions:** always qualify: `loo::loo()`, `loo::loo_compare()`,
  `loo::pareto_k_values()`.
- **ggplot2 functions:** always qualify: `ggplot2::ggplot()`,
  `ggplot2::aes()`, `ggplot2::geom_*()`, `ggplot2::scale_*()`,
  `ggplot2::facet_*()`, `ggplot2::labs()`, `ggplot2::theme_minimal()`, etc.
- **dplyr / tidyr:** qualify any function that shares a name with base R
  (`dplyr::filter()`, `dplyr::select()`, `tidyr::drop_na()`). Others may
  be left unqualified if `library(dplyr)` appears at the top.
- **RBayesflow functions** (`init_workflow`, `run_phase1`, `run_phase2`,
  `run_phase3`, `run_diagnostics`, `assess_offramps`, `export_context`,
  `guide`, `diagnose.wf_state`): do NOT qualify — they are loaded by
  `source_all.R` and have no CRAN conflict.

**Known masking conflicts to watch for:**

| Function | Masked by | Safe call |
|---|---|---|
| `rhat()` | bayesplot masks brms | use `brms::rhat()` or avoid |
| `posterior_predict()` | tidybayes shadows brms | `brms::posterior_predict()` |
| `loo()` | brms re-exports loo, but loo:: is safer | `loo::loo()` |
| `filter()` | dplyr masks stats | `dplyr::filter()` |
| `sd()`, `var()`, `mad()` | posterior masks base | avoid bare; use `base::sd()` or explicit |

---

## 5. Priors

Define all priors using `brms::prior()`. Never use the bare `prior()` form
even though brms re-exports it, because `prior` is easily confused with
R's built-in distribution functions.

Each prior block is preceded by a comment block explaining the reasoning:

```r
# Intercept: normal(250, 100)
#   Puts ~95% prior probability on mean reaction time in (50, 450) ms.
#   Generous bounds — the data determine the exact value.
# Days slope: normal(0, 20)
#   Puts ~95% probability on per-day change in (-40, 40) ms/day.
#   Rules out implausible slopes without forcing a direction.
# sigma: exponential(0.02)   [rate = 0.02, mean = 50 ms]
#   Expected residual scale from domain knowledge about reaction time.
#   Note: exponential() in brms takes the RATE, not the mean.
#         exponential(0.02) has mean 1/0.02 = 50.
priors_m1 <- c(
  brms::prior(normal(250, 100),  class = Intercept),
  brms::prior(normal(0, 20),     class = b),
  brms::prior(exponential(0.02), class = sigma)
)
```

**Parameterisation note — always state it.** Whenever an exponential prior
appears, add an inline comment: `# rate = X, mean = 1/X`. This is the most
common source of student confusion.

---

## 6. Model fitting

Every `brm()` call is wrapped in `run_phase3()` when the RBayesflow wrappers
are appropriate. Use `brm()` directly only when the chapter logic requires
bypassing the wrapper (for example, when fitting multiple models to different
datasets in a single section).

Standard call pattern with wrappers:

```r
result_mN <- run_phase3(
  wf      = wf,
  formula = <formula>,
  family  = <family>(),
  priors  = priors_mN,
  data    = <data>,
  seed    = SEED,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_mN <- result_mN$fit
wf     <- result_mN$wf
saveRDS(fit_mN, "fit_mN.rds")
```

Standard call pattern with direct `brm()`:

```r
fit_mN <- brm(
  formula   = <formula>,
  data      = <data>,
  family    = <family>(),
  prior     = priors_mN,
  save_pars = save_pars(all = TRUE),
  backend   = "cmdstanr",
  seed      = SEED,
  refresh   = 0
)
saveRDS(fit_mN, "fit_mN.rds")
```

`save_pars = save_pars(all = TRUE)` is required whenever LOO-CV will be
run on the fit. `refresh = 0` suppresses per-iteration console output,
which is noisy and unreadable in interactive line-by-line execution.

**`saveRDS()` is mandatory** immediately after every fit. The fit is the
most expensive computation in the script; losing it to a session crash is
the most common source of wasted time.

### Stan-native chapters

Some chapters require Stan models that brms cannot express — custom survival
likelihoods, imputation models, or non-standard parameterisations. When a
chapter falls into this category:

1. State it explicitly at the top of the script header: "This script uses
   `cmdstanr` directly because the chapter's likelihoods cannot be expressed
   through brms formula interfaces."
2. Name every `.stan` file in the header's "Models covered" block.
3. Use a thin `cstan()` helper for compilation and sampling:

```r
cstan <- function(stan_file, data = list(), seed = SEED, chains = 4) {
  model <- cmdstan_model(stan_file)
  model$sample(
    data            = data,
    seed            = seed,
    chains          = chains,
    parallel_chains = chains,
    refresh         = 0
  )
}
```

4. Because `run_phase3()` wraps `brm()`, it cannot be used. Set `wf_state`
   audit fields manually after each fit:

```r
wf$fit_timestamp <- Sys.time()
wf$fit_hash <- digest::digest(
  list(formula = "<description>", data_hash = digest::digest(dat, algo = "sha256")),
  algo = "sha256"
)
```

5. Use `fit$diagnostic_summary()` in place of `run_diagnostics()` and log
   the results explicitly (see §14 for the pattern). Set
   `wf$diagnostics$passed` and `wf$diagnostics$acknowledged` manually after
   inspecting the output.

6. For models with per-observation latent parameters (imputation, per-subject
   varying effects, mixture components), use `iter = 4000, warmup = 2000`
   as the default. Default iterations are insufficient for these models:
   ESS below 400 and Rhat above 1.01 are common at the defaults, and the
   marginal cost of doubling iterations is low relative to the debugging cost
   of a borderline chain. Note the elevated iteration count in the script
   header.

---

## 7. LOO-CV

Do not pass `moment_match = TRUE` to the initial `add_criterion()` call.
Moment matching triggers rstan recompilation on Windows when cmdstanr is
the backend. Use the conditional loop pattern instead:

```r
fit_mN <- add_criterion(fit_mN, "loo", save_psis = TRUE)

# Conditional moment matching for high-k observations
if (any(loo::pareto_k_values(fit_mN$criteria$loo) > 0.7)) {
  fit_mN <- loo::loo_moment_match(fit_mN, loo = fit_mN$criteria$loo)
}
```

**LOO absence must be logged with a reason.** When a chapter does not use
LOO-CV, do not simply set `wf$loo_complete <- FALSE` and move on. Add a
`log_result()` call explaining why:

```r
# LOO-CV not applicable for this chapter.
# Reason: <e.g., survival model with censoring; comparison is qualitative
#          via posterior K-M overlay; no common likelihood for loo_compare()>
log_result("LOO-CV: not performed.")
log_result("Reason: <explanation>")
wf$loo_complete <- FALSE
```

This prevents the absence from looking like a bug during review.

---

## 8. Diagnostics

### Standard diagnostic pattern (brms fits)

After every `run_diagnostics()` / `print(wf)` block:

```r
log_result("--- Diagnostics: fit_mN ---")
log_result("  passed      : ", wf$diagnostics$passed)
log_result("  Rhat_max    : ", round(wf$diagnostics$rhat_max,    4))
log_result("  bulk_ESS_min: ", round(wf$diagnostics$bulk_ess_min, 0))
log_result("  tail_ESS_min: ", round(wf$diagnostics$tail_ess_min, 0))
log_result("  divergences : ", wf$diagnostics$n_divergences)
if (length(wf$diagnostics$failed_criteria) > 0)
  log_result("  FAILED      : ", paste(wf$diagnostics$failed_criteria, collapse = ", "))
```

### Stan-native diagnostic pattern (cmdstanr fits)

When `run_diagnostics()` is not available, use `fit$diagnostic_summary()`
and log results in the same format:

```r
log_result("--- Diagnostics: fit_mN (cmdstanr) ---")
diag <- fit_mN$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag$num_divergent))
log_result("  num_max_treedepth : ", sum(diag$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag$ebfmi, 4), collapse = ", "))
smry <- fit_mN$summary(c("param1", "param2"))   # list key parameters
capture_result(smry, label = "fit_mN: parameter summary")
```

After inspecting the output, set the wf gate fields manually:

```r
wf$diagnostics$passed      <- TRUE   # or FALSE
wf$diagnostics$acknowledged <- TRUE
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "diagnostic_acknowledged",
  timestamp = Sys.time(),
  notes     = "<brief description of what was found>"
)))
```

### Comparing results to book targets

Every parameter that has a book-stated target gets an explicit comparison
block immediately after the diagnostic log. The comparison must include a
note explaining any discrepancy:

```r
log_result("Book target M2 (sim): p[1] ≈ 0.01, p[2] ≈ 0.02")
capture_result(fit2s$summary(c("p[1]", "p[2]")), label = "M2 (sim): posterior summary")
# If the posterior mean differs from the book target, add a note:
log_result("NOTE: posterior mode shifted left of true value.")
log_result("      Expected: prior shrinkage from Beta(1,10) pulling estimates")
log_result("      toward zero. Both colors correctly separated; bias is small.")
```

A discrepancy is not always a bug. Acceptable causes include prior shrinkage,
different random seeds, and updated Stan or brms versions. An unexplained
discrepancy — one where neither shrinkage nor sampling noise accounts for the
gap — should be flagged with `# INVESTIGATE:` and left for the reviewer.

---

## 9. Posterior predictive checks

### PPC method must be justified when not using bayesplot

The default PPC tools (`bayesplot::ppc_dens_overlay()`,
`bayesplot::ppc_loo_pit_overlay()`) are appropriate for most regression
models. They are not appropriate for:

- Survival / time-to-event models (use posterior K-M simulation)
- Models with censoring (same)
- Count models where the outcome is sparse (use `bayesplot::ppc_bars()`)
- Models where the estimand is a contrast or derived quantity

When a chapter's PPC uses a non-standard method, the Phase 5 section must
include a comment block before the PPC code explaining the substitution:

```r
# PPC METHOD NOTE
# Standard bayesplot ppc_dens_overlay() is not appropriate here because
# <reason: e.g., the outcome is a censored survival time and the likelihood
# is geometric, not Gaussian>. We instead simulate posterior K-M curves
# and compare their shape to the empirical K-M. This checks the same
# property — whether the model's implied distribution matches the data —
# using the correct estimand for a survival outcome.
```

---

## 10. Simulated data and generative models

When a chapter builds and tests a generative simulator before fitting real
data, the simulator functions must be:

1. Written before any Stan code is shown in the script.
2. Documented with a comment explaining what the function simulates and which
   model it corresponds to.
3. Verified by running the simulator and plotting the result before calling
   any Stan sampler.

The standard recovery check pattern is:

```r
# Fit model on simulated data to verify recovery
set.seed(SEED)
sim_dat <- sim_fn(n = 1000, true_param = c(...))
fit_sim  <- cstan("model.stan", data = sim_dat)

log_result("Book target (sim): param ≈ X")
capture_result(fit_sim$summary("param"), label = "Simulated data recovery")

# Note whether the posterior mode matches the true value and why any
# offset is expected (prior shrinkage, sampling noise) or unexpected.
log_result("NOTE: posterior mode = X.XX vs true value X.X.")
log_result("      Offset attributable to: <prior shrinkage / sampling noise / other>.")
log_result("      Model correctly <separates groups / recovers direction / etc.>.")
```

Do not claim recovery if the posterior is visibly shifted from the true value.
Instead, explain the offset and confirm that the model still answers the
scientific question correctly despite it.

---

## 11. Figures

Every figure follows this pattern:

```r
fig_N_M <- <ggplot build expression>

print(fig_N_M)
ggplot2::ggsave(
  filename = "figs/Fig-N.M.svg",
  plot     = fig_N_M,
  width    = <W>, height = <H>, device = svg
)
```

`print()` before `ggsave()` is required so the figure appears in the
interactive session. `device = svg` (not `"svg"`) avoids the svglite
dependency.

Figure titles must match book captions exactly. If the book does not provide
a caption, construct one in the form "Figure N.M: [brief description]".

Apply `scale_x_continuous(breaks = ...)` on any plot whose x-axis covers a
narrow numeric range (e.g., probabilities below 0.05), because ggplot2's
default break algorithm may produce only one or two ticks.

---

## 12. Book comparison blocks

Every parameter with a stated book target gets a comparison block:

```r
log_result("\n=== M1: <model name> ===")
log_result("Book target: <param> ≈ X.XX (CI: Y.YY, Z.ZZ)")
capture_result(brms::fixef(fit_m1), label = "fixef(fit_m1)")
```

The comparison label must name the fit object and the function used, so the
reviewer can reproduce it without reading the full script.

---

## 13. Variance components and coefficient extraction

For variance components, `brms::VarCorr()` is correct. For fixed-effects-only
models (no random effects), `brms::VarCorr()` will error. Use
`brms::posterior_summary(fit_mN)["sigma", ]` for the residual standard
deviation instead.

For `brms::coef()` — note that `brms::coef()` is not always exported cleanly
from the brms namespace. If it errors with "could not find function", use
`brms::fixef(fit_mN)` plus `brms::ranef(fit_mN)` to reconstruct subject-level
coefficients manually.

---

## 14. Save and wrap up

Every script ends with the same block:

```r
# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete <- <TRUE or FALSE>
wf$loo_complete <- <TRUE or FALSE>
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

log_result("\n=== Ch N analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: <list of .rds files>")
log_result("Figures saved to figs/: <list of Fig-N.M.svg files>")
close(results_con)
```

`close(results_con)` flushes and closes `results.txt`. It must be the last
statement in the script. After it runs, `results.txt` in the analysis
subfolder contains a complete record of all book comparison targets, diagnostic
outcomes, LOO tables, and the file inventory.

The summary lines list every `.rds` file saved and every figure saved to
`figs/`, so the student can verify their output is complete.

---

## 15. What to omit

The following items must not appear in case study scripts:

- `pause()` calls
- `theme_set()` global theme calls
- `library(ggplot2)` or similar calls that duplicate what `source_all.R`
  loads, unless a package not in `source_all.R` is genuinely needed
- `moment_match = TRUE` in `add_criterion()` calls at fit time
- Bare `prior()` without the `brms::` prefix
- Bare `loo_compare()`, `loo()`, `pareto_k_values()` without `loo::`
- Bare `ppc_*()` functions without `bayesplot::`
- Bare `ggplot()`, `aes()`, `geom_*()`, `labs()` without `ggplot2::`
- `lw = weights(...)` in `ppc_loo_pit_overlay()` calls
- `ndraws = <N>` in `posterior_predict()` calls that feed a LOO-PIT plot
- `saveRDS()` omitted after any `brm()` call
- Figure save omitted (every figure must be saved to `figs/`)
- Bare `cat()` for book comparison targets, diagnostics, or LOO results — use
  `log_result()` / `capture_result()` so output lands in `results.txt`
- `close(results_con)` omitted at the end of the script
- A discrepancy between the posterior and a book target left unexplained
- `wf$loo_complete <- FALSE` with no accompanying `log_result()` reason

---

## 16. Results log

Every case study script writes a plain-text `results.txt` to the analysis
subfolder as it runs. This file records all numerical results — book targets,
computed posteriors, LOO comparisons, diagnostics — so the script author and
reviewer can inspect a single file rather than re-running the script.

### Setup (in §3, environment block)

```r
dir.create("figs", showWarnings = FALSE)
results_con <- file("results.txt", open = "wt")

# log_result(): text lines — writes to console and results.txt
log_result <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = results_con, append = FALSE)
}

# capture_result(): multi-line objects (fixef tables, loo_compare, etc.)
capture_result <- function(x, label = NULL) {
  txt <- paste(capture.output(print(x)), collapse = "\n")
  if (!is.null(label)) {
    cat(label, "\n", sep = "")
    cat(label, "\n", sep = "", file = results_con, append = FALSE)
  }
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = results_con, append = FALSE)
}

log_result("=== Ch N results log — ", format(Sys.time()), " ===")
```

`file(..., open = "wt")` opens `results.txt` for writing in text mode,
creating or overwriting it. The `append = FALSE` in both helpers is correct:
`cat(..., file = con, append = FALSE)` appends to an already-open connection
by default in R — `append` is ignored when `file` is an open connection, and
the connection itself handles buffering. The helpers are written this way so
they also work if `file` is replaced with a file path string during debugging.

### What to log

| What | Helper to use |
|---|---|
| Section banners, book targets, scalar diagnostics | `log_result()` |
| `brms::fixef()`, `brms::VarCorr()`, `brms::ranef()` output | `capture_result()` |
| `loo::loo_compare()` output | `capture_result()` |
| `priorsense::powerscale_sensitivity()` output | `capture_result()` |
| LOO Pareto-k counts | `log_result()` |
| Per-model diagnostic pass/fail | `log_result()` |
| File inventory at wrap-up | `log_result()` |
| Reason LOO was not run (when skipped) | `log_result()` |
| Explanation of any posterior-vs-book-target discrepancy | `log_result()` |

Do not log figure objects, raw draw matrices, or intermediate R objects.
Log the summary of the result, not the computation.

### Diagnostic log pattern

After every `run_diagnostics()` / `print(wf)` block, add:

```r
log_result("--- Diagnostics: fit_mN ---")
log_result("  passed      : ", wf$diagnostics$passed)
log_result("  Rhat_max    : ", round(wf$diagnostics$rhat_max,    4))
log_result("  bulk_ESS_min: ", round(wf$diagnostics$bulk_ess_min, 0))
log_result("  tail_ESS_min: ", round(wf$diagnostics$tail_ess_min, 0))
log_result("  divergences : ", wf$diagnostics$n_divergences)
if (length(wf$diagnostics$failed_criteria) > 0)
  log_result("  FAILED      : ", paste(wf$diagnostics$failed_criteria, collapse = ", "))
```

### LOO log pattern

After `loo::loo_compare()`:

```r
log_result("\n=== LOO comparison ===")
capture_result(loo_tab)
```

### Teardown (in §14, wrap-up block)

```r
log_result("\n=== Ch N analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: <list of .rds files>")
log_result("Figures: <list of Fig-N.M.svg files>")
close(results_con)
```

`close(results_con)` must be the last line in the script. After it runs,
`results.txt` is flushed and ready to read.

### Example results.txt excerpt

```
=== Ch 21 results log — 2026-09-19 14:32:01 ===

--- Diagnostics: bfit_0 ---
  passed      : TRUE
  Rhat_max    : 1.0008
  bulk_ESS_min: 3241
  tail_ESS_min: 2897
  divergences : 0

=== M0: Bernoulli logistic (compare to book: baseline) ===
            Estimate  Est.Error       Q2.5      Q97.5
Intercept  2.1815782 0.23037750  1.7313216  2.6289001
time      -0.2809347 0.02177642 -0.3241476 -0.2389321

=== LOO comparison ===
Book targets: M0h vs M0: elpd_diff -13, se_diff 5.7
              M0h vs M2: elpd_diff -4.7, se_diff 4.3
              M0h vs M4: elpd_diff -3.8, se_diff 3.1
   model elpd_diff se_diff p_worse
 bfit_0h       0.0     0.0      NA
  bfit_4      -3.5     3.1    0.87
  bfit_2      -4.3     4.2    0.84
  bfit_0     -12.2     5.6    0.98

=== Ch 21 analysis complete — 2026-09-19 16:47:22 ===
Saved: bfit_0.rds, bfit_0h.rds, bfit_2.rds, bfit_4.rds, loo_compare.rds, wf_final.rds
Figures: Fig-21.1.svg ... Fig-21.12.svg
```

---

## 17. Quick checklist before delivery

- [ ] All required Stan files and data files downloaded to the chapter folder before script was written
- [ ] Chapter folder contents verified against files referenced in the script
- [ ] Header lists all figures and all book comparison targets
- [ ] Stan-native chapters: header states why `run_phase3()` is bypassed and names every `.stan` file
- [ ] Stan-native chapters: latent-parameter models use `iter = 4000, warmup = 2000`
- [ ] `dir.create("figs", showWarnings = FALSE)` present in §3 setup block
- [ ] `results_con`, `log_result()`, `capture_result()` defined in §3 setup block
- [ ] `log_result("=== Ch N results log ...")` banner at start of §3
- [ ] Every `brm()` call has `refresh = 0` and `seed = SEED`
- [ ] Every `brm()` call is immediately followed by `saveRDS()`
- [ ] `save_pars = save_pars(all = TRUE)` on every fit that will get LOO
- [ ] `add_criterion(..., save_psis = TRUE)` without `moment_match`
- [ ] Conditional moment-matching loop present if LOO is run
- [ ] `loo::loo_compare()` and `loo::loo()` with explicit namespace
- [ ] `brms::posterior_predict()` with explicit namespace
- [ ] `brms::prior()` everywhere (not bare `prior()`)
- [ ] `ggplot2::` prefix on every ggplot2 function
- [ ] `bayesplot::` prefix on every bayesplot function
- [ ] `diagnose.wf_state(wf)` not `diagnose(wf)` or `wf$diagnose()`
- [ ] `fit_hash` set before `run_diagnostics()` when using `brm()` directly
- [ ] Diagnostic log block after every fit (brms or cmdstanr pattern as appropriate)
- [ ] Every posterior-vs-book-target discrepancy explained in `results.txt`
- [ ] LOO absence documented with `log_result("Reason: ...")` when skipped
- [ ] Non-standard PPC method explained in a comment block at Phase 5
- [ ] Simulated data recovery: note in `results.txt` whether mode matches true value and why any offset is expected
- [ ] Every figure: `print()` + `ggsave()` to `figs/Fig-N.M.svg`
- [ ] Figure titles match book captions exactly
- [ ] `scale_x_continuous(breaks = ...)` on any narrow-range faceted plot
- [ ] Book comparison block uses `log_result()` / `capture_result()` (not bare `cat()`)
- [ ] `loo::loo_compare()` result logged with `capture_result()`
- [ ] `priorsense::powerscale_sensitivity()` result logged with `capture_result()`
- [ ] `vapply` in data-fetch loops uses `integer(1L)` and `as.integer(i)`
- [ ] Wrap-up block uses `log_result()` for file inventory
- [ ] `close(results_con)` is the last line in the script
- [ ] `moment_match = TRUE` and `update() |> add_criterion()` pipes both fail on RTools45; always split into two calls.
