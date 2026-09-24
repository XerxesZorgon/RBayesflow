# =============================================================================
# RBayesflow learn-mode walkthrough
# Bayesian Workflow book, Chapter 18: Predictive model checking and comparison
# ---------------------------------------------------------------------------
# Case study: Nabiximols vs placebo in cannabis-use-disorder clinical trial
# Source of data: Lintzeris et al. (2019, 2020); vectors reproduced in
# Vehtari's public case-study repo.
#
# What this script covers  (roughly §18.1 through the end of §18.2):
#   Phase 1  Goal declaration + off-ramp assessment          (RBayesflow)
#   Phase 2  Prior specification + prior predictive
#   Phase 3  Fit three models: normal, binomial, beta-binomial
#            plus one refinement: beta-binomial with baseline as predictor
#   Phase 4  MCMC diagnostics (learn-mode gate fires on fit_binomial)
#   Phase 5  Posterior predictive checks + LOO-PIT-ECDF
#   Phase 6  LOO comparison across all four fits
#
# What is left for you (the student) to explore, deliberately:
#   §18.2 Power-scaling prior-likelihood sensitivity (priorsense)
#   §18.3 Treatment effect on a new individual (posterior_predict at id=129)
#   §18.4 Sensitivity to the data-generating model (normal vs beta-binomial
#         treatment-effect posteriors)
#   §18.5 Does the treatment-group variable matter? (elpd comparison of
#         beta-binomial with and without `group`)
#   §18.7 Exercises 18.1 and 18.2
#
# Run in RStudio with mode = "learn" (see init_workflow() call below).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Load RBayesflow (this walkthrough assumes you are in the RBayesflow
#    project root, with source_all.R and the wf_state / diagnostics / display
#    machinery in R/)
# -----------------------------------------------------------------------------

# Run from inside data/clinical_trial/ (the analysis subfolder). The
# .Rprofile there auto-sources source_all.R when RStudio sets that folder as
# the working directory. Running manually or in batch, call it explicitly:
source("../../R/source_all.R")

# Extra libraries used by this walkthrough (all are already in DESCRIPTION):
library(brms)
library(bayesplot)
library(loo)
library(posterior)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(tidyr)

options(brms.backend = "cmdstanr", mc.cores = 4)
ggplot2::theme_set(bayesplot::bayesplot::theme_default(base_family = "sans", base_size = 12))
SEED <- 1234


# -----------------------------------------------------------------------------
# 1. Data acquisition
#
# The Lintzeris et al. (2020) individual-level trial data are NOT posted with
# the JAMA paper. They live as hard-coded vectors inside Vehtari's public
# case-study script at:
#
#   https://github.com/avehtari/Bayesian-Workflow/blob/main/nabiximols/nabiximols.R
#
# We fetch that file once, cache it locally, and extract just the four
# data-definition assignments (id, group, week, cu). This avoids running the
# rest of Vehtari's script as a side effect.
# -----------------------------------------------------------------------------

cache_path <- "data/clinical_trial/lintzeris_2020_vectors.R"
dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_path)) {
  src_url <- paste0(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/",
    "main/nabiximols/nabiximols.R"
  )
  raw <- readLines(src_url, warn = FALSE)

  # Grab the four assignment blocks. Each starts with `id <-`, `group <-`,
  # `week <-`, or `cu <-` and can span multiple lines (group has a `labels =`
  # continuation). We stop each block at the next blank line.
  starts <- grep("^(id|group|week|cu) <-", raw)
  ends   <- vapply(starts, function(i) {
    j <- as.integer(i)
    while (j < length(raw) && nzchar(trimws(raw[j + 1L])) &&
           !grepl("^(id|group|week|cu|set|cu_df) <-", raw[j + 1L])) {
      j <- j + 1L
    }
    j
  }, integer(1L))

  data_lines <- unlist(Map(function(a, b) raw[a:b], starts, ends))
  writeLines(data_lines, cache_path)
  message("Cached data-definition lines to: ", cache_path)
}

# Evaluate the cached vectors into the current environment
source(cache_path, local = TRUE)
set   <- rep(28, length(cu))
cu_df <- data.frame(id, group, week, cu, set) |>
  tidyr::drop_na(cu)

cat("Data loaded:", nrow(cu_df), "rows,",
    length(unique(cu_df$id)), "participants,",
    length(unique(cu_df$week)), "time points.\n")


# -----------------------------------------------------------------------------
# Figure 1 (matches Fig 18.1 in the book) --- observed data
# -----------------------------------------------------------------------------
fig1_data <- cu_df |>
  ggplot2::ggplot(ggplot2::aes(x = cu)) +
  ggplot2::geom_histogram(breaks = seq(-0.5, 28.5, by = 1)) +
  ggplot2::scale_x_continuous(breaks = c(0, 10, 20, 28)) +
  ggplot2::facet_grid(group ~ week, switch = "y",
             labeller = ggplot2::labeller(group = ggplot2::label_value, week = ggplot2::label_both)) +
  ggplot2::labs(
    title    = "Figure 1: Observed cannabis-use days by group and week",
    subtitle = "128 participants, 28-day recall at 0, 4, 8, 12 weeks",
    x        = "Days of cannabis use in previous 28 days (cu)",
    y        = "Count",
    caption  = "Data: Lintzeris et al. (2020); reproduced via Vehtari (2024)."
  )
print(fig1_data)

# LEARN NOTE
# Take a moment before moving on. What structural features of this distribution
# would you expect to trip up a *continuous* model? A *binomial* model?
# (Hint: look at the piles at 0 and 28 — the "boundary bunching".)


# =============================================================================
# PHASE 1  --- Goal declaration and off-ramp assessment
# =============================================================================

wf <- init_workflow(mode = "learn", stage = "explore")

wf$declared_goal   <- "Estimate treatment effect of nabiximols vs placebo on days of cannabis use, over 12 weeks."
wf$n_observations  <- nrow(cu_df)
wf$event_rate      <- NA          # not a rare-event binary problem
wf$data_hash       <- digest::digest(cu_df, algo = "sha256")

# RBayesflow's off-ramp check. For a bounded count outcome at 128 participants
# with a hierarchical goal (repeated measures per participant), the decision
# matrix in DESIGN.md §4 raises the hierarchical warning row and offers
# alternatives with equal weight. It does NOT tell us Bayesian is best.
offramps <- assess_offramps(
  data         = cu_df,
  outcome_var  = "cu",
  outcome_type = "count",
  goal         = wf$declared_goal
)
print(offramps)

# LEARN NOTE
# The off-ramps offered here (GLMM with REML, quasi-likelihood, full Stan)
# are presented at equal weight on purpose. There is no "recommended" tag.
# We are choosing the Bayesian path because Chapter 18 is about *how* to
# check and improve a Bayesian model — not because Bayesian is automatically
# the right call for count data with n=128.
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase = 1, action = "bayesian_selected",
  timestamp = Sys.time(),
  rationale = "Case study exercise; goal is model criticism, not method selection."
)))


# =============================================================================
# PHASE 2  --- Priors and prior predictive simulation
#
# We will fit three models. Priors below match the book (§18.1 for the first
# two; §18.2 wider set for the beta-binomial refit).
# =============================================================================

priors_normal <- c(
  prior(normal(14, 1.5),  class = Intercept),
  prior(normal(0, 11),    class = b),
  prior(cauchy(1, 2),     class = sd)
)

priors_binomial <- c(
  prior(normal(0, 1.5),   class = Intercept),
  prior(normal(0, 1),     class = b),
  prior(cauchy(0, 2),     class = sd)
)

# The wider priors for the beta-binomial refit (book §18.2), motivated by
# power-scaling sensitivity results — see the "student explores" section
# at the end of this script.
priors_bb_wide <- c(
  prior(normal(0, 3),     class = Intercept),
  prior(normal(0, 3),     class = b),
  prior(normal(0, 3),     class = sd)
)

wf$priors_text <- c(
  "normal: N(14,1.5) intercept, N(0,11) beta, Cauchy(1,2) sd",
  "binomial: N(0,1.5) intercept, N(0,1) beta, Cauchy(0,2) sd",
  "beta-binomial (wider): N(0,3) intercept, N(0,3) beta, N(0,3) sd"
)
wf$priors_objects <- list(
  normal      = priors_normal,
  binomial    = priors_binomial,
  betabinom_w = priors_bb_wide
)

# LEARN NOTE (prior predictive)
# Chapter 18 does not lean heavily on prior predictive checks because the
# calibration argument is made after fitting. For your own workflow, you
# should still run one here — try:
#     brm(..., sample_prior = "only", ...)
# and see whether the prior alone produces reasonable cu values in [0, 28].


# =============================================================================
# PHASE 3  --- Fit the models
#
# We fit three models straight from §18.1–18.2:
#   fit_normal        --- Gaussian likelihood (§18.1)
#   fit_binomial      --- Binomial(28, .) with logit link (§18.1)
#   fit_betabinomial  --- Beta-binomial refit with wider priors (§18.2)
#   fit_betabinomial2b --- Same but week=0 moved to a baseline predictor,
#                          removing the meaningless treatment*week=0 term
# =============================================================================

fit_normal <- brm(
  formula = cu ~ group * week + (1 | id),
  data    = cu_df, family = gaussian(),
  prior   = priors_normal,
  save_pars = save_pars(all = TRUE),
  seed = SEED, refresh = 0
)
fit_normal <- add_criterion(fit_normal, "loo", save_psis = TRUE)

fit_binomial <- brm(
  formula = cu | trials(set) ~ group * week + (1 | id),
  data    = cu_df, family = binomial(link = "logit"),
  prior   = priors_binomial,
  save_pars = save_pars(all = TRUE),
  seed = SEED, refresh = 0
)
fit_binomial <- add_criterion(fit_binomial, "loo", save_psis = TRUE)

fit_betabinomial <- brm(
  formula = cu | trials(set) ~ group * week + (1 | id),
  data    = cu_df, family = beta_binomial(),
  prior   = priors_bb_wide,
  save_pars = save_pars(all = TRUE),
  seed = SEED, refresh = 0
)
fit_betabinomial <- add_criterion(fit_betabinomial, "loo", save_psis = TRUE)

# Build the "week=0 as baseline predictor" data frame per §18.2 end.
baseline <- cu_df |>
  filter(week == 0) |>
  select(id, cu_baseline = cu)

cu_df_b <- cu_df |>
  filter(week != 0) |>
  mutate(week = droplevels(week)) |>
  left_join(baseline, by = "id")

fit_betabinomial2b <- brm(
  formula = cu | trials(set) ~ group * week + cu_baseline + (1 | id),
  data    = cu_df_b, family = beta_binomial(),
  prior   = c(prior(normal(0, 3), class = Intercept),
              prior(normal(0, 3), class = b)),
  save_pars = save_pars(all = TRUE),
  seed = SEED, refresh = 0
)
fit_betabinomial2b <- add_criterion(fit_betabinomial2b, "loo", save_psis = TRUE)

# Moment matching: apply only where needed, and only after all fits are done.
# moment_match = TRUE triggers an rstan recompile on Windows, which fails if
# the cmdstanr/rstan session state is mixed. The safe pattern (following
# Vehtari's own script) is to run plain LOO first, inspect Pareto-k values,
# and re-run with moment_match only for fits that have problematic k̂ > 0.7.
for (fit_name in c("fit_normal", "fit_binomial",
                   "fit_betabinomial", "fit_betabinomial2b")) {
  fit_obj   <- get(fit_name)
  loo_obj   <- fit_obj$criteria$loo
  n_bad_k   <- sum(loo::pareto_k_values(loo_obj) > 0.7)
  if (n_bad_k > 0) {
    message(fit_name, ": ", n_bad_k,
            " obs with Pareto-k > 0.7 — applying moment matching.")
    fit_obj  <- add_criterion(fit_obj, "loo",
                              save_psis = TRUE, moment_match = TRUE)
    assign(fit_name, fit_obj)
  } else {
    message(fit_name, ": LOO reliable (no Pareto-k > 0.7), skipping moment match.")
  }
}

# Register the "best so far" fit on the workflow state so the display
# contract and phase-4 diagnostics use it.
wf$formula         <- formula(fit_betabinomial2b)
wf$family          <- beta_binomial()
wf$fit_timestamp   <- Sys.time()
wf$stan_backend    <- "cmdstanr"
wf$parameterization <- detect_parameterization(fit_betabinomial2b)
wf$fit_hash <- digest::digest(
  list(formula   = as.character(formula(fit_betabinomial2b)),
       data_hash = digest::digest(fit_betabinomial2b$data, algo = "sha256"),
       timestamp = wf$fit_timestamp),
  algo = "sha256"
)


# =============================================================================
# PHASE 4  --- MCMC diagnostics
#
# This is where the RBayesflow diagnostic gate really shows what it is for.
# We run diagnostics on fit_binomial FIRST (the misspecified one), so you can
# see the coefficient-withholding behaviour in learn mode. Then we run them
# on fit_betabinomial2b (the good one) so the display contract flips.
# =============================================================================

# --- 4a. Diagnostics on the misspecified binomial model ---------------------
wf_bin <- wf
wf_bin$formula         <- formula(fit_binomial)
wf_bin$family          <- binomial()
wf_bin$fit_timestamp   <- fit_binomial$fit$metadata()$start_datetime  |>
  tryCatch(error = function(e) Sys.time())
# record_fit() links wf to the exact fit via SHA-256 hash (ADR-009).
# We set fit_hash directly here because record_fit() also calls
# export_context(), which needs the analysis subfolder to exist.
wf_bin$fit_hash <- digest::digest(
  list(formula   = as.character(formula(fit_binomial)),
       data_hash = digest::digest(fit_binomial$data, algo = "sha256"),
       timestamp = wf_bin$fit_timestamp),
  algo = "sha256"
)
wf_bin <- run_diagnostics(fit_binomial, wf_bin)

cat("\n--- print(wf_bin) BEFORE acknowledgment ---\n")
print(wf_bin)
# LEARN MODE, DIAGNOSTICS FAILED, NOT YET ACKNOWLEDGED:
# You should see a failure message and a prompt to call wf$diagnose().
# You should NOT see any coefficient table, and (importantly for learn mode)
# no prior-vs-posterior overlay either — the fit is untrustworthy.

# The book's Figure 18.3 covers the two objects we care about here: PPC
# histograms and the LOO-PIT-ECDF calibration curve. We draw them explicitly
# so you can see the calibration problem.
yrep_bin <- brms::posterior_predict(fit_binomial)  # all 4000 draws — must match PSIS object

# Two complementary views of the binomial model's PPC failure:
#
# Panel A — density overlay: observed (dark) vs 50 posterior predictive draws
# (light). One panel; all draws superimposed. Shows the overall shape mismatch.
fig_ppc_bin_dens <- bayesplot::ppc_dens_overlay(cu_df$cu, yrep_bin[1:50, ]) +
  ggplot2::scale_x_continuous(limits = c(-1, 29), breaks = c(0, 7, 14, 21, 28)) +
  ggplot2::labs(
    title    = "Figure 2a: Posterior predictive density overlay, binomial model",
    subtitle = "Light curves = posterior predictive draws; dark = observed data.\nThe model compresses mass away from 0 and 28.",
    x = "Days of cannabis use (cu)", y = "Density"
  )
print(fig_ppc_bin_dens)

# Panel B — boundary mass statistics: what fraction of observations are exactly
# 0 or exactly 28? The model should match these proportions.
stat_zero <- function(y) mean(y == 0)
stat_max  <- function(y) mean(y == 28)

fig_ppc_bin_zeros <- bayesplot::ppc_stat(cu_df$cu, yrep_bin,
                               stat = "stat_zero", binwidth = 0.01) +
  ggplot2::labs(
    title    = "Figure 2b: Boundary mass at 0, binomial model",
    subtitle = "T(y) = proportion of observations equal to 0.\nIf the bar (observed) falls in the tail, the model under-predicts zeros.",
    x = "Proportion of cu = 0", y = NULL
  )
print(fig_ppc_bin_zeros)

fig_ppc_bin_maxes <- bayesplot::ppc_stat(cu_df$cu, yrep_bin,
                               stat = "stat_max", binwidth = 0.01) +
  ggplot2::labs(
    title    = "Figure 2c: Boundary mass at 28, binomial model",
    subtitle = "T(y) = proportion of observations equal to 28.\nIf the bar falls in the tail, the model under-predicts 28s.",
    x = "Proportion of cu = 28", y = NULL
  )
print(fig_ppc_bin_maxes)

fig_loopit_bin <- bayesplot::ppc_loo_pit_overlay(
  y           = cu_df$cu,
  yrep        = yrep_bin,
  psis_object = fit_binomial$criteria$loo$psis_object
) +
  ggplot2::labs(
    title    = "Figure 3: LOO-PIT calibration, binomial model",
    subtitle = "Piles near 0 and 1 => predictive intervals are too narrow (underdispersion).",
    x = "LOO-PIT", y = "Density"
  )
print(fig_loopit_bin)

# Now formally acknowledge — this is what a real learn-mode session looks like.
# In an interactive session you would type "yes" at the prompt. For a batch
# run, set acknowledged directly (with an audit-trail note) so the script
# continues.
if (interactive()) {
  wf_bin <- diagnose.wf_state(wf_bin)
} else {
  wf_bin$diagnostics$acknowledged <- TRUE
  wf_bin$audit_trail <- append(wf_bin$audit_trail, list(list(
    phase = 4, action = "diagnostic_acknowledged_noninteractive",
    timestamp = Sys.time(),
    failed_criteria = wf_bin$diagnostics$failed_criteria
  )))
}

cat("\n--- print(wf_bin) AFTER acknowledgment ---\n")
print(wf_bin)
# LEARN MODE, ACKNOWLEDGED BUT STILL FAILED:
# The failure message repeats with an "(acknowledged)" note. Coefficients are
# still not shown in learn mode — because the point of learn mode is to force
# the user to *understand* that a failing model's coefficients are not what
# they look like.


# --- 4b. Diagnostics on the beta-binomial refit ------------------------------
wf <- run_diagnostics(fit_betabinomial2b, wf)

cat("\n--- print(wf) for fit_betabinomial2b ---\n")
print(wf)
# LEARN MODE, DIAGNOSTICS PASSED:
# You should now see the prior-vs-posterior overlay plot and a one-line
# health summary. No coefficient table — that is *practice* mode. Learn
# mode is teaching you to look at the whole posterior first.


# =============================================================================
# PHASE 5  --- Posterior predictive checks and LOO-PIT-ECDF
#
# We reproduce the equivalent of Fig 18.5 for the beta-binomial model.
# =============================================================================

yrep_bb <- brms::posterior_predict(fit_betabinomial2b)  # all 4000 draws — must match PSIS object

# Same two panels for the beta-binomial refit — compare directly to Figures 2a–c.
fig_ppc_bb_dens <- bayesplot::ppc_dens_overlay(cu_df_b$cu, yrep_bb[1:50, ]) +
  ggplot2::scale_x_continuous(limits = c(-1, 29), breaks = c(0, 7, 14, 21, 28)) +
  ggplot2::labs(
    title    = "Figure 4a: Posterior predictive density overlay, beta-binomial refit",
    subtitle = "Light curves = posterior predictive draws; dark = observed data.\nBoundary bunching at 0 and 28 now captured.",
    x = "Days of cannabis use (cu)", y = "Density"
  )
print(fig_ppc_bb_dens)

fig_ppc_bb_zeros <- bayesplot::ppc_stat(cu_df_b$cu, yrep_bb,
                              stat = "stat_zero", binwidth = 0.01) +
  ggplot2::labs(
    title    = "Figure 4b: Boundary mass at 0, beta-binomial refit",
    subtitle = "Observed proportion of cu = 0 should now fall inside the predictive distribution.",
    x = "Proportion of cu = 0", y = NULL
  )
print(fig_ppc_bb_zeros)

fig_ppc_bb_maxes <- bayesplot::ppc_stat(cu_df_b$cu, yrep_bb,
                              stat = "stat_max", binwidth = 0.01) +
  ggplot2::labs(
    title    = "Figure 4c: Boundary mass at 28, beta-binomial refit",
    subtitle = "Observed proportion of cu = 28 should now fall inside the predictive distribution.",
    x = "Proportion of cu = 28", y = NULL
  )
print(fig_ppc_bb_maxes)

fig_loopit_bb <- bayesplot::ppc_loo_pit_overlay(
  y           = cu_df_b$cu,
  yrep        = yrep_bb,
  psis_object = fit_betabinomial2b$criteria$loo$psis_object
) +
  ggplot2::labs(
    title    = "Figure 5: LOO-PIT calibration, beta-binomial refit",
    subtitle = "Curve close to uniform => calibrated predictive intervals.",
    x = "LOO-PIT", y = "Density"
  )
print(fig_loopit_bb)

wf$ppc_complete <- TRUE
wf$ppc_summary  <- list(
  n_draws_shown   = 50,
  plot_types      = c("ppc_dens_overlay", "ppc_stat (zero mass)", "ppc_stat (28 mass)"),
  loo_pit_used    = TRUE
)


# =============================================================================
# PHASE 6  --- LOO comparison across all four fits
#
# One subtlety worth internalising: fit_betabinomial2b uses cu_df_b (week=0
# rows removed), so its per-observation elpd is not directly comparable to
# the first three fits on their raw scale. §18.2 handles this by comparing
# the beta-binomials to each other, and comparing the first three separately.
# We do the same.
# =============================================================================

cat("\n--- LOO comparison: models fit on the full cu_df ---\n")
print(loo::loo_compare(fit_normal, fit_binomial, fit_betabinomial))

cat("\n--- LOO comparison: beta-binomial vs beta-binomial+baseline ---\n")
# Not identical to the book's approach (which restricts pointwise elpd of
# fit_betabinomial to week>0 rows before comparing). Kept simple here on
# purpose — see the exploration list below.
print(loo::loo(fit_betabinomial2b))

wf$loo_complete <- TRUE
wf$loo_table    <- loo::loo_compare(fit_normal, fit_binomial, fit_betabinomial)

# LEARN NOTE
# Two things to notice:
#   1. The binomial model loses badly to the normal, even though the outcome
#      IS a count in [0, 28]. Chapter 18's point is that log-score cares
#      about calibration, and an underdispersed count model calibrates
#      worse than a wrong-family continuous one.
#   2. The beta-binomial then beats the normal by a wide margin. This is
#      the payoff of adding overdispersion.


# =============================================================================
# HANDOFF --- what you (the student) explore from here
#
# The book takes the story four more steps, all of which are more instructive
# once you have seen learn mode work through the first three phases:
#
#   §18.2 (end) Power-scaling prior-likelihood sensitivity.
#         Install {priorsense} and run:
#             priorsense::powerscale_sensitivity(fit_betabinomial)
#         Compare to fit_betabinomial2 with wider priors. Reproduce the
#         table in Fig 18.7 and explain why "potential conflict" does not
#         mean "wrong prior".
#
#   §18.3 Treatment effect on a new individual.
#         Use posterior_predict() on a new id=129 with cu_baseline=28 and
#         reproduce Figs 18.9 and 18.11 (aleatoric-included vs expectation-
#         only). Add titles.
#
#   §18.4 Data-model sensitivity.
#         Fit fit_normal2b (Gaussian, baseline predictor, wider priors) and
#         compare its treatment-effect posterior to the beta-binomial's.
#         Reproduce Fig 18.12. Why does the normal underestimate the effect
#         magnitude AND report a narrower posterior?
#
#   §18.5 Does treatment matter?
#         Fit fit_betabinomial3b without `group`. Compare LOO. Then compare
#         mean absolute error of predictive means (dropping aleatoric
#         uncertainty). Bridge-sampling Bayes factor is optional.
#
#   §18.7 Exercises 18.1 and 18.2 — designed for this exact workflow.
#
# When you finish, run:
#     wf <- run_phase7_report(wf)   # RBayesflow's Quarto report generator
# to produce a fully reproducible bayesflow_report.pdf/html.
# =============================================================================

# Save the workflow state so a follow-up session can pick up where you stopped.
saveRDS(wf, file = "data/clinical_trial/wf_after_phase6.rds")
saveRDS(list(
  fit_normal          = fit_normal,
  fit_binomial        = fit_binomial,
  fit_betabinomial    = fit_betabinomial,
  fit_betabinomial2b  = fit_betabinomial2b
), file = "data/clinical_trial/fits_ch18.rds")

export_context(wf)   # writes wf_context.json in the current analysis folder

cat("\nEnd of Chapter 18 walkthrough (§18.1 – §18.2). ",
    "Next session picks up from §18.3.\n", sep = "")
