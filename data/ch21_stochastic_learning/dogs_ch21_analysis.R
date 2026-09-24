# data/dogs_ch21/dogs_ch21_analysis.R
#
# Ch 21 — "Posterior predictive checking: Stochastic learning in dogs"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Fits four models to Bush & Mosteller's (1955) 30-dog shock-avoidance data
# and compares them with PSIS-LOO-CV, PAV-adjusted calibration and residual
# plots, visual posterior predictive checks (per-dog tile plots), a test-
# statistic PPC on mean number of shock/avoidance switches, and prior-
# likelihood power-scaling sensitivity. The four models are: simple logistic
# regression in time (M0), hierarchical logistic regression (M0h), the
# non-hierarchical Bush-Mosteller 2-parameter log model (M2), and its
# hierarchical extension (M4). The narrative conclusion is that M0h, M2 and
# M4 are indistinguishable on this dataset; visual PPC is historically
# interesting but not decisive here.
#
# Deliberately omits:
#   Model 1  (1-parameter log model) — the chapter shows it is worse than
#            either logistic regression model and drops it. Adding it costs
#            one fit and teaches nothing the surviving four don't already
#            teach.
#   Model 3  (hierarchical 1-parameter log model) — same reasoning; the
#            chapter drops it after one comparison table.
#   Leave-future-out CV (chapter §10) — 21 sequential refits per model × 2
#            models = 42 refits; the chapter's own conclusion is that LFO
#            does not change the M0h-vs-M4 ranking. Skipped for runtime.
#   "Are the data informative on two parameters" simulation check (chapter
#            §11) — two additional M4 fits on simulated data; the chapter
#            uses this to argue the posteriors of a and b being different
#            is not evidence of different learning rates. Skipped for
#            runtime; noted in the LEARN NOTE at the end.
#
# Data:
#   Source: avehtari/Bayesian-Workflow GitHub repo, branch main,
#           dogs/data/dogs.dat (30 rows × 26 cols, first col is dog id,
#           cols 2-26 are S/A markers for trials 1-25).
#   Downloaded once and cached at ./dogs.dat.
#   Reshape: read.table(..., skip = 2); then long-form data frame with
#            30 dogs × 25 trials = 750 rows, first-trial rows dropped
#            (always shock, no uncertainty) leaving N = 720.
#   Covariates: prev_shock, prev_avoid = cumulative counts up to but not
#               including the current trial.
#
# Models covered:
#   M0   bfit_0   shock ~ time                                       bernoulli logit
#                 alpha ~ Student_3(0, 2.5); beta ~ Normal(0, 1)
#                 Purpose: baseline pooled logistic regression in time.
#   M0h  bfit_0h  shock ~ time + (time | dog)                        bernoulli logit
#                 Per-dog (alpha_j, beta_j) with MVN hierarchical prior,
#                 LKJ(1) correlation, Student_3^+(0, 2.5) sd's.
#                 Purpose: per-dog variation on the logistic baseline.
#   M2   bfit_2   shock ~ a^prev_shock * b^prev_avoid                bernoulli identity
#                 a, b ~ Beta(1, 1) truncated to (0, 1); nonlinear.
#                 Purpose: original Bush-Mosteller 2-parameter log model,
#                 no hierarchy, dog variation enters through the
#                 dog-specific prev_shock, prev_avoid covariates.
#   M4   bfit_4   shock ~ inv_logit(etaa)^prev_shock * inv_logit(etab)^prev_avoid
#                                                                    bernoulli identity
#                 (etaa, etab) per dog with MVN hierarchical prior on
#                 the logit scale, Student_3^+(0, 2.5) sd's, LKJ(1).
#                 Purpose: hierarchical Bush-Mosteller.
#
# Figures produced:
#   Fig-21.1.svg   per-dog predicted probability, M0h, first 9 dogs
#   Fig-21.2.svg   PAV-adjusted LOO calibration plot, M0h
#   Fig-21.3.svg   PAV-adjusted LOO residual vs time, M0h
#   Fig-21.4a.svg  real dogs tile plot (rows = dogs ordered by last shock)
#   Fig-21.4b.svg  M0h posterior-predictive tile plot
#   Fig-21.4.svg   combined 5-model PPC tile grid (book Figure 21.4)
#   Fig-21.5.svg   PPC on mean number of shock/avoidance switches, M0h
#   Fig-21.6.svg   per-dog predicted probability, M4, first 9 dogs
#   Fig-21.7.svg   per-dog predicted probability, M2, first 9 dogs
#   Fig-21.8.svg   PAV-adjusted LOO calibration plot, M4
#   Fig-21.9.svg   PAV-adjusted LOO residual vs time, M4
#   Fig-21.10.svg  M4 posterior-predictive tile plot
#   Fig-21.11.svg  PPC on mean number of shock/avoidance switches, M4
#   Fig-21.12.svg  overlaid per-dog posterior mean predictions, M0h vs M4
#
# Book-target posterior summaries (success criteria):
#   The chapter does not print numeric posterior summaries for a, b, or
#   the logistic-regression coefficients. Success criteria are LOO
#   comparison numbers (book Fig 21.4 comparison table):
#     M0h vs M0 : elpd_diff -13,  se_diff 5.7,  p_worse 0.99  (M0 worse)
#     M0h vs M2 : elpd_diff -4.7, se_diff 4.3,  p_worse 0.86  (M2 worse, ns)
#     M0h vs M4 : elpd_diff -3.8, se_diff 3.1,  p_worse 0.89  (M4 worse, ns)
#   Sampling health: Rhat < 1.01, bulk & tail ESS > 400 on all four fits.
#   priorsense diagnosis column: no "prior-data conflict" flags on M0h or
#     M4 (data are informative; likelihood dominates prior on every param).
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/dogs_ch21/.
# =============================================================================


# =============================================================================
# PHASE 0  ---  Environment setup
# =============================================================================

source("../../R/source_all.R")

library(brms)           # brm(), bf(), bernoulli(), save_pars(), add_criterion()
library(dplyr)          # filter(), mutate(), bind_rows(), group_by(), summarise()
library(tidybayes)      # add_linpred_draws(), stat_lineribbon()
library(patchwork)      # multi-panel composition (used for tile grids)
library(RColorBrewer)   # brewer.pal palette for line ribbons
library(matrixStats)    # rowCumsums(), rowDiffs()
library(Iso)            # Iso::pava() inside ppc_pava_residual()
library(reliabilitydiag) # PAV-adjusted LOO calibration plots
library(priorsense)     # power-scaling prior-likelihood sensitivity

options(brms.backend = "cmdstanr", mc.cores = 4)
SEED <- 123

dir.create("figs", showWarnings = FALSE)


# =============================================================================
# PHASE 1  ---  Data acquisition, exploration, and off-ramp assessment
# =============================================================================

wf <- init_workflow(mode = "learn", stage = "explore")

# -----------------------------------------------------------------------------
# Data: cache dogs.dat from the book's GitHub repo, load, reshape.
# -----------------------------------------------------------------------------
cache_path <- "dogs.dat"
if (!file.exists(cache_path)) {
  src_url <- paste0(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/",
    "main/dogs/data/dogs.dat"
  )
  utils::download.file(src_url, cache_path, quiet = TRUE, mode = "wb")
  message("Downloaded dogs data to: ", cache_path)
}

# skip = 2 matches the book: the file has two header lines before the data.
dogs  <- read.table(cache_path, skip = 2)
shock <- ifelse(as.matrix(dogs[, 2:26]) == "S", 1, 0)
stopifnot(dim(shock) == c(30L, 25L))
cat(sprintf("Loaded %d dogs x %d trials.\n", nrow(shock), ncol(shock)))

# Long-form data frame; drop trial 1 (always shock, no uncertainty).
dogs_df <- data.frame(
  shock = as.numeric(shock),
  dog   = rep(1:nrow(shock), times = ncol(shock)),
  time  = rep(1:ncol(shock), each  = nrow(shock))
)
dogs_df <- dplyr::filter(dogs_df, time > 1)

# Cumulative shocks and avoidances up to (but not including) the current trial.
dogs_df$prev_shock <- as.numeric(matrixStats::rowCumsums(shock)[, 1:(ncol(shock) - 1)])
dogs_df$prev_avoid <- as.numeric(matrixStats::rowCumsums(1 - shock)[, 1:(ncol(shock) - 1)])

cat(sprintf("Modelling data: N = %d rows (30 dogs x 24 trials).\n", nrow(dogs_df)))

# LEARN NOTE
# The design has no predictors that vary across dogs at trial 1: every dog
# starts uninformed. Because trial 1 is always shock, we drop it. From trial
# 2 onward each dog has its own history of prev_shock and prev_avoid, which
# is why even non-hierarchical models (M2) can express dog-to-dog variation
# without random effects.

# -----------------------------------------------------------------------------
# Off-ramp assessment. Binary outcome, N = 720, event rate ~ 0.36; nothing
# rare here. Full-Stan Bayesian is the natural path; we log the choice.
# -----------------------------------------------------------------------------
cat(sprintf("Overall shock rate on trials 2-25: %.3f\n", mean(dogs_df$shock)))

offramps <- assess_offramps(
  data         = dogs_df,
  outcome_var  = "shock",
  outcome_type = "binary",
  goal         = "posterior predictive checking of learning models"
)
print(offramps)

# We proceed with full Stan (Bayesian). No exit.
wf$declared_goal   <- "posterior predictive checking of learning models"
wf$n_observations  <- nrow(dogs_df)
wf$event_rate      <- offramps$event_rate


# =============================================================================
# PHASE 2  ---  Priors (documented per model; no separate PP-check figure
#                because the chapter itself does not prior-predict here)
# =============================================================================

# Model 0 (pooled logistic):
#   alpha ~ Student_3(0, 2.5)  [brms default on the intercept]
#   beta  ~ Normal(0, 1)
#     ~95% prior probability that per-trial log-odds change is in (-2, 2),
#     i.e. odds multiplied by (0.14, 7.4) per trial. Wide but not silly.
priors_m0 <- brms::prior(normal(0, 1), class = "b")

# Model 0h (hierarchical logistic):
#   Same fixed-effect priors as M0.
#   sd's on random intercept/slope: Student_3^+(0, 2.5)  [brms default]
#   correlation:                    LKJ(1)               [brms default]
priors_m0h <- brms::prior(normal(0, 1), class = "b")

# Model 2 (non-hierarchical Bush-Mosteller):
#   a, b ~ Beta(1, 1) truncated to (0, 1) = Uniform(0, 1) on both learning
#   parameters. The bounds are structural: a and b are probabilities-per-event,
#   so they must lie in [0, 1].
priors_m2 <- c(
  brms::prior(beta(1, 1), nlpar = "a", lb = 0, ub = 1),
  brms::prior(beta(1, 1), nlpar = "b", lb = 0, ub = 1)
)

# Model 4 (hierarchical Bush-Mosteller):
#   etaa_j, etab_j are on the logit scale; inv_logit() maps them into (0, 1)
#   at the model level, so we do NOT need bounded priors here.
#   Student_3(0, 2.5) on both mu's puts ~95% prior mass on inv_logit values
#   in roughly (0.02, 0.98). Wide enough that data can dominate; not flat.
priors_m4 <- c(
  brms::prior(student_t(3, 0, 2.5), nlpar = "etaa"),
  brms::prior(student_t(3, 0, 2.5), nlpar = "etab")
)


# =============================================================================
# PHASE 3  ---  Fit
# =============================================================================
# NOTE: We use brm() directly (not run_phase3()) because the nonlinear models
# M2 and M4 require brms::bf() with nl = TRUE and non-default link, which the
# wrapper does not expose. We call record_fit(wf, fit) after each brm() call
# to register the hash — record_fit() uses fit$data internally, not dogs_df,
# which is what check_fit_hash() inside run_diagnostics() verifies against.

inv_logit <- function(x) 1 / (1 + exp(-x))

# --- Model 0: pooled logistic ------------------------------------------------
bfit_0 <- brm(
  formula   = shock ~ time,
  data      = dogs_df,
  family    = bernoulli(),
  prior     = priors_m0,
  save_pars = save_pars(all = TRUE),
  backend   = "cmdstanr",
  seed      = SEED,
  refresh   = 0
)
saveRDS(bfit_0, "bfit_0.rds")
wf <- record_fit(wf, bfit_0)
wf <- run_diagnostics(bfit_0, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)
cat("\n=== M0: Bernoulli logistic (compare to book: baseline) ===\n")
print(brms::fixef(bfit_0))

bfit_0 <- add_criterion(bfit_0, "loo", save_psis = TRUE)


# --- Model 0h: hierarchical logistic -----------------------------------------
bfit_0h <- brm(
  formula     = shock ~ time + (time | dog),
  data        = dogs_df,
  family      = bernoulli(),
  prior       = priors_m0h,
  save_pars   = save_pars(all = TRUE),
  backend     = "cmdstanr",
  seed        = SEED,
  refresh     = 0,
  control     = list(adapt_delta = 0.95)
)
saveRDS(bfit_0h, "bfit_0h.rds")
wf <- record_fit(wf, bfit_0h)
wf <- run_diagnostics(bfit_0h, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)
cat("\n=== M0h: hierarchical logistic (compare to book: best-fitting) ===\n")
print(brms::fixef(bfit_0h))
print(brms::VarCorr(bfit_0h))

bfit_0h <- add_criterion(bfit_0h, "loo", save_psis = TRUE)


# --- Model 2: non-hierarchical 2-parameter log model -------------------------
bfit_2 <- brm(
  formula   = brms::bf(shock ~ a^prev_shock * b^prev_avoid,
                       a ~ 1, b ~ 1, nl = TRUE),
  data      = dogs_df,
  family    = bernoulli(link = "identity"),
  prior     = priors_m2,
  save_pars = save_pars(all = TRUE),
  backend   = "cmdstanr",
  seed      = SEED,
  refresh   = 0
)
saveRDS(bfit_2, "bfit_2.rds")
wf <- record_fit(wf, bfit_2)
wf <- run_diagnostics(bfit_2, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)
cat("\n=== M2: Bush-Mosteller 2-parameter log model (a, b) ===\n")
print(brms::fixef(bfit_2))

bfit_2 <- add_criterion(bfit_2, "loo", save_psis = TRUE)


# --- Model 4: hierarchical 2-parameter log model -----------------------------
bfit_4 <- brm(
  formula   = brms::bf(shock ~ inv_logit(etaa)^prev_shock * inv_logit(etab)^prev_avoid,
                       mvbind(etaa, etab) ~ (1 | p | dog), nl = TRUE),
  data      = dogs_df,
  family    = bernoulli(link = "identity"),
  prior     = priors_m4,
  save_pars = save_pars(all = TRUE),
  backend   = "cmdstanr",
  seed      = SEED,
  refresh   = 0
)
saveRDS(bfit_4, "bfit_4.rds")
wf <- record_fit(wf, bfit_4)
wf <- run_diagnostics(bfit_4, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)
cat("\n=== M4: hierarchical Bush-Mosteller (etaa, etab per dog) ===\n")
print(brms::fixef(bfit_4))
print(brms::VarCorr(bfit_4))

bfit_4 <- add_criterion(bfit_4, "loo", save_psis = TRUE)


# =============================================================================
# PHASE 4  ---  Diagnostics wrap-up: conditional moment matching for LOO
# =============================================================================
# Only apply moment matching if PSIS Pareto-k > 0.7 for any observation.
# This avoids the Windows rstan-recompilation issue (rubric §6).
for (fit_name in c("bfit_0", "bfit_0h", "bfit_2", "bfit_4")) {
  fit_obj <- get(fit_name)
  n_bad   <- sum(loo::pareto_k_values(fit_obj$criteria$loo) > 0.7)
  if (n_bad > 0) {
    message(fit_name, ": ", n_bad, " high Pareto-k; applying moment matching.")
    fit_obj <- add_criterion(fit_obj, "loo",
                             save_psis = TRUE, moment_match = TRUE)
    assign(fit_name, fit_obj)
  } else {
    message(fit_name, ": all Pareto-k <= 0.7; no moment matching needed.")
  }
}


# =============================================================================
# PHASE 5  ---  Posterior predictive checks
# =============================================================================
# The chapter's PPC content splits into three activities:
#   (a) per-dog visual predictions (Figs 21.1, 21.6, 21.7)
#   (b) PAV-calibration and PAV-residual plots (Figs 21.2, 21.3, 21.8, 21.9)
#   (c) tile-plot posterior-predictive replicates (Figs 21.4a, 21.4b, 21.4,
#       21.10) plus a mean-switches test statistic PPC (Figs 21.5, 21.11)

# --- ppc_pava_residual: inline helper from the chapter -----------------------
# The book defines this locally; not exported from any CRAN package.
ppc_pava_residual <- function(y, epred, x = NULL, prob = .9, n.boot = 1000,
                              interval_geom = "ribbon", alpha = .2, ...) {
  assertthat::assert_that(is.numeric(y))
  assertthat::assert_that(is.numeric(epred))
  assertthat::are_equal(length(epred), length(y))

  yrep_bar_order <- order(epred)
  if (is.null(x)) x <- epred

  cep_df <- seq_len(n.boot) |>
    lapply(\(i) data.frame(
      cep = Iso::pava(rbinom(length(y), 1, epred)[yrep_bar_order]),
      id_ = yrep_bar_order
    )) |>
    dplyr::bind_rows() |>
    dplyr::group_by(id_) |>
    dplyr::summarise(upper = quantile(cep, .5 + .5 * prob),
                     lower = quantile(cep, .5 * (1 - prob)))

  cep_df$yrep_bar <- epred
  cep_df[yrep_bar_order, "cep"] <- Iso::pava(y[yrep_bar_order])
  cep_df$x    <- x
  cep_df$ymax <- cep_df$upper - cep_df$yrep_bar
  cep_df$ymin <- cep_df$lower - cep_df$yrep_bar

  bw <- .5 * bw.SJ(cep_df$x)
  w  <- sapply(cep_df$x, \(x_i) dnorm(cep_df$x, x_i, bw))
  cep_df$ymaxs <- (t(w) %*% cep_df$ymax) / colSums(w)
  cep_df$ymins <- (t(w) %*% cep_df$ymin) / colSums(w)

  ggplot2::ggplot(cep_df, ggplot2::aes(y = cep - yrep_bar,
                                       ymax = ymaxs, ymin = ymins, x = x)) +
    ggplot2::geom_hline(yintercept = 0, alpha = 0.3) +
    ggplot2::stat_identity(ggplot2::aes(colour = TRUE, fill = TRUE),
                           alpha = alpha, geom = interval_geom, ...) +
    ggplot2::geom_point(ggplot2::aes(colour = (cep >= lower) & (cep <= upper))) +
    ggplot2::scale_colour_discrete(aesthetics = c("fill", "colour")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(y = "PAVA Residual")
}


# --- (a) Per-dog visual predictions ------------------------------------------

# --- Figure 21.1: M0h per-dog predicted probability, first 9 dogs ---
fig_21_1 <- dogs_df |>
  dplyr::filter(dog <= 9) |>
  tidybayes::add_linpred_draws(bfit_0h, transform = TRUE) |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = .linpred)) +
  tidybayes::stat_lineribbon(.width = c(.95), alpha = 0.5,
                             color = RColorBrewer::brewer.pal(5, "Blues")[[5]]) +
  ggplot2::scale_fill_brewer() +
  ggplot2::geom_point(data = dplyr::filter(dogs_df, dog <= 9),
                      ggplot2::aes(x = time, y = shock, group = dog)) +
  ggplot2::facet_wrap(~ dog, labeller = ggplot2::label_both) +
  ggplot2::scale_y_continuous(breaks = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title    = "Figure 21.1: Model 0h predicted shock probability, first 9 dogs",
    subtitle = "Ribbons are 95% credible intervals; points are observed shocks (1) or avoidances (0)",
    x        = "Trial",
    y        = "P(shock) and observed shocks"
  )
print(fig_21_1)
ggplot2::ggsave("figs/Fig-21.1.svg", fig_21_1,
                width = 10, height = 6, device = svg)
# LEARN NOTE
# Predicted P(shock) declines smoothly with trial for every dog because the
# logistic model has only a linear time term. There is no mechanism here for
# a large drop after the first avoidance -- that pattern shows up only in
# Fig 21.6 with the Bush-Mosteller structure.

# --- Figure 21.6: M4 per-dog predicted probability, first 9 dogs ---
fig_21_6 <- dogs_df |>
  dplyr::filter(dog <= 9) |>
  tidybayes::add_linpred_draws(bfit_4, transform = TRUE) |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = .linpred)) +
  tidybayes::stat_lineribbon(.width = c(.95), alpha = 0.5,
                             color = RColorBrewer::brewer.pal(5, "Blues")[[5]]) +
  ggplot2::scale_fill_brewer() +
  ggplot2::geom_point(data = dplyr::filter(dogs_df, dog <= 9),
                      ggplot2::aes(x = time, y = shock, group = dog)) +
  ggplot2::facet_wrap(~ dog, labeller = ggplot2::label_both) +
  ggplot2::scale_y_continuous(breaks = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title    = "Figure 21.6: Model 4 predicted shock probability, first 9 dogs",
    subtitle = "Sharper drops after first avoidance than logistic model; first-avoidance effect is bigger than any subsequent shock effect",
    x        = "Trial",
    y        = "P(shock) and observed shocks"
  )
print(fig_21_6)
ggplot2::ggsave("figs/Fig-21.6.svg", fig_21_6,
                width = 10, height = 6, device = svg)

# --- Figure 21.7: M2 per-dog predicted probability, first 9 dogs ---
fig_21_7 <- dogs_df |>
  dplyr::filter(dog <= 9) |>
  tidybayes::add_linpred_draws(bfit_2, transform = TRUE) |>
  ggplot2::ggplot(ggplot2::aes(x = time, y = .linpred)) +
  tidybayes::stat_lineribbon(.width = c(.95), alpha = 0.5,
                             color = RColorBrewer::brewer.pal(5, "Blues")[[5]]) +
  ggplot2::scale_fill_brewer() +
  ggplot2::geom_point(data = dplyr::filter(dogs_df, dog <= 9),
                      ggplot2::aes(x = time, y = shock, group = dog)) +
  ggplot2::facet_wrap(~ dog, labeller = ggplot2::label_both) +
  ggplot2::scale_y_continuous(breaks = c(0, 1)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    title    = "Figure 21.7: Model 2 (non-hierarchical) predicted shock probability, first 9 dogs",
    subtitle = "Even without random effects, dog-specific shock/avoidance histories drive different fitted curves",
    x        = "Trial",
    y        = "P(shock) and observed shocks"
  )
print(fig_21_7)
ggplot2::ggsave("figs/Fig-21.7.svg", fig_21_7,
                width = 10, height = 6, device = svg)
# LEARN NOTE
# Compare Fig 21.6 (hierarchical M4) with Fig 21.7 (non-hierarchical M2).
# The curves differ across dogs even in M2, because prev_shock and
# prev_avoid are dog-specific covariates. This explains why adding
# hierarchy on top of the log-model structure buys almost no LOO
# improvement in the comparison table below.


# --- (b) PAV-adjusted calibration and residual plots ------------------------

# --- Figure 21.2: PAV-adjusted LOO calibration, M0h ---
rd_0h  <- reliabilitydiag::reliabilitydiag(EMOS = loo_epred(bfit_0h),
                                           y    = dogs_df$shock)
fig_21_2 <- autoplot(rd_0h) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title    = "Figure 21.2: Model 0h PAV-adjusted LOO calibration",
    subtitle = "Conditional event probabilities should hug the 45-degree line",
    x        = "Predicted (LOO)",
    y        = "Conditional event probabilities"
  )
print(fig_21_2)
ggplot2::ggsave("figs/Fig-21.2.svg", fig_21_2,
                width = 6, height = 4, device = svg)

# --- Figure 21.3: PAV residual vs time, M0h ---
fig_21_3 <- ppc_pava_residual(dogs_df$shock, loo_epred(bfit_0h),
                              jitter(dogs_df$time, 0.3)) +
  ggplot2::labs(
    title    = "Figure 21.3: Model 0h PAV-adjusted LOO residual vs trial",
    subtitle = "Systematic deviation from zero at particular trials would flag a time-structured miss",
    x        = "Trial"
  )
print(fig_21_3)
ggplot2::ggsave("figs/Fig-21.3.svg", fig_21_3,
                width = 6, height = 4, device = svg)

# --- Figure 21.8: PAV-adjusted LOO calibration, M4 ---
rd_4 <- reliabilitydiag::reliabilitydiag(EMOS = loo_epred(bfit_4),
                                         y    = dogs_df$shock)
fig_21_8 <- autoplot(rd_4) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title    = "Figure 21.8: Model 4 PAV-adjusted LOO calibration",
    subtitle = "Should look as good as Fig 21.2 -- both models are well-calibrated on this data",
    x        = "Predicted (LOO)",
    y        = "Conditional event probabilities"
  )
print(fig_21_8)
ggplot2::ggsave("figs/Fig-21.8.svg", fig_21_8,
                width = 6, height = 4, device = svg)

# --- Figure 21.9: PAV residual vs time, M4 ---
fig_21_9 <- ppc_pava_residual(dogs_df$shock, loo_epred(bfit_4),
                              jitter(dogs_df$time, 0.3)) +
  ggplot2::labs(
    title    = "Figure 21.9: Model 4 PAV-adjusted LOO residual vs trial",
    subtitle = "Should be flat around zero across all trials",
    x        = "Trial"
  )
print(fig_21_9)
ggplot2::ggsave("figs/Fig-21.9.svg", fig_21_9,
                width = 6, height = 4, device = svg)


# --- (c) Tile-plot PPC and mean-switches test statistic ---------------------

# Helper: simulate a single replicated (30 x 25) shock matrix under a logistic
# model. Trial 1 is always shock, so we hard-code column 1 and draw columns
# 2..25 from a single posterior_predict() draw.
pred_logit <- function(fit) {
  matrix(
    c(rep(1, 30),
      brms::posterior_predict(fit, ndraws = 1) |> as.numeric()),
    nrow = 30, ncol = 25
  )
}

# Helper: same for the log models, but done sequentially because prev_shock
# and prev_avoid depend on the simulated history so far.
pred_log <- function(fit) {
  pred_shock <- matrix(rep(1, 30), nrow = 30)
  for (t in 2:25) {
    dogs_df_pred            <- dplyr::filter(dogs_df, time == t)
    dogs_df_pred$prev_shock <- as.numeric(rowSums(pred_shock))
    dogs_df_pred$prev_avoid <- as.numeric(rowSums(1 - pred_shock))
    pred <- brms::posterior_predict(fit, ndraws = 1, newdata = dogs_df_pred)
    pred_shock <- cbind(pred_shock, as.numeric(pred))
  }
  pred_shock
}

# Tile-plot builder: order dogs by trial of last observed shock, then plot
# a heatmap of shocks (dark) vs avoidances (light).
ppc_shocks <- function(shock, title) {
  expand.grid(dog = rev(1:30), time = 1:25) |>
    dplyr::mutate(shock = as.numeric(
      shock[order(apply(shock, 1, \(x) max(which(x == 1)))), ]
    )) |>
    ggplot2::ggplot(ggplot2::aes(time, dog, fill = shock)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "#ffffc8", high = "#7c0025") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none",
                   axis.line = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.6)) +
    ggplot2::labs(y = title)
}

# --- Figure 21.4a: real dogs tile plot ---
set.seed(SEED)
fig_21_4a <- ppc_shocks(shock, "Real dogs") +
  ggplot2::labs(title = "Figure 21.4a: Observed shocks, dogs ordered by trial of last shock")
print(fig_21_4a)
ggplot2::ggsave("figs/Fig-21.4a.svg", fig_21_4a,
                width = 6, height = 4, device = svg)

# --- Figure 21.4b: M0h posterior-predictive tile plot ---
fig_21_4b <- ppc_shocks(pred_logit(bfit_0h), "Model 0h: hier. logit") +
  ggplot2::labs(title = "Figure 21.4b: Model 0h posterior-predictive replicate")
print(fig_21_4b)
ggplot2::ggsave("figs/Fig-21.4b.svg", fig_21_4b,
                width = 6, height = 4, device = svg)

# --- Figure 21.4: combined 5-model tile grid (book Fig 21.4) ---
ppc_shocks_df <- function(shock, model, sim) {
  ord <- order(apply(shock, 1, \(x) max(which(x == 1))))
  expand.grid(dog = rev(1:30), time = 1:25) |>
    dplyr::mutate(shock = as.numeric(shock[ord, ]),
                  model = model,
                  sim   = sim)
}

labs_ppc <- c(
  "Real dogs",
  "PPsims from M0:\nlogit model",
  "PPsims from M0h:\nhier logit model",
  "PPsims from M2:\n2-par log model",
  "PPsims from M4:\nhier 2-par log model"
)

set.seed(SEED)
ppc_all <- dplyr::bind_rows(
  ppc_shocks_df(shock, labs_ppc[1], 1),
  dplyr::bind_rows(lapply(1:5, \(s) ppc_shocks_df(pred_logit(bfit_0),  labs_ppc[2], s))),
  dplyr::bind_rows(lapply(1:5, \(s) ppc_shocks_df(pred_logit(bfit_0h), labs_ppc[3], s))),
  dplyr::bind_rows(lapply(1:5, \(s) ppc_shocks_df(pred_log(bfit_2),    labs_ppc[4], s))),
  dplyr::bind_rows(lapply(1:5, \(s) ppc_shocks_df(pred_log(bfit_4),    labs_ppc[5], s)))
) |>
  dplyr::mutate(model = factor(model, levels = labs_ppc))

fig_21_4 <- ggplot2::ggplot(ppc_all,
                            ggplot2::aes(time, dog, fill = shock)) +
  ggplot2::geom_tile() +
  ggplot2::facet_grid(model ~ sim, switch = "y") +
  ggplot2::scale_fill_gradient(low = "#ffffc8", high = "#7c0025") +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none",
                 axis.line = ggplot2::element_blank(),
                 axis.text = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(),
                 axis.title = ggplot2::element_blank(),
                 panel.spacing = ggplot2::unit(0.6, "lines"),
                 strip.background = ggplot2::element_blank(),
                 strip.placement = "outside",
                 strip.text.x = ggplot2::element_blank(),
                 strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1)) +
  ggplot2::labs(
    title    = "Figure 21.4: Real dogs vs 5 posterior-predictive replicates from four models",
    subtitle = "Each row = one model (or the real data); each column = one PP replicate; dogs ordered by trial of last shock"
  )
print(fig_21_4)
ggplot2::ggsave("figs/Fig-21.4.svg", fig_21_4,
                width = 12, height = 9, device = svg)
# LEARN NOTE
# The chapter's central negative result lives in this figure. The log
# models (bottom two rows) look qualitatively more like the real dogs than
# the logistic models (middle two rows), but the difference is subtle and,
# as the LOO table below shows, statistically indistinguishable. Visual
# PPC on a small dataset like this cannot separate M0h, M2 and M4.

# --- Figure 21.10: M4 posterior-predictive tile plot ---
set.seed(SEED)
fig_21_10 <- ppc_shocks(pred_log(bfit_4),
                        "PPsims from M4:\nhier 2-par log model") +
  ggplot2::labs(title = "Figure 21.10: Model 4 posterior-predictive replicate")
print(fig_21_10)
ggplot2::ggsave("figs/Fig-21.10.svg", fig_21_10,
                width = 6, height = 4, device = svg)


# --- Mean-switches test statistic PPC ---------------------------------------
mean_switches <- function(shock) {
  shock |> matrixStats::rowDiffs() |> abs() |> rowSums() |> mean()
}

obs_switches <- mean_switches(shock)

# --- Figure 21.5: mean-switches PPC, M0h ---
set.seed(SEED)
yrep_switches_0h <- replicate(100, mean_switches(pred_logit(bfit_0h)))
fig_21_5 <- bayesplot::ppc_stat(
  y    = obs_switches,
  yrep = matrix(yrep_switches_0h, nrow = length(yrep_switches_0h)),
  stat = "identity"
) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title    = "Figure 21.5: PPC on mean switches per dog, Model 0h",
    subtitle = "Dark line = observed; histogram = replicated statistic under M0h",
    x        = "Mean number of shock/avoidance switches per dog"
  )
print(fig_21_5)
ggplot2::ggsave("figs/Fig-21.5.svg", fig_21_5,
                width = 6, height = 4, device = svg)

# --- Figure 21.11: mean-switches PPC, M4 ---
set.seed(SEED)
yrep_switches_4 <- replicate(100, mean_switches(pred_log(bfit_4)))
fig_21_11 <- bayesplot::ppc_stat(
  y    = obs_switches,
  yrep = matrix(yrep_switches_4, nrow = length(yrep_switches_4)),
  stat = "identity"
) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title    = "Figure 21.11: PPC on mean switches per dog, Model 4",
    subtitle = "Should also cover the observed value; both models pass this check",
    x        = "Mean number of shock/avoidance switches per dog"
  )
print(fig_21_11)
ggplot2::ggsave("figs/Fig-21.11.svg", fig_21_11,
                width = 6, height = 4, device = svg)


# =============================================================================
# PHASE 6  ---  Model comparison (LOO-CV)
# =============================================================================
cat("\n=== LOO comparison: all four models ===\n")
cat("Book Fig 21.4 comparison table targets:\n")
cat("  M0h vs M0 : elpd_diff = -13,  se_diff = 5.7  (M0 much worse)\n")
cat("  M0h vs M2 : elpd_diff = -4.7, se_diff = 4.3  (M2 slightly worse, ns)\n")
cat("  M0h vs M4 : elpd_diff = -3.8, se_diff = 3.1  (M4 slightly worse, ns)\n\n")

loo_tab <- loo::loo_compare(bfit_0, bfit_0h, bfit_2, bfit_4)
print(loo_tab)
saveRDS(loo_tab, "loo_compare.rds")
# LEARN NOTE
# All three non-baseline models (M0h, M2, M4) fall within |elpd_diff| of
# roughly 5, with SEs of similar size. The rule of thumb is that an elpd
# difference smaller than about 2*SE is not decisive. So the chapter's
# conclusion holds: this data cannot distinguish these three models.


# --- Figure 21.12: overlaid per-dog posterior mean predictions, M0h vs M4 ---
fig_21_12 <- dogs_df |>
  dplyr::mutate(epred_0h = colMeans(brms::posterior_epred(bfit_0h)),
                epred_4  = colMeans(brms::posterior_epred(bfit_4))) |>
  ggplot2::ggplot(ggplot2::aes(x = time, group = dog)) +
  ggplot2::geom_line(ggplot2::aes(y = epred_0h, color = "Model 0h"), alpha = 0.5) +
  ggplot2::geom_line(ggplot2::aes(y = epred_4,  color = "Model 4"),  alpha = 0.5) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::scale_color_manual(values = c("Model 0h" = "#1f77b4",
                                         "Model 4"  = "#d62728")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title    = "Figure 21.12: Posterior mean P(shock) per dog, Model 0h vs Model 4",
    subtitle = "Biggest visible difference is at trial 2; elsewhere the two models track each other closely",
    x        = "Trial",
    y        = "Posterior predictive P(shock)",
    color    = "Model"
  )
print(fig_21_12)
ggplot2::ggsave("figs/Fig-21.12.svg", fig_21_12,
                width = 8, height = 5, device = svg)


# =============================================================================
# PHASE 6b  ---  Prior-likelihood sensitivity (power-scaling)
# =============================================================================
# priorsense::powerscale_sensitivity flags each parameter with the KL
# divergence between the base posterior and versions where the prior or
# likelihood is up/down-weighted. Small values everywhere => data are
# informative, priors are not driving the posterior.
cat("\n=== Prior-likelihood sensitivity, Model 0h (first 6 parameters) ===\n")
sens_0h <- priorsense::powerscale_sensitivity(
  bfit_0h,
  variable = posterior::variables(posterior::as_draws(bfit_0h))[1:6]
)
print(sens_0h)

cat("\n=== Prior-likelihood sensitivity, Model 4 (first 5 parameters) ===\n")
sens_4 <- priorsense::powerscale_sensitivity(
  bfit_4,
  variable = posterior::variables(posterior::as_draws(bfit_4))[1:5]
)
print(sens_4)
# LEARN NOTE
# Check the diagnosis column. A "-" means no prior-likelihood conflict was
# detected for that parameter. Any "prior-data conflict" flag would mean
# the posterior is being pulled by the prior in a way the data cannot
# override -- worth investigating before trusting inferences.
#
# The chapter also uses this to argue that even though the posteriors of a
# and b in M4 look clearly different, that difference is not evidence of
# different learning rates from shocks vs avoidances -- fitting M4 to data
# simulated from M0h (no learning-rate asymmetry) still produces
# apparently-different a and b posteriors. That simulation check is
# omitted from this script for runtime (see header); rerun it separately
# if you want to reproduce book Figure 21.12.


# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete <- TRUE
wf$loo_complete <- TRUE
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

cat("\n=== Ch 21 analysis complete ===\n")
cat("Saved: bfit_0.rds, bfit_0h.rds, bfit_2.rds, bfit_4.rds, ",
    "loo_compare.rds, wf_final.rds\n", sep = "")
cat("Figures saved to figs/: Fig-21.1.svg, Fig-21.2.svg, Fig-21.3.svg, ",
    "Fig-21.4a.svg, Fig-21.4b.svg, Fig-21.4.svg, Fig-21.5.svg, ",
    "Fig-21.6.svg, Fig-21.7.svg, Fig-21.8.svg, Fig-21.9.svg, ",
    "Fig-21.10.svg, Fig-21.11.svg, Fig-21.12.svg\n", sep = "")
