# data/classification_ch20/classification_ch20_analysis.R
#
# Ch 20 — "Using a fitted model for decision analysis: Classification competition"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Fits a 3-component Gaussian mixture (means fixed at 0, -1, +1) to the
# per-series least-squares slope estimates from 1000 short time series, then
# uses the posterior classification probabilities p[n,k] to solve a decision
# problem: guess which series are trendless (mu=0), sloping-down (mu=-1), or
# sloping-up (mu=+1), and count the expected number of correct guesses under
# the contest's >=900 correct / $100k prize threshold.
#
# Deliberately omits:
#   Exercise 20.1 (single-stage integrated model over raw series),
#   Exercise 20.2 (per-series AR(1) errors),
#   Exercise 20.3 (unknown mu with ordered constraint, per-component sigma),
#   Phase 6 LOO comparison (chapter fits one model; nothing to compare),
#   Formal PPC with replicated datasets (chapter's "posterior prediction"
#     IS the classification probability p[n,k], not a replicated dataset).
# Includes one RBayesflow workflow figure not in the book (Fig 20.0, prior
# predictive check on slopes) because it is required by Phase 2 of the flow.
#
# Data:
#   Source: avehtari/Bayesian-Workflow GitHub repo, branch main,
#           timeseries/data/Series1000.txt (raw text, 135000 numbers).
#   Downloaded once and cached at ./Series1000.txt.
#   Reshape: matrix(scan(...), nrow = 1000, ncol = 135, byrow = TRUE).
#   Derived: slope, se from lm(y ~ time) per series, both x100 to give
#            per-century values matching the contest instructions.
#
# Models covered:
#   m1  fit_mix1  cmdstan_model("mixture.stan")    3-comp Gaussian mixture,
#                                                   mu fixed at (0, -1, +1),
#                                                   theta ~ implicit Dirichlet(1,1,1),
#                                                   sigma ~ Exponential(0.1).
#                                                   Purpose: recover theta, sigma.
#   m2  fit_mix2  cmdstan_model("mixture_2.stan")  Same sampling model, plus
#                                                   generated quantities matrix
#                                                   p[N,K] of posterior class
#                                                   probabilities per series.
#                                                   Purpose: enable decision analysis.
#
# Figures produced:
#   Fig 20.0    prior predictive check on slopes         [RBayesflow workflow fig]
#   Fig 20.1    raw 1000 time series overlaid            [book Fig 20.1]
#   Fig 20.2a   estimated slope vs SE, per series        [book Fig 20.2 (a)]
#   Fig 20.2b   histogram of 1000 estimated slopes       [book Fig 20.2 (b)]
#
# Book-target posterior summaries (success criteria for the fit):
#   theta[1] (mu =  0, null):    0.54, 95% central (0.50, 0.57)
#   theta[2] (mu = -1, down):    0.24, 95% central (0.21, 0.27)
#   theta[3] (mu = +1, up):      0.22, 95% central (0.20, 0.25)
#   sigma:                        0.40, 95% central (0.38, 0.44)
#   Rhat < 1.01 for all sampled parameters; bulk & tail ESS > 400.
#
# Book-target decision-analysis outputs:
#   table(choice):          ~(560, 231, 209)   # order 1=null, 2=down, 3=up
#   expected_correct:        854.0
#   sd_correct:              10.3
#   P(N_correct >= 900):     4.2e-6   ==>   ~1 in 240,000
#   P(N_correct >= 899.5):   5.3e-6   ==>   ~1 in 190,000  (continuity-corrected)
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/classification_ch20/.
# =============================================================================


# =============================================================================
# PHASE 0  ---  Environment setup
# =============================================================================

source("../../R/source_all.R")

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

options(brms.backend = "cmdstanr", mc.cores = 4)
SEED <- 1123

dir.create("figs", showWarnings = FALSE)


# diagnose_cmdstan() is sourced from R/diagnose_cmdstan.R via source_all.R.
# See that file for the full implementation and argument documentation.


# =============================================================================
# PHASE 1  ---  Data acquisition, exploration, and off-ramp assessment
# =============================================================================

wf <- init_workflow(mode = "learn", stage = "explore")

# -----------------------------------------------------------------------------
# Data: 1000 series of length 135, cached from the book's GitHub repo.
# -----------------------------------------------------------------------------
cache_path <- "Series1000.txt"
if (!file.exists(cache_path)) {
  src_url <- paste0(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/",
    "main/timeseries/data/Series1000.txt"
  )
  utils::download.file(src_url, cache_path, quiet = TRUE, mode = "wb")
  message("Downloaded time series data to: ", cache_path)
}

T <- 135L
N <- 1000L
series <- matrix(
  scan(cache_path, quiet = TRUE),
  nrow  = N,
  ncol  = T,
  byrow = TRUE
)
stopifnot(dim(series) == c(N, T))
cat(sprintf("Loaded %d series of length %d.\n", N, T))


# --- Figure 20.1: raw time series --------------------------------------------
# ggplot2 with 1000 * 135 = 135,000 line-segment points takes ~5-10 s to
# render; alpha 0.05 is chosen so the fanning-out pattern is visible without
# the overlaid lines going to solid black in dense regions.
series_df <- data.frame(
  series_id = rep(seq_len(N), each = T),
  time      = rep(seq_len(T), times = N),
  y         = as.vector(t(series))
)

fig_20_1 <- ggplot2::ggplot(
  series_df,
  ggplot2::aes(x = time, y = y, group = series_id)
) +
  ggplot2::geom_line(alpha = 0.05, linewidth = 0.3) +
  ggplot2::labs(
    title    = "Figure 20.1: Raw data from the time series problem",
    subtitle = "1000 series of length 135 fanning out from a common start; some are trendless, some carry a +/- 1 deg C per-century trend.",
    x        = "Time",
    y        = "y"
  ) +
  ggplot2::theme_minimal()
print(fig_20_1)
ggplot2::ggsave("figs/Fig-20.1.svg", fig_20_1,
                width = 6, height = 4, device = svg)
# LEARN NOTE
# The common starting point at t = 1 is a feature of how the series were
# generated. The visible spread at t = T is a mix of stochastic wander and
# any added trend. The eye cannot separate the two -- that is the whole
# problem this chapter is trying to solve.


# --- Per-series least-squares slope estimates --------------------------------
# Two-stage strategy from the book: reduce each length-135 series to a single
# slope estimate, then fit a mixture over the 1000 slopes. Not the most
# efficient use of the data (the series are autocorrelated, so LS is not
# minimum-variance), but the chapter's purpose is decision analysis, not
# efficient estimation. Ex 20.2 explores the autocorrelation-corrected version.
slope <- numeric(N)
se    <- numeric(N)
time  <- seq_len(T)
for (n in seq_len(N)) {
  fit_n    <- lm(series[n, ] ~ time)
  coefs    <- summary(fit_n)$coefficients
  slope[n] <- 100 * coefs[2, "Estimate"]     # per-century scale
  se[n]    <- 100 * coefs[2, "Std. Error"]   # per-century scale
}
slopes_df <- data.frame(slope = slope, se = se)


# --- Figure 20.2a: estimated slope vs SE, per series -------------------------
fig_20_2a <- ggplot2::ggplot(slopes_df, ggplot2::aes(x = slope, y = se)) +
  ggplot2::geom_point(alpha = 0.4, size = 0.7, colour = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1.05 * max(se)),
                              expand = c(0, 0)) +
  ggplot2::labs(
    title    = "Figure 20.2 (a): Least-squares estimates and standard errors of slopes",
    subtitle = "SEs are essentially constant across series -- they carry no information for classification.",
    x        = "Estimated slope (per century)",
    y        = "SE (per century)"
  ) +
  ggplot2::theme_minimal()
print(fig_20_2a)
ggplot2::ggsave("figs/Fig-20.2a.svg", fig_20_2a,
                width = 6, height = 4, device = svg)
# LEARN NOTE
# All 1000 SEs sit in a narrow band -- what changes across series is the
# slope, not its precision. This is why the mixture model treats each slope
# as a single scalar observation with a shared sigma, ignoring the individual
# SEs. If the SEs did vary meaningfully we would want a measurement-error
# mixture, which is Exercise 20.1 (integrate the two stages).


# --- Figure 20.2b: histogram of estimated slopes -----------------------------
fig_20_2b <- ggplot2::ggplot(slopes_df, ggplot2::aes(x = slope)) +
  ggplot2::geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
  ggplot2::labs(
    title    = "Figure 20.2 (b): Histogram of the 1000 estimated slopes",
    subtitle = "Three visible bumps near -1, 0, +1 motivate the 3-component Gaussian mixture with mu = (0, -1, 1) below.",
    x        = "Estimated slope (per century)",
    y        = "Frequency"
  ) +
  ggplot2::theme_minimal()
print(fig_20_2b)
ggplot2::ggsave("figs/Fig-20.2b.svg", fig_20_2b,
                width = 6, height = 4, device = svg)
# LEARN NOTE
# The three bumps at roughly -1, 0, +1 are the visual signature of the
# generating mixture. They are visible only because we already know from the
# contest description that these are the three component means. Without that
# hint the middle bump alone would look like a symmetric noise cloud.


# --- Off-ramp assessment -----------------------------------------------------
# assess_offramps() has no row for "mixture classification". Passing
# outcome_type = "continuous" yields the standard continuous menu (t-test,
# LM+bootstrap, full Stan). None of these actually solve the decision problem.
# The Bayesian path is chosen because we need P(component | slope) per series,
# not a single point estimate of a mean.
offramps <- assess_offramps(
  data         = slopes_df,
  outcome_var  = "slope",
  outcome_type = "continuous",
  goal         = "coefficient estimation"
)
print_offramps(offramps)
# LEARN NOTE
# The off-ramp menu here is a poor fit -- none of t-test, LM+bootstrap, or
# even "full Stan" as an unstructured regression addresses classification.
# The correct Bayesian off-ramp for this problem would be a frequentist
# mixture such as mixtools::normalmixEM. The current assess_offramps() does
# not include mixture alternatives; consider extending its decision matrix
# in a future ADR. For this chapter the Bayesian path is the intended demo
# so we proceed anyway.

wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 1,
  action    = "bayesian_selected",
  timestamp = Sys.time(),
  notes     = "Off-ramp taxonomy does not cover mixture classification; going Bayesian for full P(component | y_n)."
)))
wf$declared_goal  <- "coefficient estimation"
wf$n_observations <- N


# =============================================================================
# PHASE 2  ---  Prior specification and prior predictive check
# =============================================================================

# --- Priors (hand-written Stan, so declared in comments, not brms::prior) ----
#
# theta (mixture weights): implicit Dirichlet(1, 1, 1) via the
#   `simplex[K] theta` declaration. This is uniform over the 2-simplex --
#   no component is a priori preferred, and no component is ruled out.
#
# sigma (within-component slope SD): exponential(0.1).
#   rate = 0.1, mean = 1/0.1 = 10. Intentionally very diffuse on the
#   per-century slope scale. Slopes are numerically <~ 3 in absolute value
#   so a prior mean of 10 puts essentially no shape on the posterior -- the
#   data are expected to concentrate sigma tightly (posterior ~0.4).
#
# mu (component means): FIXED as data at (0, -1, +1) per the contest
#   description. Not sampled. Exercise 20.3 relaxes this to a sampled
#   ordered vector.
#
# NOTE on exponential(): in Stan's `exponential(0.1)` the argument is the
# RATE (not the mean). exponential(rate) has mean = 1/rate. This is a
# frequent source of student confusion.


# --- Prior predictive check on slopes ----------------------------------------
# For each of a small number of prior draws, simulate an entire vector of N
# slopes and overlay its density against the observed slopes. Uses
# bayesplot::ppc_dens_overlay so the panel matches the standard PPC visual
# even though what we are checking is prior-predictive, not posterior.
set.seed(SEED)
n_prior_draws <- 20L
yrep_prior    <- matrix(NA_real_, nrow = n_prior_draws, ncol = N)
mu_vec        <- c(0, -1, 1)
for (d in seq_len(n_prior_draws)) {
  raw     <- rexp(3L)
  theta_d <- raw / sum(raw)                    # Dirichlet(1,1,1) via norm-Exp
  sigma_d <- rexp(1L, rate = 0.1)              # Exponential(0.1)
  z_d     <- sample.int(3L, N, replace = TRUE, prob = theta_d)
  yrep_prior[d, ] <- rnorm(N, mean = mu_vec[z_d], sd = sigma_d)
}

# NOTE on xlim: no coord_cartesian clip. The prior-predictive draws span a
# much wider range than the observed slopes (that is the point of the check),
# so we let ggplot2 auto-scale to show the mismatch.
fig_20_0 <- bayesplot::ppc_dens_overlay(y = slope, yrep = yrep_prior) +
  ggplot2::labs(
    title    = "Figure 20.0: Prior predictive check on estimated slopes",
    subtitle = "Observed slopes (dark) vs 20 prior-predictive draws (light). Wide light lines reflect the diffuse Exp(0.1) prior on sigma (mean = 10).",
    x        = "Slope (per century)",
    y        = "Density"
  ) +
  ggplot2::theme_minimal()
print(fig_20_0)
ggplot2::ggsave("figs/Fig-20.0.svg", fig_20_0,
                width = 6, height = 4, device = svg)
# LEARN NOTE
# The prior-predictive spread is much wider than the observed slopes.
# That is fine -- with N = 1000 slopes the likelihood will overwhelm the
# diffuse Exp(0.1) prior on sigma. The check exists to confirm that no
# prior draw is qualitatively impossible (e.g. all draws clumped near zero),
# not to make the prior match the data.

wf$prior_pred_draws <- yrep_prior


# =============================================================================
# PHASE 3  ---  Model fitting
# =============================================================================

# --- Model 1: base mixture (fit_mix1) ----------------------------------------
# Purpose: recover theta and sigma; establish that convergence is clean before
# adding the generated quantities block. Model 2 below adds p[n,k] and refits.
# This matches the book's exposition sequence in section 20.3.

stan_mixture <- "
data {
  int<lower=1> K;
  int<lower=1> N;
  array[N] real y;
  array[K] real mu;
}
parameters {
  simplex[K] theta;
  real<lower=0> sigma;
}
model {
  array[K] real ps;
  sigma ~ exponential(0.1);           // rate = 0.1, mean = 10 (weak)
  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma);
    }
    target += log_sum_exp(ps);        // marginalise the latent z_n
  }
}
"
writeLines(stan_mixture, "mixture.stan")
mod_mix1 <- cmdstanr::cmdstan_model("mixture.stan")

data_mix <- list(K = 3L, N = N, y = slope, mu = mu_vec)

fit_mix1 <- mod_mix1$sample(
  data            = data_mix,
  seed            = SEED,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 0
)
saveRDS(fit_mix1, "fit_mix1.rds")

# Phase 4: diagnostics for m1
wf <- diagnose_cmdstan(fit_mix1, wf, params = c("theta", "sigma"))
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)

cat("\n=== Model 1: Posterior summary (compare to book Ch 20 sect 20.3) ===\n")
cat("Book: theta[1] 0.54, 95% (0.50, 0.57)\n")
cat("Book: theta[2] 0.24, 95% (0.21, 0.27)\n")
cat("Book: theta[3] 0.22, 95% (0.20, 0.25)\n")
cat("Book: sigma    0.40, 95% (0.38, 0.44)\n\n")
print(posterior::summarise_draws(
  fit_mix1$draws(variables = c("theta", "sigma")),
  posterior::default_summary_measures()
))


# --- Model 2: mixture with generated quantities for p[n, k] (fit_mix2) -------
# Purpose: same posterior over (theta, sigma) as m1, but each draw also
# produces the posterior classification probabilities p[n, k] for every one
# of the 1000 series. Posterior mean of p is what the decision analysis uses.
#
# LEARN NOTE
# We re-fit rather than editing the existing Stan file in place because
# generated quantities are recomputed only if the fit is re-run. In a real
# workflow you would put the GQ block into mixture.stan from the start and
# fit once; the two-fit structure here follows the book's exposition, which
# introduces generated quantities as a follow-on refinement.

stan_mixture_2 <- "
data {
  int<lower=1> K;
  int<lower=1> N;
  array[N] real y;
  array[K] real mu;
}
parameters {
  simplex[K] theta;
  real<lower=0> sigma;
}
model {
  array[K] real ps;
  sigma ~ exponential(0.1);
  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma);
    }
    target += log_sum_exp(ps);
  }
}
generated quantities {
  // p[n, k] = P(z_n = k | y_n, theta, sigma) at this posterior draw.
  matrix[N, K] p;
  for (n in 1:N) {
    vector[K] p_raw;
    for (k in 1:K) {
      p_raw[k] = theta[k] * exp(normal_lpdf(y[n] | mu[k], sigma));
    }
    for (k in 1:K) {
      p[n, k] = p_raw[k] / sum(p_raw);
    }
  }
}
"
writeLines(stan_mixture_2, "mixture_2.stan")
mod_mix2 <- cmdstanr::cmdstan_model("mixture_2.stan")

fit_mix2 <- mod_mix2$sample(
  data            = data_mix,
  seed            = SEED,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 0
)
saveRDS(fit_mix2, "fit_mix2.rds")

# Phase 4: diagnostics for m2 (sampled params only; skip p[.,.] GQ)
wf <- diagnose_cmdstan(fit_mix2, wf, params = c("theta", "sigma"))
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)

cat("\n=== Model 2 (with GQ): sanity check theta and sigma match m1 ===\n")
cat("Book: theta[1] 0.54  theta[2] 0.24  theta[3] 0.22  sigma 0.40\n\n")
print(posterior::summarise_draws(
  fit_mix2$draws(variables = c("theta", "sigma")),
  posterior::default_summary_measures()
))
# LEARN NOTE
# theta and sigma should match m1 within Monte Carlo error since the sampling
# model is identical. If they differ meaningfully, check whether the seed,
# iterations, or Stan program body drifted between the two files.


# =============================================================================
# PHASE 5  ---  Decision analysis
#
# The chapter's contribution: use the posterior classification probabilities
# to answer "how many series would we get right if we submitted this entry?"
# and "what is the chance we clear the 900-correct threshold?".
# =============================================================================

# --- Posterior mean classification matrix, N x K -----------------------------
# posterior::mean() on an rvar averages across posterior draws and returns
# a plain numeric array of the same [N, K] shape.
prob_sims <- posterior::as_draws_rvars(fit_mix2$draws(variables = "p"))
prob      <- mean(prob_sims$p)
stopifnot(is.numeric(prob), identical(dim(prob), c(N, 3L)))

cat("\n=== First 10 rows of posterior mean classification matrix ===\n")
cat("Book (approx, table on p.320):\n")
cat("  [,1]=null  [,2]=down  [,3]=up   -- row 1 ~ (0.09, 0.00, 0.91) etc.\n\n")
print(round(prob[1:10, ], 2))
# LEARN NOTE
# Each row sums to 1 and gives the posterior probability that series n
# belongs to the null, sloping-down, or sloping-up component. A row like
# (0.83, 0.17, 0.00) means "probably null, possibly sloping-down, almost
# certainly not sloping-up".


# --- Point classification: argmax over columns -------------------------------
max_prob <- apply(prob, 1, max)
choice   <- apply(prob, 1, which.max)

cat("\n=== Guess counts by component (compare to book) ===\n")
cat("Book: choice table ~ (560, 231, 209)   # 1=null, 2=down, 3=up\n")
print(table(choice))
# LEARN NOTE
# The guesses are not in the 500 / 250 / 250 proportion of the estimated
# theta -- we guess "null" a bit more often, because the mid-mass of the
# distribution is closer to mu = 0 than to the other two means. That is the
# decision-theoretic right answer for "maximise expected count correct":
# when uncertain, guess the most common category.


# --- Expected number correct and its SD -------------------------------------
# E[N_correct] = sum over n of max_p_n
# Var[N_correct] ~ sum over n of max_p_n * (1 - max_p_n)   (independent-Bernoulli approx)
expected_correct <- sum(max_prob)
sd_correct       <- sqrt(sum(max_prob * (1 - max_prob)))

cat("\n=== Expected number correct (compare to book) ===\n")
cat(sprintf("Book: expected_correct = 854.0   sd_correct = 10.3\n"))
cat(sprintf("Ours: expected_correct = %6.1f   sd_correct = %4.1f\n",
            expected_correct, sd_correct))
# LEARN NOTE
# 854 is the *best case* for the naive strategy of guessing the argmax
# category per series. The contest threshold is 900. The gap of ~46 is
# ~4.5 standard deviations away, so under this method the winning event is
# far in the tail.


# --- Probability of clearing the >=900 threshold ----------------------------
# Normal approximation to sum of independent Bernoulli(max_p_n) indicators.
p_ge_900       <- pnorm(expected_correct, mean = 900,   sd = sd_correct)
p_ge_899_5     <- pnorm(expected_correct, mean = 899.5, sd = sd_correct)
recip_900      <- 1 / p_ge_900
recip_899_5    <- 1 / p_ge_899_5

cat("\n=== Probability of winning the contest ===\n")
cat(sprintf("Book: P(N_correct >= 900)   ~ 4.2e-6  (1 in ~240,000)\n"))
cat(sprintf("Ours: P(N_correct >= 900)   = %.2e  (1 in ~%s)\n",
            p_ge_900,   format(round(recip_900), big.mark = ",")))
cat(sprintf("Book: P(N_correct >= 899.5) ~ 5.3e-6  (1 in ~190,000, continuity-corrected)\n"))
cat(sprintf("Ours: P(N_correct >= 899.5) = %.2e  (1 in ~%s)\n",
            p_ge_899_5, format(round(recip_899_5), big.mark = ",")))
# LEARN NOTE
# For the $10 entry / $100,000 prize to be a fair bet on expected monetary
# value alone, the win probability would have to be >= 1 in 10,000. It is
# ~20x too small for that. The chapter's whole point is that Bayesian
# workflow tells you *not to play* -- and the same machinery tells you
# exactly how far from playable the bet is.


# =============================================================================
# Save and wrap up
# =============================================================================
#
# ppc_complete: FALSE. The chapter has no traditional posterior predictive
#   check with replicated y_rep. The "posterior prediction" here is p[n, k],
#   which is used directly by Phase 5 decision analysis.
# loo_complete: FALSE. Only one model is fitted (m1 and m2 have the same
#   sampling posterior). No loo_compare() applies.

wf$ppc_complete <- FALSE
wf$loo_complete <- FALSE
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

cat("\n=== Ch 20 analysis complete ===\n")
cat("Saved: fit_mix1.rds, fit_mix2.rds, wf_final.rds, wf_context.json\n")
cat("       mixture.stan, mixture_2.stan\n")
cat("Cached data: Series1000.txt\n")
cat("Figures in figs/:\n")
cat("  Fig-20.0.svg    prior predictive check on slopes   [RBayesflow workflow]\n")
cat("  Fig-20.1.svg    raw 1000 time series overlaid      [book Fig 20.1]\n")
cat("  Fig-20.2a.svg   estimated slope vs SE, per series  [book Fig 20.2 (a)]\n")
cat("  Fig-20.2b.svg   histogram of estimated slopes      [book Fig 20.2 (b)]\n")
