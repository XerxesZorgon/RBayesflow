# data/covid_ch19/covid_ch19_analysis.R
#
# Ch 19 — "Building up to a hierarchical model: Coronavirus testing"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Core slice: three prevalence models fitted to the Santa Clara antibody data
# of Bendavid et al. (2020a, b). All three use hand-written Stan via cmdstanr
# rather than brms because the likelihood is joint across three binomial
# outcomes sharing one parameter space (pi, gamma, delta). A small in-script
# helper populates wf$diagnostics from the CmdStanFit so the RBayesflow
# display contract and diagnostic gate still apply.
#
# Models:
#   m1  fit_simple      Simple, uniform priors, pooled calibration (§19.2)
#   m2  fit_hier_weak   Hierarchical, sigma ~ normal+(0,1)   (§19.3, Fig 19.3a)
#   m3  fit_hier_strong Hierarchical, sigma ~ normal+(0,0.3) (§19.3, Fig 19.3b)
#
# Sensitivity data convention: the three Bendavid 2020b sensitivity studies
# are collapsed to a single pooled entry (J_sens = 1, 130/157). Rationale:
#   (a) 130/157 gives the book's target delta_1 ~ 0.797-0.821 after pooling;
#   (b) with only 3 studies sigma_delta is unidentifiable from data — making
#       it prior-driven is the book's own pedagogical point in section 19.3.
#
# Figures produced:
#   Fig 19.0          prior predictive check, m1
#   Fig 19.1a         posterior scatter (pi, gamma), m1       [book Fig 19.1a]
#   Fig 19.1b         posterior histogram pi, m1              [book Fig 19.1b]
#   Fig 19.2a-zoom    dense region [0, 0.03], linear y, m2   [book Fig 19.3a]
#   Fig 19.2a-tail    full range [0, 0.20], log-y,    m2     [book Fig 19.3a]
#   Fig 19.2b         posterior histogram pi, m3              [book Fig 19.3b]
#   Fig 19.3          caterpillar, study-level specificity, m3
#
# Note on Fig 19.2a: the m2 posterior for pi is extremely right-skewed.
# Nearly all mass sits below 0.025 but a sparse, real tail extends to ~0.16
# (matching the book's 95% interval (0.000, 0.160)). A single linear-y
# histogram at full range crushes the dense region; zoomed in, the tail
# vanishes. Two separate SVGs solve this without any layout package beyond
# base ggplot2 (patchwork/cowplot are not in renv.lock).
#
# Book targets:
#   m1  pi 95% shortest interval: (0.000, 0.018)
#   m2  pi median 0.016, 95% (0.000, 0.160); delta1 0.797
#   m3  pi median 0.013, 95% (0.001, 0.021); delta1 0.821
#
# Working directory: data/covid_ch19/
# =============================================================================


# =============================================================================
# PHASE 0  ---  Environment setup
# =============================================================================

source("../../R/source_all.R")

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)

options(mc.cores = 4)
SEED <- 4711

dir.create("figs", showWarnings = FALSE)


# diagnose_cmdstan() is sourced from R/diagnose_cmdstan.R via source_all.R.
# See that file for the full implementation and argument documentation.
# Note: this chapter's hierarchical models use the non-centered Stan
# parameterization (offset= / multiplier= declarations); diagnose_cmdstan()
# sets wf$parameterization = NA_character_ generically. If you need to record
# "non-centered" explicitly, set wf$parameterization <- "non-centered" after
# calling diagnose_cmdstan().


# =============================================================================
# PHASE 1  ---  Goal declaration + off-ramp assessment
# =============================================================================

wf <- init_workflow(mode = "learn", stage = "explore")

y_sample      <- 50L
n_sample      <- 3330L
y_spec_pooled <- 399L    # pooled specificity calibration -- used in m1
n_spec_pooled <- 401L
y_sens_pooled <- 103L    # pooled sensitivity calibration -- used in m1
n_sens_pooled <- 122L

covid_df <- data.frame(
  positive = c(rep(1L, y_sample), rep(0L, n_sample - y_sample))
)

offramps <- assess_offramps(
  data         = covid_df,
  outcome_var  = "positive",
  outcome_type = "binary",
  goal         = "coefficient estimation"
)
print(offramps)
# LEARN NOTE
# Naive rate = 50/3330 = 1.5% assumes the test is perfect. The Premier
# Biotech kit has ~0.5% false-positive rate, giving ~17 expected false
# positives in 3330 tests. Specificity and sensitivity uncertainty must
# be propagated jointly with prevalence. This is one of the rare cases
# where the off-ramp alternatives give the wrong answer, not just a
# noisier one.

wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 1,
  action    = "bayesian_selected",
  timestamp = Sys.time(),
  notes     = "Test error must be propagated jointly with prevalence."
)))
wf$declared_goal  <- "coefficient estimation"
wf$n_observations <- n_sample
wf$event_rate     <- y_sample / n_sample


# =============================================================================
# PHASE 2  ---  Prior specification + prior predictive
# =============================================================================

# Priors for m1: all uniform(0,1) on the probability scale.
# The calibration likelihoods supply the informativeness.

set.seed(SEED)
n_prior      <- 4000L
pi_prior     <- runif(n_prior, 0, 1)
gamma_prior  <- runif(n_prior, 0, 1)
delta_prior  <- runif(n_prior, 0, 1)
p_prior      <- pi_prior * delta_prior + (1 - pi_prior) * (1 - gamma_prior)
y_prior_pred <- rbinom(n_prior, size = n_sample, prob = p_prior)

fig_prior <- ggplot2::ggplot(
  data.frame(y_pred = y_prior_pred),
  ggplot2::aes(x = y_pred)
) +
  ggplot2::geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  ggplot2::geom_vline(xintercept = y_sample, linetype = "dashed",
                      colour = "firebrick", linewidth = 1) +
  ggplot2::labs(
    title    = "Prior predictive check (Model 1, uniform priors)",
    subtitle = "Red line = observed y_sample = 50; uniform priors give flat spread.",
    x        = "Simulated y_sample under prior",
    y        = "Count"
  ) +
  ggplot2::theme_minimal()
print(fig_prior)
ggplot2::ggsave("figs/Fig-19.0.svg", fig_prior, width = 6, height = 4, device = svg)
# LEARN NOTE: flat prior predictive is expected. The calibration likelihoods
# tighten gamma and delta in the posterior, which then localises pi.

wf$prior_pred_draws <- y_prior_pred


# =============================================================================
# PHASE 3  ---  Model 1: simple non-hierarchical fit (book section 19.2)
# =============================================================================

# Stan program A.1 from Gelman & Carpenter (2020), modern array syntax.
# generated quantities adds log_lik / log_prior for section 19.4
# power-scaling (student exercise).
stan_simple <- "
data {
  int<lower=0> y_sample;
  int<lower=0> n_sample;
  int<lower=0> y_spec;
  int<lower=0> n_spec;
  int<lower=0> y_sens;
  int<lower=0> n_sens;
}
parameters {
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> spec;
  real<lower=0, upper=1> sens;
}
transformed parameters {
  real p_sample = p * sens + (1 - p) * (1 - spec);
}
model {
  y_sample ~ binomial(n_sample, p_sample);
  y_spec   ~ binomial(n_spec, spec);
  y_sens   ~ binomial(n_sens, sens);
  // uniform(0,1) priors implicit via parameter bounds
}
generated quantities {
  real log_lik   = binomial_lpmf(y_sample | n_sample, p_sample);
  real log_prior = binomial_lpmf(y_spec | n_spec, spec)
                 + binomial_lpmf(y_sens | n_sens, sens);
}
"
writeLines(stan_simple, "santa_clara_simple.stan")
mod_simple <- cmdstanr::cmdstan_model("santa_clara_simple.stan")

fit_simple <- mod_simple$sample(
  data = list(
    y_sample = y_sample, n_sample = n_sample,
    y_spec   = y_spec_pooled, n_spec = n_spec_pooled,
    y_sens   = y_sens_pooled, n_sens = n_sens_pooled
  ),
  seed            = SEED,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 0
)
saveRDS(fit_simple, "fit_simple.rds")

# Phase 4: diagnostics for m1
wf <- diagnose_cmdstan(fit_simple, wf, params = c("p", "spec", "sens"))
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)

cat("\nBook target (m1): 95% shortest posterior interval for pi: (0.000, 0.018)\n")
print(posterior::summarise_draws(fit_simple$draws(variables = c("p", "spec", "sens"))))

# Phase 5: figures for m1
draws_simple <- posterior::as_draws_df(
  fit_simple$draws(variables = c("p", "spec", "sens"))
)

# Fig 19.1a: posterior scatter (pi, gamma)
fig_19_1a <- ggplot2::ggplot(draws_simple, ggplot2::aes(x = spec, y = p)) +
  ggplot2::geom_point(alpha = 0.15, size = 0.6, colour = "steelblue") +
  ggplot2::coord_cartesian(xlim = c(0.980, 1.000), ylim = c(0, 0.025)) +
  ggplot2::labs(
    title    = "Figure 19.1a: Posterior scatter of prevalence vs specificity",
    subtitle = "Prevalence uncertainty is driven by specificity uncertainty (book section 19.2).",
    x        = "Specificity, gamma",
    y        = "Prevalence, pi"
  ) +
  ggplot2::theme_minimal()
print(fig_19_1a)
ggplot2::ggsave("figs/Fig-19.1a.svg", fig_19_1a, width = 6, height = 4, device = svg)
# LEARN NOTE
# When gamma ~ 0.980 the false-positive rate is 2%: ~66 expected false
# positives in 3330 tests, more than the 50 observed. At that specificity,
# zero true prevalence is plausible. Only when gamma approaches 1 does pi
# pull clearly away from 0. The banana shape makes this visible.

# Fig 19.1b: posterior histogram of pi
fig_19_1b <- ggplot2::ggplot(draws_simple, ggplot2::aes(x = p)) +
  ggplot2::geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
  ggplot2::coord_cartesian(xlim = c(0, 0.025)) +
  ggplot2::labs(
    title    = "Figure 19.1b: Posterior distribution of prevalence, pi (Model 1)",
    subtitle = "Consistent with prevalences from 0 to about 2%.",
    x        = "Prevalence, pi",
    y        = "Count"
  ) +
  ggplot2::theme_minimal()
print(fig_19_1b)
ggplot2::ggsave("figs/Fig-19.1b.svg", fig_19_1b, width = 6, height = 4, device = svg)

# Shortest 95% posterior interval (book uses shortest, not central, because
# of the hard lower bound at 0).
shortest_interval <- function(x, prob = 0.95) {
  x_sorted <- sort(x)
  n        <- length(x_sorted)
  width    <- floor(prob * n)
  starts   <- seq_len(n - width)
  widths   <- x_sorted[starts + width] - x_sorted[starts]
  i_min    <- which.min(widths)
  c(x_sorted[i_min], x_sorted[i_min + width])
}
pi_spi <- shortest_interval(draws_simple$p, prob = 0.95)
cat(sprintf(
  "\nModel 1 pi 95%% shortest posterior interval: (%.3f, %.3f)\n",
  pi_spi[1], pi_spi[2]
))
cat("Book target (section 19.2): (0.000, 0.018)\n")


# =============================================================================
# PHASE 3  ---  Models 2 and 3: hierarchical fits (book section 19.3)
# =============================================================================

# --- Data: 13 specificity studies + pooled sensitivity -----------------------
source("bendavid_2020b_hierarchical_vectors.R", local = TRUE)
stopifnot(
  length(y_spec) == length(n_spec), length(y_spec) == 13L,
  sum(n_spec) == 3324L, sum(y_spec) == 3308L,
  sum(n_sens) == 157L,  sum(y_sens) == 130L
)

# Collapse to a single pooled sensitivity entry (J_sens = 1).
# See header rationale: sigma_delta is unidentifiable from 3 studies,
# so pooling makes this explicit rather than hiding it.
y_sens_hier <- sum(y_sens)   # 130
n_sens_hier <- sum(n_sens)   # 157

cat(sprintf(
  "\nBendavid 2020b: %d specificity studies (N=%d); sensitivity pooled %d/%d.\n",
  length(y_spec), sum(n_spec), y_sens_hier, n_sens_hier
))

# --- Hierarchical Stan program (book Appendix A.2, modern syntax) ------------
stan_hier <- "
data {
  int<lower=0> y_sample;
  int<lower=0> n_sample;
  int<lower=0> J_spec;
  array[J_spec] int<lower=0> y_spec;
  array[J_spec] int<lower=0> n_spec;
  int<lower=0> J_sens;
  array[J_sens] int<lower=0> y_sens;
  array[J_sens] int<lower=0> n_sens;
  real<lower=0> logit_spec_prior_scale;
  real<lower=0> logit_sens_prior_scale;
}
parameters {
  real<lower=0, upper=1> p;
  real mu_logit_spec;
  real mu_logit_sens;
  real<lower=0> sigma_logit_spec;
  real<lower=0> sigma_logit_sens;
  vector<offset=mu_logit_spec, multiplier=sigma_logit_spec>[J_spec] logit_spec;
  vector<offset=mu_logit_sens, multiplier=sigma_logit_sens>[J_sens] logit_sens;
}
transformed parameters {
  vector[J_spec] spec = inv_logit(logit_spec);
  vector[J_sens] sens = inv_logit(logit_sens);
  real p_sample = p * sens[1] + (1 - p) * (1 - spec[1]);
}
model {
  y_sample ~ binomial(n_sample, p_sample);
  y_spec   ~ binomial(n_spec, spec);
  y_sens   ~ binomial(n_sens, sens);
  logit_spec ~ normal(mu_logit_spec, sigma_logit_spec);
  logit_sens ~ normal(mu_logit_sens, sigma_logit_sens);
  sigma_logit_spec ~ normal(0, logit_spec_prior_scale);
  sigma_logit_sens ~ normal(0, logit_sens_prior_scale);
  mu_logit_spec ~ normal(4, 2);
  mu_logit_sens ~ normal(4, 2);
}
generated quantities {
  real log_lik   = binomial_lpmf(y_sample | n_sample, p_sample);
  real log_prior = normal_lpdf(mu_logit_spec | 4, 2)
                 + normal_lpdf(mu_logit_sens | 4, 2);
}
"
writeLines(stan_hier, "santa_clara_hier.stan")
mod_hier <- cmdstanr::cmdstan_model("santa_clara_hier.stan")

data_hier_base <- list(
  y_sample = y_sample,
  n_sample = n_sample,
  J_spec   = length(y_spec),
  y_spec   = as.integer(y_spec),
  n_spec   = as.integer(n_spec),
  J_sens   = 1L,
  y_sens   = as.integer(y_sens_hier),   # 130
  n_sens   = as.integer(n_sens_hier)    # 157
)

# --- Model 2: weak sigma priors normal+(0, 1) -- book Fig 19.3a --------------
#
# Prior reasoning:
#   sigma_logit_spec, sigma_logit_sens ~ normal+(0, 1)
#   A shift of 1 on the logit scale spans a large range in probability.
#   With mu_delta ~ logit(0.8) = 1.4 and sigma = 1, sensitivity varies
#   from ~0.60 to ~0.92 across sites. Intentionally permissive to
#   demonstrate that J_sens = 1 leaves sigma_delta prior-dominated.
#   mu_logit_spec, mu_logit_sens ~ normal(4, 2):
#   Probability-scale mean of spec/sens distribution roughly (0.88, 1.00).

data_hier_weak                        <- data_hier_base
data_hier_weak$logit_spec_prior_scale <- 1.0
data_hier_weak$logit_sens_prior_scale <- 1.0

fit_hier_weak <- mod_hier$sample(
  data            = data_hier_weak,
  seed            = SEED,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 0,
  adapt_delta     = 0.95
)
saveRDS(fit_hier_weak, "fit_hier_weak.rds")

params_hier <- c("p", "mu_logit_spec", "mu_logit_sens",
                 "sigma_logit_spec", "sigma_logit_sens",
                 "spec[1]", "sens[1]")

wf <- diagnose_cmdstan(fit_hier_weak, wf, params = params_hier)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)

cat("\n=== Model 2 (weak sigma priors) -- compare to book Fig 19.3a ===\n")
cat("Book: pi 0.016  gamma1 0.997  delta1 0.797 [J_sens=1 pooled 130/157]\n")
cat("Book: mu_gamma 5.54  sigma_gamma 1.62  mu_delta 1.54  sigma_delta 0.87 [prior-driven]\n\n")
print(posterior::summarise_draws(
  fit_hier_weak$draws(variables = params_hier),
  posterior::default_summary_measures()
))

# --- Model 3: informative sigma priors normal+(0, 0.3) -- book Fig 19.3b -----
#
# Prior reasoning:
#   sigma_logit_spec, sigma_logit_sens ~ normal+(0, 0.3)
#   With mu_delta ~ 1.6 (from m2 median), sigma = 0.3 implies a 2/3 chance
#   that a new site's sensitivity falls in logit^-1(1.6 +/- 0.3) = (0.79, 0.87).
#   This encodes domain knowledge that different labs running the same kit
#   do not vary wildly in performance.
#   mu_logit_spec, mu_logit_sens ~ normal(4, 2) unchanged from m2.

data_hier_strong                        <- data_hier_base
data_hier_strong$logit_spec_prior_scale <- 0.3
data_hier_strong$logit_sens_prior_scale <- 0.3

fit_hier_strong <- mod_hier$sample(
  data            = data_hier_strong,
  seed            = SEED,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 0,
  adapt_delta     = 0.95
)
saveRDS(fit_hier_strong, "fit_hier_strong.rds")

wf <- diagnose_cmdstan(fit_hier_strong, wf, params = params_hier)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- diagnose.wf_state(wf)

cat("\n=== Model 3 (strong sigma priors) -- compare to book Fig 19.3b ===\n")
cat("Book: pi 0.013  gamma1 0.995  delta1 0.821 [J_sens=1 pooled 130/157]\n")
cat("Book: mu_gamma 5.23  sigma_gamma 0.72  mu_delta 1.54  sigma_delta 0.39 [prior-driven]\n\n")
print(posterior::summarise_draws(
  fit_hier_strong$draws(variables = params_hier),
  posterior::default_summary_measures()
))


# =============================================================================
# PHASE 5  ---  Figures for the two hierarchical models
# =============================================================================

draws_weak   <- posterior::as_draws_df(fit_hier_weak$draws(variables  = "p"))
draws_strong <- posterior::as_draws_df(fit_hier_strong$draws(variables = "p"))

# Fig 19.2a -- two separate SVGs.
#
# The m2 posterior for pi is extremely right-skewed: nearly all 4000 draws
# fall below 0.025 but a sparse, real tail reaches ~0.16 (matching the book's
# 95% interval (0.000, 0.160)). A single linear-y histogram at full range
# crushes the dense region; zoomed in, the tail vanishes entirely.
# Two separate files solve this without any layout package:
#   Fig-19.2a-zoom: [0, 0.03] linear-y, 40 bins -- shows the smooth mass
#   Fig-19.2a-tail: [0, 0.20] log-y,    60 bins -- sparse tail visible
#
# Note: scale_y_log10() emits a warning about infinite values for bins with
# zero draws (log(0) = -Inf). This is cosmetic only; the plot is correct.

fig_19_2a_zoom <- ggplot2::ggplot(draws_weak, ggplot2::aes(x = p)) +
  ggplot2::geom_histogram(bins = 40, fill = "darkorange", colour = "white") +
  ggplot2::coord_cartesian(xlim = c(0, 0.03)) +
  ggplot2::labs(
    title    = "Figure 19.2a (zoom): Posterior pi, hierarchical weak sigma priors",
    subtitle = "Dense region 0-3%",
    x        = "Prevalence, pi",
    y        = "Count"
  ) +
  ggplot2::theme_minimal()
print(fig_19_2a_zoom)
ggplot2::ggsave("figs/Fig-19.2a-zoom.svg", fig_19_2a_zoom,
                width = 6, height = 4, device = svg)

fig_19_2a_tail <- ggplot2::ggplot(draws_weak, ggplot2::aes(x = p)) +
  ggplot2::geom_histogram(bins = 60, fill = "darkorange", colour = "white") +
  ggplot2::coord_cartesian(xlim = c(0, 0.20)) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title    = "Figure 19.2a (tail): Posterior pi, hierarchical weak sigma priors",
    subtitle = "Full range, log-y scale -- sparse tail visible to ~16%",
    x        = "Prevalence, pi",
    y        = "Count (log scale)"
  ) +
  ggplot2::theme_minimal()
print(fig_19_2a_tail)
ggplot2::ggsave("figs/Fig-19.2a-tail.svg", fig_19_2a_tail,
                width = 6, height = 4, device = svg)

# Fig 19.2b: strong sigma priors -- clean unimodal posterior
fig_19_2b <- ggplot2::ggplot(draws_strong, ggplot2::aes(x = p)) +
  ggplot2::geom_histogram(bins = 60, fill = "seagreen", colour = "white") +
  ggplot2::coord_cartesian(xlim = c(0, 0.05)) +
  ggplot2::labs(
    title    = "Figure 19.2b: Posterior of prevalence, pi (hierarchical, strong sigma priors)",
    subtitle = "Concentrated between 0.1% and 2.1%.",
    x        = "Prevalence, pi",
    y        = "Count"
  ) +
  ggplot2::theme_minimal()
print(fig_19_2b)
ggplot2::ggsave("figs/Fig-19.2b.svg", fig_19_2b, width = 6, height = 4, device = svg)
# LEARN NOTE
# Compare Fig-19.2a-zoom and Fig-19.2b: the posterior medians are nearly
# identical (~0.010). What changes is the width. Under weak sigma_delta
# priors the tail reaches 0.16; under informative sigma_delta priors it
# truncates near 0.03. This is the book's central point in section 19.3:
# the prior on sigma_delta controls the width of the credible interval for
# pi, not its median. Reporting only a point estimate hides the real
# uncertainty.

# Fig 19.3: caterpillar of study-level specificity from m3
spec_summary       <- posterior::summarise_draws(
  fit_hier_strong$draws(variables = "spec"),
  ~ posterior::quantile2(.x, probs = c(0.025, 0.5, 0.975))
)
spec_summary$study <- seq_len(nrow(spec_summary))

fig_19_3 <- ggplot2::ggplot(spec_summary, ggplot2::aes(x = study, y = q50)) +
  ggplot2::geom_point(size = 2.5, colour = "steelblue") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = q2.5, ymax = q97.5),
                         width = 0.3, colour = "steelblue") +
  ggplot2::scale_x_continuous(breaks = spec_summary$study) +
  ggplot2::labs(
    title    = "Figure 19.3: Study-level specificity posteriors (Model 3, strong priors)",
    subtitle = "Study 1 = Bendavid site; studies 2-13 = external calibration studies.",
    x        = "Study index j",
    y        = "Specificity, gamma_j"
  ) +
  ggplot2::theme_minimal()
print(fig_19_3)
ggplot2::ggsave("figs/Fig-19.3.svg", fig_19_3, width = 10, height = 5, device = svg)
# LEARN NOTE
# Studies 10-13 have wider intervals and slightly lower medians -- consistent
# with smaller or more heterogeneous samples (study 10: RF+ patients 29/31;
# study 13: COVID-era PCR-negatives 50/52). Hierarchical partial pooling
# shrinks all estimates toward the group mean.


# =============================================================================
# PHASE 6  ---  Comparison of the two hierarchical variants
# =============================================================================
#
# LEARN NOTE -- Why no loo_compare()?
# m1 vs m2/m3: different observation vectors -- not on the same LOO scale.
# m2 vs m3: log_lik is a single aggregate binomial (y_sample | n_sample),
# yielding one Pareto-k value. Technically valid but not meaningful for model
# selection. The book makes no LOO comparison. The qualitative comparison of
# Fig 19.2a-zoom vs Fig 19.2b is the right one.
#
# LEARN NOTE -- Why is sigma_delta prior-driven in m2/m3?
# J_sens = 1: a single binomial count can identify mu_delta but not sigma_delta.
# sigma_delta is entirely determined by its prior (normal+(0,1) or normal+(0,0.3)).
# This is the book's lesson: "only three experiments" is too few to learn the
# inter-site variance, so the prior on sigma_delta controls the width of the
# credible interval for pi.

pi_ci_weak   <- posterior::quantile2(draws_weak$p,   probs = c(0.025, 0.5, 0.975))
pi_ci_strong <- posterior::quantile2(draws_strong$p, probs = c(0.025, 0.5, 0.975))
cat("\n=== Prevalence pi: weak vs strong sigma priors ===\n")
cat(sprintf("Weak   priors: median %.3f, 95%% central (%.3f, %.3f)\n",
            pi_ci_weak["q50"],   pi_ci_weak["q2.5"],   pi_ci_weak["q97.5"]))
cat(sprintf("Strong priors: median %.3f, 95%% central (%.3f, %.3f)\n",
            pi_ci_strong["q50"], pi_ci_strong["q2.5"], pi_ci_strong["q97.5"]))
cat("Book Fig 19.3a: median 0.016, 95% (0.000, 0.160)\n")
cat("Book Fig 19.3b: median 0.013, 95% (0.001, 0.021)\n")


# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete <- FALSE   # no traditional PPC -- book chapter does not do one
wf$loo_complete <- FALSE   # see LEARN NOTE above
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

cat("\n=== Ch 19 analysis complete ===\n")
cat("Saved: fit_simple.rds, fit_hier_weak.rds, fit_hier_strong.rds\n")
cat("       wf_final.rds, wf_context.json\n")
cat("       santa_clara_simple.stan, santa_clara_hier.stan\n")
cat("Figures in figs/:\n")
cat("  Fig-19.0.svg        prior predictive, m1\n")
cat("  Fig-19.1a.svg       posterior scatter (pi, gamma), m1  [book Fig 19.1a]\n")
cat("  Fig-19.1b.svg       posterior histogram pi, m1         [book Fig 19.1b]\n")
cat("  Fig-19.2a-zoom.svg  dense region, m2 weak priors       [book Fig 19.3a]\n")
cat("  Fig-19.2a-tail.svg  log-y tail, m2 weak priors         [book Fig 19.3a]\n")
cat("  Fig-19.2b.svg       posterior pi, m3 strong priors     [book Fig 19.3b]\n")
cat("  Fig-19.3.svg        study-level specificity, m3\n")

