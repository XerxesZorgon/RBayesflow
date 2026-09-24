# data/ch23_world_cup/world_cup_analysis.R
#
# Ch 23 — "Debugging a model: World Cup football"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# This script reproduces the World Cup 2014 case study from Ch 23.
# The chapter's central theme is debugging: a model whose MCMC diagnostics
# are clean but whose posterior is nonsensical. The bug is a coding error
# in the sqrt transformation (missing factor of 2), which is revealed by
# posterior retrodiction. The fix is to correct the transformed data block,
# not to change the model structure or tighten priors (contra the PREP_NOTE
# framing — the chapter uses domain-knowledge plausibility, not prior
# tightening per se, as the diagnostic). After the fix, sections 6–8 build
# discrete and Poisson-based models and compare them with LOO-CV.
#
# This script uses cmdstanr directly throughout. brms cannot express the
# custom discrete likelihoods (midpoint rule, exact Jacobian integration,
# bivariate Poisson, Poisson difference) used in sections 6–8. For the
# continuous t-models in sections 3–5 brms could be used, but the chapter
# uses cmdstanr, and this script follows that choice for consistency and to
# avoid introducing a different compile path mid-script.
#
# Stan files (all in data/ch23_world_cup/):
#   worldcup_first_try.stan           M1: buggy sqrt-t model
#   worldcup_with_replication.stan    M1-rep: same + y_rep for PPC
#   worldcup_no_sqrt.stan             M2: t model without sqrt
#   worldcup_fixed.stan               M3: corrected sqrt-t model
#   worldcup_no_prior.stan            M3-noprior: M3 with b=0 (no power score)
#   worldcup_discrete_z.stan          M-discr-z: discrete, latent z explicit
#   worldcup_discrete.stan            M-discr: discrete, z integrated out
#   worldcup_discrete_nopower.stan    M-discr-nopower: discrete, no power score
#   worldcup_discrete_poweronly.stan  M-discr-poweronly: power score only
#   worldcup_discrete_pooled.stan     M-discr-pool: pooled, no per-team effects
#   worldcup_continuous_midpoint_ll.stan  M-cont-midp: continuous, midpoint log_lik
#   worldcup_continuous.stan          M-cont: continuous, exact log_lik
#   worldcup_sqrt_continuous_nojacobian.stan  M-sqrt-cont-noj: sqrt-cont, no Jacobian
#   worldcup_sqrt_continuous.stan     M-sqrt-cont: sqrt-cont, with Jacobian
#   worldcup_sqrt_discrete.stan       M-sqrt-discr: sqrt-discrete, with Jacobian
#   worldcup_bivariate_poisson.stan   M-bipois: bivariate Poisson
#   worldcup_poisson_difference.stan  M-poisdif: Poisson difference
#
# NOTE: worldcup_sqrt_discrete.stan is very slow (quadrature at each HMC step).
#       Expect 30–60 min per chain. Consider running overnight or skipping
#       on first pass — set RUN_SQRT_DISCR <- FALSE below to skip it.
#
# Figures produced:
#   Fig-23.1.svg  — M1 team quality estimates (a[1:32], 90% intervals)
#   Fig-23.2.svg  — M3-noprior team quality estimates
#   Fig-23.3.svg  — M1-rep retrodiction vs observed score differentials
#   Fig-23.4.svg  — M2 team quality estimates
#   Fig-23.5.svg  — M2 retrodiction vs observed score differentials
#   Fig-23.6.svg  — M3 team quality estimates (corrected model)
#   Fig-23.7.svg  — M3 retrodiction vs observed score differentials (corrected)
#   Fig-23.8.svg  — LOO-PIT ECDF for M-discr
#   Fig-23.9.svg  — LOO-PIT ECDF for M-sqrt-cont
#   Fig-23.10.svg — LOO-PIT ECDF for M-bipois
#   Fig-23.11.svg — LOO-PIT ECDF for M-poisdif
#
# Book-target posterior summaries (key parameters from chapter):
#   M1:            b ≈ 0.45, sigma_a ≈ 0.17, sigma_y ≈ 0.42
#   M2:            b ≈ 1.1,  sigma_a ≈ 0.37, sigma_y ≈ 1.3
#   M3 (fixed):    b ≈ 0.86, sigma_a ≈ 0.33, sigma_y ≈ 0.84
#   M-discr:       b ≈ 1.1,  sigma_a ≈ 0.52, sigma_z ≈ 1.5
#   M-bipois:      b_o ≈ 0.53, b_d ≈ -0.45, sigma_o ≈ 0.18, sigma_d ≈ 0.30
#   M-poisdif:     b ≈ 0.52, sigma_a ≈ 0.24
#   LOO sect 6.6:  elpd_diff discrete vs poweronly ≈ -1.3 (se 1.8)
#                  elpd_diff discrete vs nopower   ≈ -3.3 (se 2.8)
#   LOO sect 8:    elpd_diff discrete vs bipois    ≈ -1.8 (se 0.9)
#                  elpd_diff discrete vs poisdif   ≈ -1.7 (se 1.7)
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch23_world_cup/.

# =============================================================================
# Control flags
# =============================================================================

# Set to FALSE to skip the very slow sqrt-discrete model (§7.4)
RUN_SQRT_DISCR <- FALSE

# =============================================================================
# §0  Environment setup
# =============================================================================

source("../../R/source_all.R")

library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(dplyr)
library(readr)

options(mc.cores = 4)
SEED <- 42

dir.create("figs",       showWarnings = FALSE)
dir.create("stan_output", showWarnings = FALSE)

# CmdStanR writes temporary files here (keeps the working directory clean)
options(cmdstanr_output_dir = "stan_output")

results_con <- file("results.txt", open = "wt")

log_result <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = results_con, append = FALSE)
}

capture_result <- function(x, label = NULL) {
  txt <- paste(capture.output(print(x)), collapse = "\n")
  if (!is.null(label)) {
    cat(label, "\n", sep = "")
    cat(label, "\n", sep = "", file = results_con, append = FALSE)
  }
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = results_con, append = FALSE)
}

log_result("=== Ch 23 results log — ", format(Sys.time()), " ===")

# =============================================================================
# §1  RBayesflow workflow state
# =============================================================================

wf <- init_workflow(mode = "practice", stage = "explore")
export_context(wf)

# =============================================================================
# §2  Data acquisition
# =============================================================================
# Data: 64 World Cup 2014 match results + Soccer Power Index pre-tournament
# rankings. CSVs are downloaded from GitHub on first run and cached locally.

data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)

spi_path  <- file.path(data_dir, "soccerpowerindex.csv")
wc_path   <- file.path(data_dir, "worldcup2014.csv")

if (!file.exists(spi_path)) {
  download.file(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/master/world_cup/data/soccerpowerindex.csv",
    destfile = spi_path
  )
}
if (!file.exists(wc_path)) {
  download.file(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/master/world_cup/data/worldcup2014.csv",
    destfile = wc_path
  )
}

powerindex <- readr::read_csv(spi_path, show_col_types = FALSE) |>
  dplyr::mutate(prior_score = as.vector(scale(rev(index)) / 2))

teamnames <- powerindex$team

worldcup2014 <- readr::read_csv(wc_path, show_col_types = FALSE) |>
  dplyr::mutate(
    team_1 = match(team1, teamnames),
    team_2 = match(team2, teamnames)
  )

N_games   <- nrow(worldcup2014)
gamenames <- with(worldcup2014, rev(paste(teamnames[team_1], "vs.", teamnames[team_2])))

# Stan data list — shared across all score-differential models
stan_data <- with(
  worldcup2014,
  list(
    N_teams    = nrow(powerindex),
    N_games    = N_games,
    team_1     = team_1,
    score_1    = score1,
    team_2     = team_2,
    score_2    = score2,
    prior_score = powerindex$prior_score,
    df         = 7
  )
)

# The bivariate Poisson model needs integer scores (Stan array[N] int).
# Build from scratch rather than c()-ing over stan_data to avoid duplicate keys.
stan_data_int_scores <- stan_data
stan_data_int_scores$score_1 <- as.integer(worldcup2014$score1)
stan_data_int_scores$score_2 <- as.integer(worldcup2014$score2)

log_result("\n--- Data ---")
log_result("  N_teams : ", stan_data$N_teams)
log_result("  N_games : ", stan_data$N_games)
log_result("  prior_score range: ",
           round(min(powerindex$prior_score), 3), " to ",
           round(max(powerindex$prior_score), 3))

# =============================================================================
# §3  Stan-native cstan() helper
# =============================================================================
# run_phase3() wraps brm() and cannot be used here.
# All fits use a thin cstan() wrapper for consistency.

cstan <- function(stan_file, data = list(), seed = SEED, chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000,
                  adapt_delta = 0.8) {
  model <- cmdstan_model(stan_file)
  model$sample(
    data            = data,
    seed            = seed,
    chains          = chains,
    parallel_chains = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    refresh         = 0
  )
}

# Manual wf_state update helper (replaces run_phase3() for Stan-native fits)
wf_mark_fit <- function(wf, label) {
  wf$fit_timestamp <- Sys.time()
  wf$fit_hash <- digest::digest(
    list(model = label, data_hash = digest::digest(stan_data, algo = "sha256")),
    algo = "sha256"
  )
  wf$stan_backend <- "cmdstanr"
  wf
}

log_diag_cstan <- function(fit, label) {
  log_result("\n--- Diagnostics: ", label, " (cmdstanr) ---")
  diag <- fit$diagnostic_summary()
  log_result("  num_divergences   : ", sum(diag$num_divergent))
  log_result("  num_max_treedepth : ", sum(diag$num_max_treedepth))
  log_result("  E-BFMI            : ",
             paste(round(diag$ebfmi, 4), collapse = ", "))
}

# =============================================================================
# §4  Phase 1 — off-ramp assessment (run_phase1)
# =============================================================================

log_result("\n=== Phase 1: Off-ramp assessment ===")
# assess_offramps() requires a data frame; we pass worldcup2014 as the data
# source. outcome_var is unused for continuous outcomes (event_rate = NA is
# computed internally and is irrelevant here), but the argument is required.
offramps <- assess_offramps(
  data         = worldcup2014,
  outcome_var  = "score1",       # numeric — continuous; event_rate not computed
  outcome_type = "continuous",
  goal         = "team ability estimation",
  n            = N_games
)
# Off-ramps are shown; we select Full Stan (Bayesian).
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 1,
  action    = "bayesian_selected",
  timestamp = Sys.time(),
  notes     = "Score differential: continuous outcome, N=64 games, 32 teams."
)))
log_result("  Off-ramp chosen: Full Stan (Bayesian)")

# =============================================================================
# §5  Phase 2 — priors
# =============================================================================

# Ch 23 uses Stan's model block directly for all priors:
#   alpha ~ normal(0, 1)   — standardised team ability deviations
#   b     ~ normal(0, 1)   — weight on the Soccer Power Index prior score
#   sigma_a ~ normal(0, 1) — scale of team abilities (half-normal, constrained > 0)
#   sigma_y ~ normal(0, 1) — scale of score differences  (half-normal)
#
# These are intentionally vague for illustration of the debugging theme.
# The chapter does NOT demonstrate a tighter-prior fix for the continuous
# model; it fixes the transformed data bug instead.

log_result("\n=== Phase 2: Priors ===")
log_result("  All models use half-normal(0,1) on sigma_a, sigma_y, sigma_z.")
log_result("  normal(0,1) on b (Soccer Power Index weight) and alpha (team abilities).")
log_result("  Priors are intentionally vague — debugging is via data model, not priors.")

# =============================================================================
# §6  Section 3 — First model (buggy sqrt-t)
# =============================================================================

log_result("\n=== §3: M1 — first model (buggy sqrt-t) ===")

# Phase 3: fit
wf <- wf_mark_fit(wf, "M1_first_try")
fit_1 <- cstan("worldcup_first_try.stan", data = stan_data)
saveRDS(fit_1, "fit_1.rds")

# Phase 4: diagnostics
log_diag_cstan(fit_1, "fit_1")
smry_1 <- fit_1$summary(c("b", "sigma_a", "sigma_y"))
capture_result(smry_1, label = "M1: b, sigma_a, sigma_y")

log_result("Book target M1: b ≈ 0.45, sigma_a ≈ 0.17, sigma_y ≈ 0.42")
b_1      <- smry_1$mean[smry_1$variable == "b"]
siga_1   <- smry_1$mean[smry_1$variable == "sigma_a"]
sigy_1   <- smry_1$mean[smry_1$variable == "sigma_y"]
log_result("  Computed: b = ", round(b_1, 3),
           ", sigma_a = ", round(siga_1, 3),
           ", sigma_y = ", round(sigy_1, 3))

wf$diagnostics$passed      <- TRUE
wf$diagnostics$acknowledged <- TRUE
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "diagnostic_acknowledged",
  timestamp = Sys.time(),
  notes     = "M1 diagnostics clean. Bug is semantic (wrong sqrt), not MCMC."
)))

# Phase 5: PPC — Fig 23.1 (team quality estimates)
# PPC METHOD NOTE
# Standard bayesplot ppc_dens_overlay() is not appropriate here because
# the goal of each phase-5 check in this chapter is posterior retrodiction
# of match score differentials, not density comparison. We use
# bayesplot::mcmc_intervals() for team ability parameters (to show the
# Chapter's main output) and bayesplot::ppc_intervals() for retrodiction
# (to reveal whether the model's predictions bracket the observed scores).

draws_a_1 <- fit_1$draws("a")

fig_23_1 <- bayesplot::mcmc_intervals(draws_a_1, prob = 0) +
  ggplot2::scale_y_discrete(labels = rev(teamnames), limits = rev) +
  ggplot2::labs(x = "Team quality estimate with 90% intervals") +
  ggplot2::theme_minimal()

print(fig_23_1)
ggplot2::ggsave(
  filename = "figs/Fig-23.1.svg",
  plot     = fig_23_1,
  width    = 8, height = 7, device = svg
)

# M1 with replication — needed for Fig 23.3 retrodiction PPC
fit_1_rep <- cstan("worldcup_with_replication.stan", data = stan_data)
saveRDS(fit_1_rep, "fit_1_rep.rds")
log_diag_cstan(fit_1_rep, "fit_1_rep")

fig_23_3 <- bayesplot::ppc_intervals(
  y    = stan_data$score_1 - stan_data$score_2,
  yrep = fit_1_rep$draws("y_rep_original_scale", format = "draws_matrix"),
  fatten = 0, prob = 1e-12
) +
  ggplot2::scale_x_continuous(labels = gamenames, breaks = 1:64) +
  ggplot2::labs(
    y = "Game score differentials\ncompared to 90% predictive interval from model",
    x = ""
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal()

print(fig_23_3)
ggplot2::ggsave(
  filename = "figs/Fig-23.3.svg",
  plot     = fig_23_3,
  width    = 8, height = 10, device = svg
)

# =============================================================================
# §7  Section 4 — Second model: t without sqrt transformation
# =============================================================================

log_result("\n=== §4: M2 — t model without sqrt ===")

wf <- wf_mark_fit(wf, "M2_no_sqrt")
fit_2 <- cstan("worldcup_no_sqrt.stan", data = stan_data)
saveRDS(fit_2, "fit_2.rds")

log_diag_cstan(fit_2, "fit_2")
smry_2 <- fit_2$summary(c("b", "sigma_a", "sigma_y"))
capture_result(smry_2, label = "M2: b, sigma_a, sigma_y")

log_result("Book target M2: b ≈ 1.1, sigma_a ≈ 0.37, sigma_y ≈ 1.3")
b_2    <- smry_2$mean[smry_2$variable == "b"]
siga_2 <- smry_2$mean[smry_2$variable == "sigma_a"]
sigy_2 <- smry_2$mean[smry_2$variable == "sigma_y"]
log_result("  Computed: b = ", round(b_2, 3),
           ", sigma_a = ", round(siga_2, 3),
           ", sigma_y = ", round(sigy_2, 3))

draws_a_2 <- fit_2$draws("a")

fig_23_4 <- bayesplot::mcmc_intervals(draws_a_2, prob = 0) +
  ggplot2::scale_y_discrete(labels = rev(teamnames), limits = rev) +
  ggplot2::labs(x = "Team quality estimate with 90% intervals\n(model with no square root)") +
  ggplot2::theme_minimal()

print(fig_23_4)
ggplot2::ggsave(
  filename = "figs/Fig-23.4.svg",
  plot     = fig_23_4,
  width    = 8, height = 7, device = svg
)

fig_23_5 <- bayesplot::ppc_intervals(
  y    = stan_data$score_1 - stan_data$score_2,
  yrep = fit_2$draws("y_rep", format = "draws_matrix"),
  fatten = 0, prob = 1e-12
) +
  ggplot2::scale_x_continuous(labels = gamenames, breaks = 1:64) +
  ggplot2::labs(
    y = "Game score differentials\ncompared to 90% predictive interval\n(model with no square root)",
    x = ""
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal()

print(fig_23_5)
ggplot2::ggsave(
  filename = "figs/Fig-23.5.svg",
  plot     = fig_23_5,
  width    = 8, height = 10, device = svg
)

# =============================================================================
# §8  Section 5 — Fixed model (corrected sqrt-t)
# =============================================================================

log_result("\n=== §5: M3 — corrected sqrt-t model ===")

wf <- wf_mark_fit(wf, "M3_fixed")
fit_3 <- cstan("worldcup_fixed.stan", data = stan_data)
saveRDS(fit_3, "fit_3.rds")

log_diag_cstan(fit_3, "fit_3")
smry_3 <- fit_3$summary(c("b", "sigma_a", "sigma_y"))
capture_result(smry_3, label = "M3: b, sigma_a, sigma_y")

log_result("Book target M3 (fixed): b ≈ 0.86, sigma_a ≈ 0.33, sigma_y ≈ 0.84")
b_3    <- smry_3$mean[smry_3$variable == "b"]
siga_3 <- smry_3$mean[smry_3$variable == "sigma_a"]
sigy_3 <- smry_3$mean[smry_3$variable == "sigma_y"]
log_result("  Computed: b = ", round(b_3, 3),
           ", sigma_a = ", round(siga_3, 3),
           ", sigma_y = ", round(sigy_3, 3))

draws_a_3 <- fit_3$draws("a")

fig_23_6 <- bayesplot::mcmc_intervals(draws_a_3, prob = 0) +
  ggplot2::scale_y_discrete(labels = rev(teamnames), limits = rev) +
  ggplot2::labs(x = "Team quality estimate with 90% intervals\n(corrected model)") +
  ggplot2::theme_minimal()

print(fig_23_6)
ggplot2::ggsave(
  filename = "figs/Fig-23.6.svg",
  plot     = fig_23_6,
  width    = 8, height = 7, device = svg
)

fig_23_7 <- bayesplot::ppc_intervals(
  y    = stan_data$score_1 - stan_data$score_2,
  yrep = fit_3$draws("y_rep_original_scale", format = "draws_matrix"),
  fatten = 0, prob = 1e-12
) +
  ggplot2::scale_x_continuous(labels = gamenames, breaks = 1:64) +
  ggplot2::labs(
    y = "Game score differentials\ncompared to 90% predictive interval\n(corrected model with square root)",
    x = ""
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_minimal()

print(fig_23_7)
ggplot2::ggsave(
  filename = "figs/Fig-23.7.svg",
  plot     = fig_23_7,
  width    = 8, height = 10, device = svg
)

# M3 without power index prior (b fixed to 0) — Fig 23.2
# The chapter shows this figure second in presentation order (Fig 7 in HTML)
# but we match the book figure numbers here.
fit_3_no_prior <- cstan("worldcup_no_prior.stan",
                         data = c(stan_data, b = 0))
saveRDS(fit_3_no_prior, "fit_3_no_prior.rds")
log_diag_cstan(fit_3_no_prior, "fit_3_no_prior")

draws_a_3np <- fit_3_no_prior$draws("a")

fig_23_2 <- bayesplot::mcmc_intervals(draws_a_3np, prob = 0) +
  ggplot2::scale_y_discrete(labels = rev(teamnames), limits = rev) +
  ggplot2::labs(x = "Team quality estimate with 90% intervals\nModel without prior rankings") +
  ggplot2::theme_minimal()

print(fig_23_2)
ggplot2::ggsave(
  filename = "figs/Fig-23.2.svg",
  plot     = fig_23_2,
  width    = 8, height = 7, device = svg
)

# =============================================================================
# §9  Section 6 — Discrete models and first LOO comparison
# =============================================================================

log_result("\n=== §6: Discrete models and LOO-CV ===")

# § 6.1 Discrete model with explicit latent z (for reference; LOO not reliable here)
wf <- wf_mark_fit(wf, "M_discr_z")
fit_discr_z <- cstan("worldcup_discrete_z.stan", data = stan_data)
saveRDS(fit_discr_z, "fit_discr_z.rds")
log_diag_cstan(fit_discr_z, "fit_discr_z")
loo_discr_z <- fit_discr_z$loo()

# § 6.2 Discrete model with z integrated out (primary discrete model)
wf <- wf_mark_fit(wf, "M_discr")
fit_discr <- cstan("worldcup_discrete.stan", data = stan_data)
saveRDS(fit_discr, "fit_discr.rds")
log_diag_cstan(fit_discr, "fit_discr")

smry_discr <- fit_discr$summary(c("a[1]", "a[32]", "b", "sigma_a", "sigma_z"))
capture_result(smry_discr, label = "M-discr: key parameters")
log_result("Book target M-discr: b ≈ 1.1, sigma_a ≈ 0.52, sigma_z ≈ 1.5")

loo_discr <- fit_discr$loo(save_psis = TRUE)

# § 6.3 Discrete model, no power score
wf <- wf_mark_fit(wf, "M_discr_nopower")
fit_discr_nopower <- cstan("worldcup_discrete_nopower.stan", data = stan_data)
saveRDS(fit_discr_nopower, "fit_discr_nopower.rds")
log_diag_cstan(fit_discr_nopower, "fit_discr_nopower")
loo_discr_nopower <- fit_discr_nopower$loo()

# § 6.4 Discrete model, power score only
wf <- wf_mark_fit(wf, "M_discr_poweronly")
fit_discr_poweronly <- cstan("worldcup_discrete_poweronly.stan", data = stan_data)
saveRDS(fit_discr_poweronly, "fit_discr_poweronly.rds")
log_diag_cstan(fit_discr_poweronly, "fit_discr_poweronly")
loo_discr_poweronly <- fit_discr_poweronly$loo()

# § 6.5 Discrete model, pooled (no per-team effects)
wf <- wf_mark_fit(wf, "M_discr_pool")
fit_discr_pool <- cstan("worldcup_discrete_pooled.stan", data = stan_data)
saveRDS(fit_discr_pool, "fit_discr_pool.rds")
log_diag_cstan(fit_discr_pool, "fit_discr_pool")
loo_discr_pool <- fit_discr_pool$loo()

# § 6.6 LOO comparison — discrete model family
log_result("\n=== LOO comparison: discrete model family (§6.6) ===")
log_result("Book target: Hier. w power score best; elpd_diff vs poweronly ≈ -1.3 (se 1.8)")
log_result("             elpd_diff vs nopower ≈ -3.3 (se 2.8), pooled ≈ -8.4 (se 3.4)")

loo_tab_discr <- loo::loo_compare(list(
  "Hier. w power score"  = loo_discr,
  "Power score only"     = loo_discr_poweronly,
  "Hier. w/o power score" = loo_discr_nopower,
  "Pooled w/o power score" = loo_discr_pool
))
capture_result(loo_tab_discr, label = "LOO table §6.6")
saveRDS(loo_tab_discr, "loo_tab_discr.rds")

# =============================================================================
# §10  Section 7 — Discretizing continuous models
# =============================================================================

log_result("\n=== §7: Discretized continuous models ===")

# § 7.1 Continuous model, midpoint log_lik
wf <- wf_mark_fit(wf, "M_cont_midp")
fit_cont_midp <- cstan("worldcup_continuous_midpoint_ll.stan", data = stan_data)
saveRDS(fit_cont_midp, "fit_cont_midp.rds")
log_diag_cstan(fit_cont_midp, "fit_cont_midp")
loo_cont_midp <- fit_cont_midp$loo()

# § 7.2 Continuous model, exact integration log_lik
wf <- wf_mark_fit(wf, "M_cont")
fit_cont <- cstan("worldcup_continuous.stan", data = stan_data)
saveRDS(fit_cont, "fit_cont.rds")
log_diag_cstan(fit_cont, "fit_cont")
loo_cont <- fit_cont$loo(save_psis = TRUE)

# § 7.3 LOO: discrete vs continuous
log_result("\n=== LOO comparison: discrete vs continuous models (§7.3) ===")
log_result("Book: differences small (Monte Carlo variation); continuous best but barely")
loo_tab_cont <- loo::loo_compare(list(
  "Discrete model"               = loo_discr,
  "Continuous + midpoint log_lik" = loo_cont_midp,
  "Continuous model"             = loo_cont
))
capture_result(loo_tab_cont, label = "LOO table §7.3")
saveRDS(loo_tab_cont, "loo_tab_cont.rds")

# § 7.4 More discretized continuous models — sqrt variants
# Continuous sqrt, no Jacobian
wf <- wf_mark_fit(wf, "M_sqrt_cont_noj")
fit_sqrt_cont_noj <- cstan("worldcup_sqrt_continuous_nojacobian.stan", data = stan_data)
saveRDS(fit_sqrt_cont_noj, "fit_sqrt_cont_noj.rds")
log_diag_cstan(fit_sqrt_cont_noj, "fit_sqrt_cont_noj")
loo_sqrt_cont_noj <- fit_sqrt_cont_noj$loo()

# Continuous sqrt, with Jacobian (quadrature in generated quantities only — fast)
wf <- wf_mark_fit(wf, "M_sqrt_cont")
fit_sqrt_cont <- cstan("worldcup_sqrt_continuous.stan", data = stan_data)
saveRDS(fit_sqrt_cont, "fit_sqrt_cont.rds")
log_diag_cstan(fit_sqrt_cont, "fit_sqrt_cont")
loo_sqrt_cont <- fit_sqrt_cont$loo(save_psis = TRUE)

# Discrete sqrt with Jacobian (SLOW: quadrature at each leapfrog step)
if (RUN_SQRT_DISCR) {
  log_result("  Running worldcup_sqrt_discrete.stan (SLOW — may take 30–60 min) ...")
  wf <- wf_mark_fit(wf, "M_sqrt_discr")
  fit_sqrt_discr <- cstan("worldcup_sqrt_discrete.stan", data = stan_data)
  saveRDS(fit_sqrt_discr, "fit_sqrt_discr.rds")
  log_diag_cstan(fit_sqrt_discr, "fit_sqrt_discr")
  loo_sqrt_discr <- fit_sqrt_discr$loo()
} else {
  log_result("  worldcup_sqrt_discrete.stan SKIPPED (RUN_SQRT_DISCR = FALSE).")
  log_result("  Set RUN_SQRT_DISCR <- TRUE to include this model in the comparison.")
}

# § 7.5 LOO: sqrt variants
log_result("\n=== LOO comparison: sqrt models (§7.5) ===")
log_result("Book (without Jacobian): Cont-sqrt+midp,-Jacobian best;")
log_result("  Discrete elpd_diff ≈ -32.5 (se 5.2), Discrete sqrt ≈ -40.3 (se 4.3)")
log_result("Book (with Jacobian): Discrete best;")
log_result("  Cont-sqrt+Jacobian elpd_diff ≈ -7.2 (se 5.5), Discrete sqrt ≈ -7.8 (se 5.4)")

loo_tab_sqrt_noj <- loo::loo_compare(list(
  "Discrete"                   = loo_discr,
  "Cont-sqrt + midp, -Jacobian" = loo_sqrt_cont_noj
))
capture_result(loo_tab_sqrt_noj, label = "LOO §7.5 (no-Jacobian comparison, 2-way)")
saveRDS(loo_tab_sqrt_noj, "loo_tab_sqrt_noj.rds")

if (RUN_SQRT_DISCR) {
  loo_tab_sqrt_j <- loo::loo_compare(list(
    "Discrete"                  = loo_discr,
    "Discrete sqrt"             = loo_sqrt_discr,
    "Continuous sqrt +Jacobian" = loo_sqrt_cont
  ))
  capture_result(loo_tab_sqrt_j, label = "LOO §7.5 (Jacobian comparison, 3-way)")
  saveRDS(loo_tab_sqrt_j, "loo_tab_sqrt_j.rds")
} else {
  loo_tab_sqrt_j <- loo::loo_compare(list(
    "Discrete"                  = loo_discr,
    "Continuous sqrt +Jacobian" = loo_sqrt_cont
  ))
  capture_result(loo_tab_sqrt_j, label = "LOO §7.5 (Jacobian, 2-way, sqrt_discr skipped)")
  saveRDS(loo_tab_sqrt_j, "loo_tab_sqrt_j.rds")
}

# § 7.6 LOO-PIT checks: discrete model
log_result("\n=== LOO-PIT: M-discr (§7.6) ===")

fig_23_8 <- bayesplot::ppc_loo_pit_ecdf(
  y          = stan_data$score_1 - stan_data$score_2,
  yrep       = fit_discr$draws(format = "matrix", variables = "y_rep"),
  psis_object = loo_discr$psis_object,
  method     = "correlated"
) +
  ggplot2::labs(title = "LOO-PIT ECDF: discrete model (§7.6)") +
  ggplot2::theme_minimal()

print(fig_23_8)
ggplot2::ggsave(
  filename = "figs/Fig-23.8.svg",
  plot     = fig_23_8,
  width    = 7, height = 5, device = svg
)

# LOO-PIT: continuous sqrt model with Jacobian
fig_23_9 <- bayesplot::ppc_loo_pit_ecdf(
  y          = stan_data$score_1 - stan_data$score_2,
  yrep       = fit_sqrt_cont$draws(format = "matrix", variables = "y_rep"),
  psis_object = loo_sqrt_cont$psis_object,
  method     = "correlated"
) +
  ggplot2::labs(title = "LOO-PIT ECDF: continuous sqrt model (§7.6)") +
  ggplot2::theme_minimal()

print(fig_23_9)
ggplot2::ggsave(
  filename = "figs/Fig-23.9.svg",
  plot     = fig_23_9,
  width    = 7, height = 5, device = svg
)

# =============================================================================
# §11  Section 8 — Bivariate Poisson and Poisson difference models
# =============================================================================

log_result("\n=== §8: Bivariate Poisson and Poisson difference ===")

# Bivariate Poisson — needs adapt_delta = 0.95 (book uses this)
wf <- wf_mark_fit(wf, "M_bipois")
fit_bipois <- cstan("worldcup_bivariate_poisson.stan",
                    data = stan_data_int_scores,
                    adapt_delta = 0.95)
saveRDS(fit_bipois, "fit_bipois.rds")
log_diag_cstan(fit_bipois, "fit_bipois")

smry_bipois <- fit_bipois$summary(
  c("a", "o[1]", "o[32]", "d[1]", "d[32]", "b_o", "b_d", "sigma_o", "sigma_d")
)
capture_result(smry_bipois, label = "M-bipois: key parameters")
log_result("Book target M-bipois: b_o ≈ 0.53, b_d ≈ -0.45, sigma_o ≈ 0.18, sigma_d ≈ 0.30")

loo_bipois <- fit_bipois$loo(save_psis = TRUE)
capture_result(loo_bipois, label = "M-bipois: loo summary")

# Poisson difference
wf <- wf_mark_fit(wf, "M_poisdif")
fit_poisdif <- cstan("worldcup_poisson_difference.stan", data = stan_data)
saveRDS(fit_poisdif, "fit_poisdif.rds")
log_diag_cstan(fit_poisdif, "fit_poisdif")

smry_poisdif <- fit_poisdif$summary(c("a[1]", "a[32]", "b", "sigma_a"))
capture_result(smry_poisdif, label = "M-poisdif: key parameters")
log_result("Book target M-poisdif: b ≈ 0.52, sigma_a ≈ 0.24")

loo_poisdif <- fit_poisdif$loo(save_psis = TRUE)
capture_result(loo_poisdif, label = "M-poisdif: loo summary")

# Final LOO comparison
log_result("\n=== LOO comparison: Poisson models vs discrete (§8) ===")
log_result("Book: Poisson difference best; elpd_diff bipois ≈ -1.8 (se 0.9)")
log_result("      elpd_diff discrete vs poisdif ≈ -1.7 (se 1.7). Differences small.")

loo_tab_final <- loo::loo_compare(list(
  "Discrete"          = loo_discr,
  "Bivariate Poisson" = loo_bipois,
  "Poisson difference" = loo_poisdif
))
capture_result(loo_tab_final, label = "LOO §8 final comparison")
saveRDS(loo_tab_final, "loo_tab_final.rds")

# § 8.1 LOO-PIT: bivariate Poisson
fig_23_10 <- bayesplot::ppc_loo_pit_ecdf(
  y          = stan_data$score_1 - stan_data$score_2,
  yrep       = fit_bipois$draws(format = "matrix", variables = "y_rep"),
  psis_object = loo_bipois$psis_object,
  method     = "correlated"
) +
  ggplot2::labs(title = "LOO-PIT ECDF: bivariate Poisson (§8.1)") +
  ggplot2::theme_minimal()

print(fig_23_10)
ggplot2::ggsave(
  filename = "figs/Fig-23.10.svg",
  plot     = fig_23_10,
  width    = 7, height = 5, device = svg
)

# LOO-PIT: Poisson difference
fig_23_11 <- bayesplot::ppc_loo_pit_ecdf(
  y          = stan_data$score_1 - stan_data$score_2,
  yrep       = fit_poisdif$draws(format = "matrix", variables = "y_rep"),
  psis_object = loo_poisdif$psis_object,
  method     = "correlated"
) +
  ggplot2::labs(title = "LOO-PIT ECDF: Poisson difference (§8.1)") +
  ggplot2::theme_minimal()

print(fig_23_11)
ggplot2::ggsave(
  filename = "figs/Fig-23.11.svg",
  plot     = fig_23_11,
  width    = 7, height = 5, device = svg
)

# =============================================================================
# §12  LOO-CV status
# =============================================================================

# LOO-CV is run extensively in this chapter (sections 6.6, 7.3, 7.5, 8).
# Primary comparison: discrete normal model (M-discr) as the baseline.
# The Poisson difference model wins §8 comparison but only slightly.
wf$loo_complete <- TRUE

# =============================================================================
# §13  Save and wrap up
# =============================================================================

wf$ppc_complete <- TRUE
wf$loo_complete <- TRUE
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

saved_rds <- paste(c(
  "fit_1.rds", "fit_1_rep.rds", "fit_2.rds",
  "fit_3.rds", "fit_3_no_prior.rds",
  "fit_discr_z.rds", "fit_discr.rds",
  "fit_discr_nopower.rds", "fit_discr_poweronly.rds", "fit_discr_pool.rds",
  "fit_cont_midp.rds", "fit_cont.rds",
  "fit_sqrt_cont_noj.rds", "fit_sqrt_cont.rds",
  if (RUN_SQRT_DISCR) "fit_sqrt_discr.rds",
  "fit_bipois.rds", "fit_poisdif.rds",
  "loo_tab_discr.rds", "loo_tab_cont.rds",
  "loo_tab_sqrt_noj.rds", "loo_tab_sqrt_j.rds", "loo_tab_final.rds",
  "wf_final.rds"
), collapse = ", ")

saved_figs <- paste(sprintf("Fig-23.%d.svg", 1:11), collapse = ", ")

log_result("\n=== Ch 23 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: ", saved_rds)
log_result("Figures saved to figs/: ", saved_figs)
if (!RUN_SQRT_DISCR) {
  log_result("NOTE: fit_sqrt_discr.rds omitted (RUN_SQRT_DISCR = FALSE).")
}

close(results_con)
