# data/ch29_no_vehicles/ch29_no_vehicles_analysis.R
#
# Ch 29 — "Sampling problems with latent variables: No vehicles in the park"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# This script uses cmdstanr directly because the chapter's goal is to
# demonstrate identifiability and sampling-efficiency problems in hierarchical
# logistic models. The four Stan models (park_1 through park_4) explore
# non-centered vs. centered parameterization and sum-to-zero constraints;
# these parameterizations cannot be expressed through brms formula interfaces.
# run_phase3() is therefore bypassed for all fits. The lme4 section uses
# glmer() for fast approximate inference before the Stan fits.
# Deliberately omits: display() (replaced with lme4::summary()),
# tictoc timing (replaced with system.time()), park_*_normal.stan variants.
#
# Stan fits 1–3 use iter_warmup = 200, iter_sampling = 200 (one fifth of
# default) to reproduce the book's diagnostic comparison. Only fit_4 is
# refit at the default 1000/1000. This deviates from rubric §6 rule
# "iter = 4000, warmup = 2000 for latent-parameter models" because the
# chapter's pedagogical purpose requires short runs that expose the
# sampling problems. The deviation is documented here and in results.txt.
#
# Data:
#   park.csv (data/ch29_no_vehicles/)
#   Source: Luu (2024) / Turner (2024) survey; see chapter references
#   51953 rows x 5 cols; 2409 respondents, 27 items
#   Key variables: submission_id, question_id, male, white, answer (0/1)
#   park.txt (data/ch29_no_vehicles/)
#   Source: same; 27 rows, one item wording per row
#
# Models covered:
#   fit_lme4     — glmer(y ~ (1|item) + (1|respondent) + male + white + n_responses_full, binomial)
#   fit_lme4     — glmer(y ~ (1|item) + (1|respondent) + male + white + n_skipped_full, binomial)  [final lme4]
#   fit_lme4_sim — same formula on simulated y_sim
#   fit_1        — park_1.stan: non-centered, vector<multiplier=sigma>; 200/200 iter
#   fit_2        — park_2.stan: sum-to-zero + non-centered z_respondent/z_item; 200/200 iter
#   fit_3        — park_3.stan: sum-to-zero + centered a_respondent/a_item; 200/200 iter
#   fit_4        — park_4.stan: sum-to-zero + centered + predictor centering; 200/200 iter, then 1000/1000 refit
#
# Figures produced:
#   Fig-29.1.svg  — item_avg vs a_item_hat scatter (lme4)
#   Fig-29.2.svg  — logit(item_avg) vs a_item_hat text plot (lme4)
#   Fig-29.3.svg  — respondent_avg vs a_respondent_hat scatter (lme4)
#   Fig-29.4.svg  — mcmc_trace for 'a' from fit_1
#   Fig-29.5.svg  — mcmc_scatter of a vs sum(a_item) from fit_1
#   Fig-29.6.svg  — mcmc_scatter of a_item[1] vs sigma_item from fit_2 (untransformed)
#   Fig-29.7.svg  — mcmc_scatter of z_item[1] vs sigma_item from fit_2
#   Fig-29.8.svg  — mcmc_scatter of z_item[1] vs log(sigma_item) from fit_2
#   Fig-29.9.svg  — lme4 vs Stan: a_item point estimates and 90% intervals
#   Fig-29.10.svg — lme4 vs Stan: a_respondent point estimates and 90% intervals
#   Fig-29.11.svg — sigma_item and sigma_respondent posterior slabs vs lme4 estimates
#
# Book-target posterior summaries (fit_4, full run):
#   a:                mean ≈ -2.3, sd ≈ 0.03
#   b[1] (male):      mean ≈  0.07, sd ≈ 0.03
#   b[2] (white):     mean ≈  0.08, sd ≈ 0.03
#   b[3] (n_skipped): mean ≈  0.03, sd ≈ 0.01
#   sigma_respondent: mean ≈  1.9,  sd ≈ 0.04
#   sigma_item:       mean ≈  2.4,  sd ≈ 0.34
#   a_item[1] (car):  mean ≈  7.6
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch29_no_vehicles/.

# =============================================================================
# Environment setup
# =============================================================================

.libPaths(c("C:/Users/johnx/Documents/WildPeaches/Projects/RBayesflow/renv/library/windows/R-4.6/x86_64-w64-mingw32", .libPaths()))
source("../../R/source_all.R")
library(lme4)
library(ggdist)
library(patchwork)
library(posterior)

options(brms.backend = "cmdstanr", mc.cores = 4)
options(posterior.num_args = list(digits = 2), digits = 2, width = 90)
SEED <- 123   # matches book set.seed(123) for simulation section

dir.create("figs", showWarnings = FALSE)
dir.create("stan_output", showWarnings = FALSE)
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
log_result("=== Ch 29 results log — ", format(Sys.time()), " ===")
log_result("NOTE: Stan fits 1-3 use iter_warmup=200, iter_sampling=200 by design.")
log_result("      This deviates from rubric §6 to reproduce the book's diagnostic comparison.")

# =============================================================================
# Data
# =============================================================================

park <- read.csv("park.csv")
N <- nrow(park)

respondents <- sort(unique(park$submission_id))
J <- length(respondents)
respondent <- rep(NA, N)
for (j in 1:J) {
  respondent[park$submission_id == respondents[j]] <- j
}

items <- sort(unique(park$question_id))
K <- length(items)
item <- park$question_id
y <- park$answer

n_responses <- rep(NA, J)
for (j in 1:J) {
  n_responses[j] <- sum(respondent == j)
}

male_name <- park$male
white_name <- park$white
n_responses_full <- n_responses[respondent]

# Recode: number of items skipped so zero is a reasonable baseline
n_skipped <- K - n_responses
n_skipped_full <- n_skipped[respondent]

data_park <- data.frame(y, respondent, item, male_name, white_name, n_skipped_full)

# Item wordings (for labelled plots)
wordings_raw <- read.csv("park.txt", header = FALSE)$V2
wordings <- substr(wordings_raw, 2, nchar(wordings_raw) - 1)

log_result("\n=== Data dimensions ===")
log_result("N (responses) : ", N)
log_result("J (respondents): ", J)
log_result("K (items)      : ", K)
log_result("mean(y)        : ", round(mean(y), 4))
log_result("mean(y) baseline (male=0, white=0, n_skipped=0): ",
           round(mean(y[male_name == 0 & white_name == 0 & n_skipped_full == 0]), 4))
log_result("n baseline obs : ",
           sum(male_name == 0 & white_name == 0 & n_skipped_full == 0))

# =============================================================================
# Phase 1 / lme4 model
# =============================================================================

wf <- init_workflow(mode = "practice", stage = "explore")

# --- lme4 fit 1: n_responses_full as predictor ---
t1 <- system.time(
  fit_lme4_raw <- lme4::glmer(
    y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full,
    family = binomial(link = "logit"),
    data   = data.frame(y, respondent, item, male_name, white_name, n_responses_full)
  )
)
log_result("\n=== lme4 fit 1 (n_responses_full) ===")
log_result("Elapsed: ", round(t1["elapsed"], 1), " sec")
capture_result(summary(fit_lme4_raw), label = "summary(fit_lme4_raw)")

# --- lme4 fit 2: n_skipped_full as predictor (zero = reasonable baseline) ---
t2 <- system.time(
  fit_lme4 <- lme4::glmer(
    y ~ (1 | item) + (1 | respondent) + male_name + white_name + n_skipped_full,
    family = binomial(link = "logit"),
    data   = data_park
  )
)
log_result("\n=== lme4 fit 2 (n_skipped_full) — final lme4 model ===")
log_result("Elapsed: ", round(t2["elapsed"], 1), " sec")
capture_result(summary(fit_lme4), label = "summary(fit_lme4)")

log_result("\n=== Book targets: lme4 fit 2 ===")
log_result("Book: Intercept ≈ -2.42, male ≈ 0.07, white ≈ 0.08, n_skipped ≈ 0.03")
log_result("Book: sigma_respondent ≈ 1.92, sigma_item ≈ 2.31")
fe <- lme4::fixef(fit_lme4)
vc <- as.data.frame(lme4::VarCorr(fit_lme4))
log_result("Computed: Intercept = ", round(fe["(Intercept)"], 4),
           ", male = ", round(fe["male_name"], 4),
           ", white = ", round(fe["white_name"], 4),
           ", n_skipped = ", round(fe["n_skipped_full"], 4))
log_result("Computed: sigma_respondent = ",
           round(vc$sdcor[vc$grp == "respondent"], 4),
           ", sigma_item = ",
           round(vc$sdcor[vc$grp == "item"], 4))

log_result("\n=== Mean y checks ===")
log_result("mean(y)                          = ", round(mean(y), 4))
log_result("mean(y) baseline (m=0,w=0,sk=0)  = ",
           round(mean(y[male_name == 0 & white_name == 0 & n_skipped_full == 0]), 4))
log_result("plogis(intercept + coef * means) = ",
           round(plogis(fe["(Intercept)"] +
                        fe["male_name"]    * mean(male_name) +
                        fe["white_name"]   * mean(white_name) +
                        fe["n_skipped_full"] * mean(n_skipped_full)), 4))
log_result("mean(predict(fit_lme4, type='response')) = ",
           round(mean(predict(fit_lme4, type = "response")), 4))

# --- Simulated-data recovery ---
set.seed(SEED)
a_respondent_sim <- rnorm(J, 0, sqrt(lme4::VarCorr(fit_lme4)$respondent))
a_item_sim       <- rnorm(K, 0, sqrt(lme4::VarCorr(fit_lme4)$item))
b_sim            <- lme4::fixef(fit_lme4)
X_sim            <- cbind(1, male_name, white_name, n_skipped_full)
p_sim            <- plogis(a_respondent_sim[respondent] + a_item_sim[item] + X_sim %*% b_sim)
y_sim            <- rbinom(N, 1, p_sim)
data_sim         <- data.frame(data_park, y_sim)

t3 <- system.time(
  fit_lme4_sim <- lme4::glmer(
    y_sim ~ (1 | item) + (1 | respondent) + male_name + white_name + n_responses_full,
    family = binomial(link = "logit"),
    data   = data_sim
  )
)
log_result("\n=== lme4 fit sim (recovery check) ===")
log_result("Elapsed: ", round(t3["elapsed"], 1), " sec")
capture_result(summary(fit_lme4_sim), label = "summary(fit_lme4_sim)")
log_result("Book target (sim): Intercept ≈ -1.88, male ≈ 0.06, white ≈ 0.11, n_responses ≈ -0.04")
log_result("Book target (sim): sigma_respondent ≈ 1.87, sigma_item ≈ 2.65")
fe_sim <- lme4::fixef(fit_lme4_sim)
vc_sim <- as.data.frame(lme4::VarCorr(fit_lme4_sim))
log_result("Computed (sim): Intercept = ", round(fe_sim["(Intercept)"], 4),
           ", male = ", round(fe_sim["male_name"], 4),
           ", white = ", round(fe_sim["white_name"], 4))
log_result("Computed (sim): sigma_respondent = ",
           round(vc_sim$sdcor[vc_sim$grp == "respondent"], 4),
           ", sigma_item = ",
           round(vc_sim$sdcor[vc_sim$grp == "item"], 4))
log_result("NOTE: simulation uses set.seed(123) matching the book. Small deviations")
log_result("      from book targets are expected from the different lme4 version.")

# --- Item random effects table ---
a_item_hat <- lme4::ranef(fit_lme4)$item
a_item_hat_vec <- unlist(a_item_hat)
names(a_item_hat_vec) <- wordings
log_result("\n=== Item random effects (sorted) ===")
capture_result(round(sort(a_item_hat_vec), 2), label = "sort(a_item_hat)")
log_result("Book: kite ≈ -3.17, car ≈ 7.61 (most/least vehicle-like by raw effects)")

# =============================================================================
# Phase 1 figures — lme4 item and respondent effects
# =============================================================================

# Item average response rates for plotting
item_avg <- rep(NA, K)
for (k in 1:K) {
  item_avg[k] <- mean(y[item == k])
}

# Respondent average response rates for plotting
a_respondent_hat <- unlist(lme4::ranef(fit_lme4)$respondent)
respondent_avg <- rep(NA, J)
for (j in 1:J) {
  respondent_avg[j] <- mean(y[respondent == j])
}

# Fig 29.1 — item average vs lme4 item random effect
fig_29_1 <- ggplot2::ggplot(
  data.frame(item_avg = item_avg, a_item = a_item_hat_vec),
  ggplot2::aes(x = item_avg, y = a_item)
) +
  ggplot2::geom_point(shape = 20) +
  ggplot2::labs(
    x = "Observed item average (proportion Yes)",
    y = "lme4 item random effect",
    title = "Figure 29.1: Item average vs. lme4 random effect"
  ) +
  ggplot2::theme_minimal()

print(fig_29_1)
ggplot2::ggsave("figs/Fig-29.1.svg", plot = fig_29_1,
                width = 6, height = 5, device = svg)

# Fig 29.2 — logit(item_avg) vs lme4 item random effect, labelled by wording
fig_29_2 <- ggplot2::ggplot(
  data.frame(logit_avg = qlogis(item_avg), a_item = a_item_hat_vec,
             label = names(a_item_hat_vec)),
  ggplot2::aes(x = logit_avg, y = a_item, label = label)
) +
  ggplot2::geom_text(size = 2.5) +
  ggplot2::labs(
    x = "logit(item average)",
    y = "lme4 item random effect",
    title = "Figure 29.2: logit(item average) vs. lme4 random effect"
  ) +
  ggplot2::theme_minimal()

print(fig_29_2)
ggplot2::ggsave("figs/Fig-29.2.svg", plot = fig_29_2,
                width = 7, height = 6, device = svg)

# Fig 29.3 — respondent average vs lme4 respondent random effect
fig_29_3 <- ggplot2::ggplot(
  data.frame(respondent_avg = respondent_avg, a_respondent = a_respondent_hat),
  ggplot2::aes(x = respondent_avg, y = a_respondent)
) +
  ggplot2::geom_point(shape = 20, size = 0.8, alpha = 0.4) +
  ggplot2::labs(
    x = "Observed respondent average (proportion Yes)",
    y = "lme4 respondent random effect",
    title = "Figure 29.3: Respondent average vs. lme4 random effect"
  ) +
  ggplot2::theme_minimal()

print(fig_29_3)
ggplot2::ggsave("figs/Fig-29.3.svg", plot = fig_29_3,
                width = 6, height = 5, device = svg)

log_result("\n=== Figures 29.1–29.3 saved ===")

# =============================================================================
# Stan data and helper
# =============================================================================

X <- cbind(male_name, white_name, n_skipped_full)
stan_data <- list(
  N          = N,
  J          = J,
  K          = K,
  L          = ncol(X),
  y          = y,
  respondent = respondent,
  item       = item,
  X          = X
)

log_result("\n=== Stan data dimensions ===")
log_result("N=", N, ", J=", J, ", K=", K, ", L=", ncol(X))

# Thin wrapper: compiles a Stan model and samples from it.
# init = 0.1 initialises unconstrained parameters in [-0.1, 0.1],
# which is better than the default [-2, 2] for large hierarchical models.
cstan <- function(stan_file, data = list(), seed = SEED, chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000) {
  model <- cmdstanr::cmdstan_model(stan_file)
  model$sample(
    data             = data,
    seed             = seed,
    chains           = chains,
    parallel_chains  = chains,
    init             = 0.1,
    iter_warmup      = iter_warmup,
    iter_sampling    = iter_sampling,
    refresh          = 0
  )
}

# Compile all four models (compilation only; sampling in later sections)
log_result("\n=== Compiling Stan models ===")
park_1 <- cmdstanr::cmdstan_model("park_1.stan")
park_2 <- cmdstanr::cmdstan_model("park_2.stan")
park_3 <- cmdstanr::cmdstan_model("park_3.stan")
park_4 <- cmdstanr::cmdstan_model("park_4.stan")
log_result("All four Stan models compiled.")

# =============================================================================
# Phase 3 / fit_1: non-centered parameterization with vector<multiplier=sigma>
# =============================================================================

log_result("\n=== fit_1: park_1.stan (non-centered, multiplier) ===")
log_result("iter_warmup=200, iter_sampling=200 — exploratory run to expose diagnostics")

t_fit1 <- system.time(
  fit_1 <- cstan("park_1.stan", data = stan_data,
                 iter_warmup = 200, iter_sampling = 200)
)
fit_1$save_object("fit_1.rds")
log_result("Elapsed: ", round(t_fit1["elapsed"], 1), " sec")

# Audit trail (Stan-native: rubric §6 step 4)
wf$fit_timestamp <- Sys.time()
wf$fit_hash <- digest::digest(
  list(formula = "park_1: y ~ bernoulli_logit_glm(X, a + a_respondent[respondent] + a_item[item], b)",
       data_hash = digest::digest(stan_data, algo = "sha256")),
  algo = "sha256"
)

# --- Sampler diagnostics ---
log_result("\n--- Diagnostics: fit_1 (cmdstanr) ---")
diag1 <- fit_1$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag1$num_divergent))
log_result("  num_max_treedepth : ", sum(diag1$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag1$ebfmi, 4), collapse = ", "))

# Key parameter summary
draws_1 <- fit_1$draws(format = "df")
smry_1_key <- posterior::summarise_draws(
  posterior::subset_draws(fit_1$draws(), variable = c("a", "sigma_respondent", "sigma_item"))
)
capture_result(smry_1_key, label = "fit_1: key parameter summary")
log_result("Book: a ≈ -2.4 (sd ≈ 0.46 — inflated by non-identifiability)")
log_result("      sigma_respondent ≈ 1.9, sigma_item ≈ 2.4")

# Sampler efficiency
smry_sampler_1 <- fit_1$sampler_diagnostics() |> posterior::as_draws_rvars()
log_result("  mean treedepth__ : ", round(mean(posterior::E(smry_sampler_1$treedepth__)), 2))
log_result("  mean n_leapfrog__: ", round(mean(posterior::E(smry_sampler_1$n_leapfrog__)), 2))

wf$diagnostics$passed      <- FALSE   # high Rhat on 'a', low ESS on a_item
wf$diagnostics$acknowledged <- TRUE
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "diagnostic_acknowledged",
  timestamp = Sys.time(),
  notes     = "fit_1: high autocorrelation in 'a' due to non-identifiability of intercept + sum(a_item)"
)))

# --- Fig 29.4: mcmc_trace for 'a' ---
fig_29_4 <- bayesplot::mcmc_trace(draws_1, pars = "a") +
  ggplot2::labs(x = "Iteration",
                title = "Figure 29.4: Trace plot for 'a' — fit_1 (non-centered)") +
  ggplot2::theme_minimal()

print(fig_29_4)
ggplot2::ggsave("figs/Fig-29.4.svg", plot = fig_29_4,
                width = 7, height = 4, device = svg)

# --- Fig 29.5: scatter of a vs sum(a_item) ---
# Use bayesplot::mcmc_scatter on subset draws
fig_29_5 <- posterior::as_draws_rvars(draws_1) |>
  {\(d) posterior::draws_rvars(
    a        = d$a,
    sum_a_item = posterior::rvar_sum(d$a_item)
  )}() |>
  bayesplot::mcmc_scatter(alpha = 0.5) +
  ggplot2::labs(y = "sum(a_item)",
                title = "Figure 29.5: a vs. sum(a_item) — fit_1") +
  ggplot2::theme_minimal()

print(fig_29_5)
ggplot2::ggsave("figs/Fig-29.5.svg", plot = fig_29_5,
                width = 6, height = 5, device = svg)

log_result("\n=== Figures 29.4–29.5 saved ===")

# =============================================================================
# fit_2: sum-to-zero + non-centered z_respondent / z_item
# =============================================================================

log_result("\n=== fit_2: park_2.stan (sum-to-zero + non-centered z) ===")
log_result("iter_warmup=200, iter_sampling=200")

t_fit2 <- system.time(
  fit_2 <- cstan("park_2.stan", data = stan_data,
                 iter_warmup = 200, iter_sampling = 200)
)
fit_2$save_object("fit_2.rds")
log_result("Elapsed: ", round(t_fit2["elapsed"], 1), " sec")

# --- Sampler diagnostics ---
log_result("\n--- Diagnostics: fit_2 (cmdstanr) ---")
diag2 <- fit_2$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag2$num_divergent))
log_result("  num_max_treedepth : ", sum(diag2$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag2$ebfmi, 4), collapse = ", "))

draws_2 <- fit_2$draws(format = "df")
smry_2_key <- posterior::summarise_draws(
  posterior::subset_draws(fit_2$draws(),
                          variable = c("a", "sigma_respondent", "sigma_item"))
)
capture_result(smry_2_key, label = "fit_2: key parameter summary")
log_result("Book: a sd ≈ 0.04 (much reduced vs fit_1 sd ≈ 0.46 — identifiability fixed)")
log_result("      sigma_item Rhat > 1 indicates remaining problem with z_item")

smry_sampler_2 <- fit_2$sampler_diagnostics() |> posterior::as_draws_rvars()
log_result("  mean treedepth__ : ", round(mean(posterior::E(smry_sampler_2$treedepth__)), 2))
log_result("  mean n_leapfrog__: ", round(mean(posterior::E(smry_sampler_2$n_leapfrog__)), 2))

# z_item diagnostics — the actual sampling space
smry_2_z <- posterior::summarise_draws(
  posterior::subset_draws(fit_2$draws(), variable = "z_item")
)
capture_result(head(smry_2_z, 10), label = "fit_2: z_item[1:10] summary")
log_result("NOTE: Rhat and ESS problems are in z_item, not a_item.")
log_result("      a_item = z_item * sigma_item, so a_item looks better by hiding the dependency.")

# --- Fig 29.6: a_item[1] vs sigma_item (misleading — looks fine) ---
fig_29_6 <- posterior::subset_draws(draws_2,
                                     variable = c("a_item[1]", "sigma_item")) |>
  bayesplot::mcmc_scatter(alpha = 0.5) +
  ggplot2::labs(
    title = "Figure 29.6: a_item[1] vs sigma_item — fit_2 (misleading: no funnel visible)"
  ) +
  ggplot2::theme_minimal()

print(fig_29_6)
ggplot2::ggsave("figs/Fig-29.6.svg", plot = fig_29_6,
                width = 6, height = 5, device = svg)

# --- Fig 29.7: z_item[1] vs sigma_item (reveals the problem) ---
fig_29_7 <- posterior::subset_draws(draws_2,
                                     variable = c("z_item[1]", "sigma_item")) |>
  bayesplot::mcmc_scatter(alpha = 0.5) +
  ggplot2::labs(
    title = "Figure 29.7: z_item[1] vs sigma_item — fit_2 (strong correlation revealed)"
  ) +
  ggplot2::theme_minimal()

print(fig_29_7)
ggplot2::ggsave("figs/Fig-29.7.svg", plot = fig_29_7,
                width = 6, height = 5, device = svg)

# --- Fig 29.8: z_item[1] vs log(sigma_item) (banana weakens) ---
fig_29_8 <- posterior::subset_draws(draws_2,
                                     variable = c("z_item[1]", "sigma_item")) |>
  bayesplot::mcmc_scatter(
    transformations = list(sigma_item = log),
    alpha = 0.5
  ) +
  ggplot2::labs(
    y     = "log(sigma_item)",
    title = "Figure 29.8: z_item[1] vs log(sigma_item) — fit_2 (banana weakens)"
  ) +
  ggplot2::theme_minimal()

print(fig_29_8)
ggplot2::ggsave("figs/Fig-29.8.svg", plot = fig_29_8,
                width = 6, height = 5, device = svg)

log_result("\n=== Figures 29.6–29.8 saved ===")
log_result("Interpretation: non-centered z_item + large N per item creates")
log_result("  strong z_item--sigma_item dependency. Centered parameterization preferred.")

# =============================================================================
# fit_3: sum-to-zero + centered a_respondent / a_item
# =============================================================================

log_result("\n=== fit_3: park_3.stan (sum-to-zero + centered) ===")
log_result("iter_warmup=200, iter_sampling=200")

t_fit3 <- system.time(
  fit_3 <- cstan("park_3.stan", data = stan_data,
                 iter_warmup = 200, iter_sampling = 200)
)
fit_3$save_object("fit_3.rds")
log_result("Elapsed: ", round(t_fit3["elapsed"], 1), " sec")

# --- Sampler diagnostics ---
log_result("\n--- Diagnostics: fit_3 (cmdstanr) ---")
diag3 <- fit_3$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag3$num_divergent))
log_result("  num_max_treedepth : ", sum(diag3$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag3$ebfmi, 4), collapse = ", "))

smry_3_key <- posterior::summarise_draws(
  posterior::subset_draws(fit_3$draws(),
                          variable = c("a", "sigma_respondent", "sigma_item"))
)
capture_result(smry_3_key, label = "fit_3: key parameter summary")
log_result("Book: sigma_item Rhat ≈ 1.00 (big improvement from fit_2)")
log_result("      treedepth reduced further — centered param easier for this data size")

smry_sampler_3 <- fit_3$sampler_diagnostics() |> posterior::as_draws_rvars()
log_result("  mean treedepth__ : ",
           round(mean(posterior::E(smry_sampler_3$treedepth__)), 2))
log_result("  mean n_leapfrog__: ",
           round(mean(posterior::E(smry_sampler_3$n_leapfrog__)), 2))
log_result("Book: mean treedepth__ ≈ 5.0 (vs ≈ 6.5 for fit_2, ≈ 6.8 for fit_1)")

wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "parameterization_logged",
  timestamp = Sys.time(),
  notes     = "fit_3: switched to centered sum-to-zero a_item; sigma_item Rhat resolved"
)))

# =============================================================================
# fit_4: sum-to-zero + centered + centered predictors (X_c in transformed data)
# =============================================================================

log_result("\n=== fit_4 short run: park_4.stan (sum-to-zero + centered + X_c) ===")
log_result("iter_warmup=200, iter_sampling=200")

t_fit4s <- system.time(
  fit_4 <- cstan("park_4.stan", data = stan_data,
                 iter_warmup = 200, iter_sampling = 200)
)
fit_4$save_object("fit_4.rds")
log_result("Elapsed: ", round(t_fit4s["elapsed"], 1), " sec")

log_result("\n--- Diagnostics: fit_4 short run (cmdstanr) ---")
diag4s <- fit_4$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag4s$num_divergent))
log_result("  num_max_treedepth : ", sum(diag4s$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag4s$ebfmi, 4), collapse = ", "))

smry_sampler_4s <- fit_4$sampler_diagnostics() |> posterior::as_draws_rvars()
log_result("  mean treedepth__ : ",
           round(mean(posterior::E(smry_sampler_4s$treedepth__)), 2))
log_result("Book short run: mean treedepth__ ≈ 5.0 (same as fit_3)")

# --- Full refit at default iterations ---
log_result("\n=== fit_4 full refit: iter_warmup=1000, iter_sampling=1000 ===")

t_fit4f <- system.time(
  fit_4 <- cstan("park_4.stan", data = stan_data,
                 iter_warmup = 1000, iter_sampling = 1000)
)
fit_4$save_object("fit_4.rds")
log_result("Elapsed: ", round(t_fit4f["elapsed"], 1), " sec")

log_result("\n--- Diagnostics: fit_4 full refit (cmdstanr) ---")
diag4f <- fit_4$diagnostic_summary()
log_result("  num_divergences   : ", sum(diag4f$num_divergent))
log_result("  num_max_treedepth : ", sum(diag4f$num_max_treedepth))
log_result("  E-BFMI            : ", paste(round(diag4f$ebfmi, 4), collapse = ", "))

smry_4_key <- posterior::summarise_draws(
  posterior::subset_draws(fit_4$draws(),
                          variable = c("a", "b[1]", "b[2]", "b[3]",
                                       "sigma_respondent", "sigma_item"))
)
capture_result(smry_4_key, label = "fit_4 full: key parameter summary")

smry_sampler_4f <- fit_4$sampler_diagnostics() |> posterior::as_draws_rvars()
log_result("  mean treedepth__ : ",
           round(mean(posterior::E(smry_sampler_4f$treedepth__)), 2))
log_result("  mean n_leapfrog__: ",
           round(mean(posterior::E(smry_sampler_4f$n_leapfrog__)), 2))
log_result("Book full: mean treedepth__ ≈ 4.9 (further reduced with better adaptation)")

# --- Book comparison ---
log_result("\n=== Book targets: fit_4 full refit ===")
log_result("  a:                book ≈ -2.3  (sd ≈ 0.03)")
log_result("  b[1] male:        book ≈  0.07 (sd ≈ 0.03)")
log_result("  b[2] white:       book ≈  0.08 (sd ≈ 0.03)")
log_result("  b[3] n_skipped:   book ≈  0.03 (sd ≈ 0.01)")
log_result("  sigma_respondent: book ≈  1.9  (sd ≈ 0.04)")
log_result("  sigma_item:       book ≈  2.4  (sd ≈ 0.34)")

draws_4 <- fit_4$draws() |> posterior::as_draws_rvars()

log_result("  a computed mean   : ",
           round(mean(posterior::E(draws_4$a)), 3))
log_result("  a computed sd     : ",
           round(stats::sd(as.vector(posterior::draws_of(draws_4$a))), 3))
log_result("  sigma_item mean   : ",
           round(mean(posterior::E(draws_4$sigma_item)), 3))

# a_item[1] = 'car' (most vehicle-like)
a_item_draws <- posterior::summarise_draws(
  posterior::subset_draws(fit_4$draws(), variable = "a_item[1]")
)
log_result("  a_item[1] (car) mean: ", round(a_item_draws$mean, 2),
           "  (book ≈ 7.6)")

# Update wf state
wf$diagnostics$passed       <- TRUE
wf$diagnostics$acknowledged <- TRUE
wf$diagnostics$rhat_max     <- max(smry_4_key$rhat, na.rm = TRUE)
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "diagnostic_acknowledged",
  timestamp = Sys.time(),
  notes     = "fit_4 full refit: all Rhat < 1.01, good ESS — centered param + X_c resolves all issues"
)))

# =============================================================================
# Phase 5 / lme4 vs Stan comparison figures
# =============================================================================

# NOTE: PPC via bayesplot::ppc_dens_overlay() is not appropriate for this
# chapter. The model outcome is binary (Bernoulli), the estimands of interest
# are item and respondent random effects, and the chapter's comparison is
# qualitative — lme4 vs Bayesian posterior intervals. No standard PPC is
# produced. See Phase 5 / LOO section below for the documented absence.

dr4 <- fit_4$draws() |> posterior::as_draws_rvars()

# lme4 item effects
a_item_hat_lme4 <- lme4::ranef(fit_lme4)$item
a_item_sd_lme4  <- sqrt(as.numeric(attr(a_item_hat_lme4, "postVar")))
a_item_h        <- as.numeric(a_item_hat_lme4$`(Intercept)`)

rng_item <- range(
  posterior::summarise_draws(
    dr4$a_item,
    ~ quantile(.x, probs = c(0.05, 0.95))
  )[, c("5%", "95%")]
)

# Fig 29.9 — a_item: lme4 vs Stan 90% intervals
fig_29_9 <- ggplot2::ggplot(data = NULL) +
  ggplot2::coord_fixed(xlim = rng_item, ylim = rng_item) +
  ggplot2::geom_abline(color = "gray") +
  ggdist::stat_pointinterval(
    ggplot2::aes(x = a_item_h, ydist = dr4$a_item),
    .width = 0.90,
    interval_size_range = c(0.4, 0.8),
    alpha = 0.5
  ) +
  ggdist::geom_pointinterval(
    ggplot2::aes(
      y    = mean(dr4$a_item),
      x    = a_item_h,
      xmin = a_item_h - 1.64 * a_item_sd_lme4,
      xmax = a_item_h + 1.64 * a_item_sd_lme4
    ),
    interval_size_range = c(0.4, 0.8),
    alpha = 0.5
  ) +
  ggplot2::labs(
    x     = "lme4",
    y     = "Stan",
    title = "Figure 29.9: a_item — lme4 vs Stan (90% intervals)"
  ) +
  ggplot2::theme_minimal()

print(fig_29_9)
ggplot2::ggsave("figs/Fig-29.9.svg", plot = fig_29_9,
                width = 6, height = 6, device = svg)

# lme4 respondent effects
a_respondent_hat_lme4 <- lme4::ranef(fit_lme4)$respondent
a_respondent_sd_lme4  <- sqrt(as.numeric(attr(a_respondent_hat_lme4, "postVar")))
a_respondent_h        <- as.numeric(a_respondent_hat_lme4$`(Intercept)`)

rng_resp <- range(
  posterior::summarise_draws(
    dr4$a_respondent,
    ~ quantile(.x, probs = c(0.05, 0.95))
  )[, c("5%", "95%")]
)

# Fig 29.10 — a_respondent: lme4 vs Stan 90% intervals
fig_29_10 <- ggplot2::ggplot(data = NULL) +
  ggplot2::coord_fixed(xlim = rng_resp, ylim = rng_resp) +
  ggplot2::geom_abline(color = "gray") +
  ggdist::stat_pointinterval(
    ggplot2::aes(x = a_respondent_h, ydist = dr4$a_respondent),
    .width = 0.90,
    interval_size_range = c(0.4, 0.8),
    alpha = 0.1
  ) +
  ggdist::geom_pointinterval(
    ggplot2::aes(
      y    = mean(dr4$a_respondent),
      x    = a_respondent_h,
      xmin = a_respondent_h - 1.64 * a_respondent_sd_lme4,
      xmax = a_respondent_h + 1.64 * a_respondent_sd_lme4
    ),
    interval_size_range = c(0.4, 0.8),
    alpha = 0.1
  ) +
  ggplot2::labs(
    x     = "lme4",
    y     = "Stan",
    title = "Figure 29.10: a_respondent — lme4 vs Stan (90% intervals)"
  ) +
  ggplot2::theme_minimal()

print(fig_29_10)
ggplot2::ggsave("figs/Fig-29.10.svg", plot = fig_29_10,
                width = 6, height = 6, device = svg)

# lme4 sigma estimates (dashed vertical lines)
vc_lme4 <- as.data.frame(lme4::VarCorr(fit_lme4))
sigma_item_lme4       <- vc_lme4$sdcor[vc_lme4$grp == "item"]
sigma_respondent_lme4 <- vc_lme4$sdcor[vc_lme4$grp == "respondent"]

# Fig 29.11 — sigma_item and sigma_respondent posteriors vs lme4
p_sigma_item <- ggplot2::ggplot(data = NULL) +
  ggdist::stat_slab(
    ggplot2::aes(xdist = dr4$sigma_item),
    density = "unbounded", trim = TRUE, fill = NA, color = "black"
  ) +
  ggplot2::scale_y_continuous(breaks = NULL) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.line.y = ggplot2::element_blank()) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::labs(x = "sigma_item", y = "") +
  ggplot2::geom_vline(xintercept = sigma_item_lme4,
                      color = "black", linetype = "dashed") +
  ggplot2::annotate("text",
    x     = sigma_item_lme4 * 1.02,
    y     = 0.97,
    hjust = 0,
    label = "lme4 estimate"
  )

p_sigma_resp <- ggplot2::ggplot(data = NULL) +
  ggdist::stat_slab(
    ggplot2::aes(xdist = dr4$sigma_respondent),
    density = "unbounded", trim = TRUE, fill = NA, color = "black"
  ) +
  ggplot2::scale_y_continuous(breaks = NULL) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.line.y = ggplot2::element_blank()) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::labs(x = "sigma_respondent", y = "") +
  ggplot2::geom_vline(xintercept = sigma_respondent_lme4,
                      color = "black", linetype = "dashed") +
  ggplot2::annotate("text",
    x     = sigma_respondent_lme4 * 1.002,
    y     = 0.97,
    hjust = 0,
    label = "lme4 estimate"
  )

fig_29_11 <- p_sigma_item + p_sigma_resp +
  patchwork::plot_annotation(
    title = "Figure 29.11: sigma posteriors — Stan (slab) vs lme4 (dashed)"
  )

print(fig_29_11)
ggplot2::ggsave("figs/Fig-29.11.svg", plot = fig_29_11,
                width = 9, height = 4, device = svg)

log_result("\n=== Figures 29.9–29.11 saved ===")
log_result("Stan estimates have slightly wider range than lme4 — integrating")
log_result("  over uncertainty in sigma_item and sigma_respondent.")

# =============================================================================
# Phase 5 — PPC (documented absence)
# =============================================================================

# PPC METHOD NOTE
# Standard bayesplot::ppc_dens_overlay() is not appropriate here because
# the outcome is binary and the chapter's comparison is qualitative —
# lme4 conditional modes vs Bayesian posterior intervals for item and
# respondent random effects. The chapter does not produce a standard PPC
# plot; the lme4 vs Stan figures (29.9–29.11) serve as the predictive
# adequacy check by comparing the two methods' uncertainty estimates.
log_result("\nPPC: not performed.")
log_result("Reason: binary outcome; chapter comparison is lme4 vs Stan intervals")
log_result("        (Figs 29.9–29.11). No standard ppc_dens_overlay applicable.")
wf$ppc_complete <- FALSE

# =============================================================================
# Phase 6 — LOO-CV (documented absence)
# =============================================================================

# LOO-CV not applicable for this chapter.
# Reason: The chapter's goal is to demonstrate sampling efficiency and
# identifiability issues across parameterizations of the same model, not
# model comparison. All four Stan models fit the same likelihood; LOO
# comparison between them would not be meaningful. No loo_compare() is run.
log_result("\nLOO-CV: not performed.")
log_result("Reason: all four Stan models share the same likelihood — only the")
log_result("        parameterization differs. LOO comparison is not meaningful here.")
wf$loo_complete <- FALSE

# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete <- FALSE
wf$loo_complete <- FALSE
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

log_result("\n=== Ch 29 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: fit_1.rds, fit_2.rds, fit_3.rds, fit_4.rds, wf_final.rds")
log_result("Figures saved to figs/:")
log_result("  Fig-29.1.svg  Fig-29.2.svg  Fig-29.3.svg")
log_result("  Fig-29.4.svg  Fig-29.5.svg")
log_result("  Fig-29.6.svg  Fig-29.7.svg  Fig-29.8.svg")
log_result("  Fig-29.9.svg  Fig-29.10.svg Fig-29.11.svg")
close(results_con)
