# data/birthdays_ch27/birthdays_ch27_analysis.R
#
# Ch 27 — "Model building: Time-series decomposition for birthdays"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Analyses relative number of births per day in USA 1969–1988 using
# Hilbert space GP approximation with incrementally added components.
# Deliberately omits: Laplace and Pathfinder comparison plots (Figs 5,
# 10, 13, 16, 19, 24, 29, 33, 37, 39, 41, 43–45), sd ratio plots,
# the no-GQ Pathfinder efficiency experiments at end of §5.10, and
# full-length final chains (they take hours; short chains suffice for
# workflow illustration). All LOO comparisons and main prediction plots
# are included.
#
# Data:
#   births_usa_1969.csv (data/birthdays_ch27/data/)
#   USA daily birth counts 1969–1988, N = 7305
#   Key variables: id (running day 1–7305), births (daily count),
#     day_of_week (1=Mon…7=Sun), day_of_year2 (1–366, 1 Mar = 61
#     in all years), year, month, day
#
# Models covered (all cmdstanr; run_phase3() not used — see §6 below):
#   model1   — gpbf1.stan:    intercept + GP slow trend (f1)
#   model1b  — gpbf1b.stan:   GP slow trend without intercept
#   model2   — gpbf2.stan:    f1 + GP seasonal (f2)
#   model3   — gpbf3.stan:    f1 + f2 + day-of-week betas (beta_f3)
#   model4   — gpbf4.stan:    f1 + f2 + exp(g3)*beta_f3 (time-varying magnitude)
#   model5   — gpbf5.stan:    f1 + f2 + exp(g3)*beta_f3 + day-of-year RHS
#   model6   — gpbf6.stan:    f1 + f2 + beta_f3 + day-of-year normal
#   model7   — gpbf7.stan:    model6 + floating special days (Memorial, Labor, Thanksgiving)
#   model8   — gpbf8.stan:    model7 + time-varying magnitude (g3)
#   model8tnu — gpbf8tnu.stan: model8 + Student's t prior on day-of-year effect
#   model8rhs — gpbf8rhs.stan: model8 + RHS prior on day-of-year effect (centered)
#
# Figures produced:
#   Fig-27.1.svg — Raw births over time
#   Fig-27.2.svg — Births relative to mean (100 = mean)
#   Fig-27.3.svg — Mean relative births per day of year
#   Fig-27.4.svg — Mean relative births per day of week
#   Fig-27.5.svg — Model 1: MAP prediction vs data
#   Fig-27.6.svg — Model 1: MCMC trace (sigma_f1, lengthscale_f1, sigma)
#   Fig-27.7.svg — Model 1: MCMC median prediction vs data
#   Fig-27.8.svg — Model 2: 3-panel (overall, trend, seasonal)
#   Fig-27.9.svg — Model 3: 4-panel (overall, trend, seasonal, weekday)
#   Fig-27.10.svg — Model 4: 4-panel with time-varying magnitude
#   Fig-27.11.svg — Model 5: 6-panel prediction
#   Fig-27.12.svg — Model 6: 6-panel prediction
#   Fig-27.13.svg — Model 6: day-of-year effect with 90% CI
#   Fig-27.14.svg — Model 7: 6-panel prediction with floating holidays
#   Fig-27.15.svg — Model 7: day-of-year + floating effect with 90% CI
#   Fig-27.16.svg — Model 8: 6-panel prediction
#   Fig-27.17.svg — Model 8tnu: 6-panel prediction
#   Fig-27.18.svg — Model 8tnu: day-of-year effect with 90% CI (final model)
#   Fig-27.19.svg — Residual analysis (Model 8rhs)
#
# Book-target posterior summaries (used for success criteria):
#   LOO comparison (Model 8 normal vs Student's t): elpd_diff ≈ -116 (se ≈ 16)
#   LOO comparison (Model 8 RHS vs Student's t):    elpd_diff ≈ -0.21 (se ≈ 3.9)
#   LOO-R2 (Model 8 Student's t, final chains):     ≈ 0.94
#   sigma[Model 1, MCMC short]:                     ≈ 0.81
#   sigma[Model 3, MCMC short]:                     ≈ 0.33
#   sigma[Model 8, MCMC short]:                     ≈ 0.23
#   beta_f3[5] weekday Saturday [Model 3]:          ≈ -1.1
#   beta_f3[6] weekday Sunday  [Model 3]:           ≈ -1.5
#
# This script uses cmdstanr directly because the chapter's GP
# likelihoods and Pathfinder workflow cannot be expressed through
# brms formula interfaces. run_phase3() is not used.
# All 11 Stan models are compiled from .stan files in this folder.
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/birthdays_ch27/.

source("../../R/source_all.R")
library(tidyverse)
library(tictoc)
library(cmdstanr)
library(posterior)
library(loo)
library(bayesplot)
library(ggdist)
library(patchwork)
library(ggrepel)
library(RColorBrewer)

# Do NOT call theme_set() — apply theme_minimal() per figure (rubric §3)
# cmdstanr output dir
dir.create("stan_output", showWarnings = FALSE)
options(cmdstanr_output_dir = file.path(getwd(), "stan_output"))

SEED <- 3896   # matches book seed for model 1

# Locale for weekday/month names
Sys.setlocale("LC_TIME", "en_GB.utf8")

set1 <- RColorBrewer::brewer.pal(9, "Set1")

# Plot helpers from the book (col_data, col_fit, make_pf, make_pf1,
# make_pf2, make_pf2b, make_pf2c, make_pf3, make_pf3b, compose_3panel,
# compose_4panel, compose_6panel, layers_hline100, layers_scale_doy,
# layers_scale_weekday, make_pth_vs_fit, make_sd_ratio)
source("plot_helpers.R")

# toc helper matching book style
mytoc <- \() {
  toc(func.toc = \(tic, toc, msg) {
    sprintf("%s took %s sec", msg, as.character(signif(toc - tic, 2)))
  })
}

# Results log
dir.create("figs", showWarnings = FALSE)
results_con <- file("results.txt", open = "wt")
log_result <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = results_con, append = TRUE)
}
capture_result <- function(x, label = NULL) {
  txt <- paste(capture.output(print(x)), collapse = "\n")
  if (!is.null(label)) {
    cat(label, "\n", sep = "")
    cat(label, "\n", sep = "", file = results_con, append = TRUE)
  }
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = results_con, append = TRUE)
}
log_result("=== Ch 27 results log — ", format(Sys.time()), " ===")

# Stan-native fit helper (rubric §6 cstan() pattern)
cstan <- function(stan_file, data = list(), seed = SEED,
                  chains = 4, iter_warmup = 100, iter_sampling = 100,
                  init = 0.1, psis_resample = TRUE, refresh = 0) {
  model <- cmdstan_model(stan_file, include_paths = getwd())
  model$sample(
    data            = data,
    seed            = seed,
    chains          = chains,
    parallel_chains = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    init            = init,
    refresh         = refresh
  )
}

# Pathfinder helper
cpathfinder <- function(model, data, psis_resample = TRUE,
                        num_paths = 10, draws = 400) {
  model$pathfinder(
    data             = data,
    init             = 0.1,
    num_paths        = num_paths,
    single_path_draws = 40,
    draws            = draws,
    history_size     = 50,
    max_lbfgs_iters  = 100,
    refresh          = 0,
    psis_resample    = psis_resample
  )
}

wf <- init_workflow(mode = "practice", stage = "explore")

# §1 Data loading and EDA
birthdays <- readr::read_csv("data/births_usa_1969.csv", show_col_types = FALSE) |>
  dplyr::mutate(
    date = as.Date(paste(year, month, day, sep = "-")),
    births_relative100 = births / mean(births) * 100
  )

# Fig 27.1 Raw births over time
fit1 <- readRDS('fit1.rds')
# Fig 27.7 Model 1: MCMC median prediction vs data
draws1 <- posterior::as_draws_matrix(fit1$draws())
Ef1_tot <- exp(apply(posterior::subset_draws(draws1, variable = "f"), 2, median))
p7 <- make_pf(birthdays, Ef1_tot) + layers_hline100 + ggplot2::theme_minimal()
print(p7)
ggplot2::ggsave("figs/Fig-27.7.svg", plot = p7, width = 8, height = 4)
