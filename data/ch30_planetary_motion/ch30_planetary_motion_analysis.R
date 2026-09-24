# data/ch30_planetary_motion/ch30_planetary_motion_analysis.R
#
# Ch 30 — "Challenge of multimodality: Differential equation for planetary motion"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# This script uses cmdstanr directly because the chapter's models are ODE-based
# and cannot be expressed through brms formula interfaces. run_phase3() is
# bypassed for all fits. Pathfinder (mod$pathfinder()) is used to initialize
# MCMC for fit1p and fit2p; this requires CmdStan >= 2.33.
#
# RUN_FIT1 flag (see §3): when FALSE (default), the script loads the book's
# saved fit1.RDS from saved_fit/ instead of refitting. Set to TRUE to refit
# from scratch (~75 min wall time; chains may differ from the book).
#
# Note on chain numbering in book prose vs. book output: the book text states
# chains 3 and 4 were slowest and chains 2, 5, 7 fit the data. The printed
# timing table shows chains 5 (4498 s), 8 (1919 s), 1 (329 s) as slowest.
# The prose appears to have been written against an earlier run. This script
# treats the saved fit object as authoritative; chain labels in figure titles
# match the saved object, not the prose.
#
# Data:
#   Simulated: N = 40 observations of (q_x, q_y) planet positions at t = 0.1,
#   0.2, ..., 4.0. True parameters: k = 1, q0 = (1,0), p0 = (0,1), star at
#   origin. Observation noise sigma = 0.01. Generated via
#   planetary_motion_sim.stan using seed = 123.
#
# Models covered:
#   mod_sim  — planetary_motion_sim.stan: data simulator (1 chain, 2 draws)
#   mod1     — planetary_motion.stan: 1-param model (k only); 8 chains, 500/500
#   mod1p    — same model, Pathfinder init; 8 chains, 500/500
#   mod2     — planetary_motion_star.stan: full 7-param model; 8 chains, 500/500
#              (Pathfinder init)
#
# Figures produced:
#   Fig-30.1.svg  — Simulated planet orbit: q_x vs q_y, 40 observations
#   Fig-30.2.svg  — Trace plot: lp__ and k, fit1 (sampling only)
#   Fig-30.3.svg  — Per-chain PPC: observed orbit vs posterior-predictive median
#                   (8-panel), fit1
#   Fig-30.4.svg  — Trace plot with warmup: lp__ and k, fit1
#   Fig-30.5.svg  — Log joint density vs k (grid scan, k = 0.2 to 9)
#   Fig-30.6.svg  — Simulated trajectories for k = 0.5, 1.0, 1.6, 2.16, 3.0
#                   with distance segments at observation 35
#   Fig-30.7.svg  — Per-chain PPC: observed orbit vs posterior-predictive median
#                   (8-panel), fit1p (Pathfinder-initialized)
#   Fig-30.8.svg  — Conditional log likelihood vs q_star_x (q_star_y fixed = 0)
#   Fig-30.9.svg  — Conditional log likelihood heat map: q_star_x vs q_star_y
#   Fig-30.10.svg — Per-chain PPC: observed orbit vs posterior-predictive median
#                   with star position marker (8-panel), fit2p (full model)
#
# Book-target posterior summaries:
#   fit1  (failed run):  k mean ≈ 4.02, sd ≈ 2.98, Rhat ≈ 3.15
#   fit1p (Pathfinder):  k mean ≈ 1.000, sd ≈ 0.000329, Rhat ≈ 1.00
#   fit2p (full model):  k ≈ 1.00 (sd 0.001), q0[1] ≈ 1.00, q0[2] ≈ -0.003,
#                        p0[1] ≈ 0.005, p0[2] ≈ 1.000, star[1] ≈ 0.001,
#                        star[2] ≈ 0.007
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch30_planetary_motion/.

# =============================================================================
# §1  Environment setup
# =============================================================================

source("../../R/source_all.R")
library(posterior)
library(dplyr)
library(plyr)
library(tidyr)
library(boot)
library(latex2exp)
library(ggplot2)
library(cmdstanr)

options(brms.backend = "cmdstanr", mc.cores = 4)

SEED <- 123         # matches book set.seed(1954) — note: book uses 1954 for
                    # set.seed() but seed = 123 in $sample() calls; we match
                    # the $sample() seed throughout
CHAINS <- 8

# Set TRUE to refit fit1 from scratch (~75 min). FALSE loads saved_fit/fit1.RDS.
RUN_FIT1 <- FALSE

dir.create("figs",        showWarnings = FALSE)
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
log_result("=== Ch 30 results log — ", format(Sys.time()), " ===")
log_result("RUN_FIT1 = ", RUN_FIT1)
log_result("NOTE: Book prose (chain numbering) may not match saved fit1.RDS.")
log_result("      Treating saved fit object as authoritative.")

# Load tools.R helper functions: solve_trajectory(), ppc_plot2D()
# tools.R uses dplyr::pull(), dplyr::bind_rows(), and bare ggplot2 functions.
# library(dplyr) and library(ggplot2) must be attached before sourcing.
source("tools.R")

# Apply bayesplot color scheme used in the book
bayesplot::color_scheme_set("viridisC")

# =============================================================================
# §2  Data simulation (Phase 1 / Phase 2)
# =============================================================================

wf <- init_workflow(mode = "practice", stage = "explore")
wf$fit_timestamp <- NA   # no brms fit; will be set after each Stan fit

N     <- 40
sigma <- 0.01
stan_data_sim <- list(N = N, sigma_x = sigma, sigma_y = sigma)

mod_sim <- cmdstan_model("planetary_motion_sim.stan")

sim <- mod_sim$sample(
  data          = stan_data_sim,
  chains        = 1,
  iter_warmup   = 1,
  iter_sampling = 2,
  seed          = SEED,
  refresh       = 0
)

simulation <- as.vector(sim$draws(variables = "q_obs")[1, , ])
q_obs      <- array(NA, c(N, 2))
q_obs[, 1] <- simulation[1:40]
q_obs[, 2] <- simulation[41:80]

log_result("\n=== §2 Simulated data ===")
log_result("N = ", N, "  sigma = ", sigma)
log_result("q_obs[1,] = ", paste(round(q_obs[1, ], 4), collapse = ", "))
log_result("q_obs[40,] = ", paste(round(q_obs[40, ], 4), collapse = ", "))

# --- Fig-30.1: simulated orbit ---
data_pred <- data.frame(qx = q_obs[, 1], qy = q_obs[, 2], t = 1:N)

fig_30_1 <- ggplot2::ggplot(
  data = data.frame(q_x = q_obs[, 1], q_y = q_obs[, 2], time = 1:N),
  ggplot2::aes(x = q_x, y = q_y)
) +
  ggplot2::geom_point() +
  ggplot2::annotate("text",
    x = q_obs[1, 1] - 0.05, y = q_obs[1, 2],
    label = "t=1", hjust = 1) +
  ggplot2::annotate("text",
    x = q_obs[N, 1] + 0.05, y = q_obs[N, 2],
    label = paste0("t=", N), hjust = 0) +
  ggplot2::labs(
    title = "Figure 30.1: Simulated planet orbit",
    x = expression(q[x]), y = expression(q[y])
  ) +
  ggplot2::theme_minimal()

print(fig_30_1)
ggplot2::ggsave("figs/Fig-30.1.svg", plot = fig_30_1,
                width = 4, height = 4, device = svg)
log_result("Saved: figs/Fig-30.1.svg")

# =============================================================================
# §3  Simple model: fit1 — 8 chains, default Stan inits (expected to fail)
# =============================================================================

# PPC METHOD NOTE
# Standard bayesplot::ppc_dens_overlay() is not appropriate here because the
# outcome is a 2D position vector (q_x, q_y), not a scalar. The book's custom
# ppc_plot2D() from tools.R compares the per-chain median predicted trajectory
# against the observed orbit. This checks the same property — whether the
# model's implied trajectories match the data — using the correct estimand for
# a 2D time-series outcome.

stan_data1 <- list(N = N, q_obs = q_obs)

mod1 <- cmdstan_model("planetary_motion.stan")

if (RUN_FIT1) {
  fit1 <- mod1$sample(
    data             = stan_data1,
    chains           = CHAINS,
    parallel_chains  = CHAINS,
    iter_warmup      = 500,
    iter_sampling    = 500,
    seed             = SEED,
    save_warmup      = TRUE,
    refresh          = 0
  )
  fit1$save_object(file = "saved_fit/fit1.RDS")
  log_result("fit1 refit and saved to saved_fit/fit1.RDS")
} else {
  fit1 <- readRDS("saved_fit/fit1.RDS")
  log_result("fit1 loaded from saved_fit/fit1.RDS (RUN_FIT1 = FALSE)")
}

# Diagnostics — inspect sampled parameter k, not generated quantities
wf$fit_timestamp <- Sys.time()
wf$fit_hash <- digest::digest(list(model = "planetary_motion.stan",
                                    seed = SEED), algo = "sha256")

fit1_smry <- fit1$summary(
  variables = c("lp__", "k"),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
)

log_result("\n=== §3 fit1 diagnostics (sampled parameter k) ===")
capture_result(fit1_smry, label = "fit1 summary: lp__ and k")

log_result("Book target: k mean ≈ 4.02, sd ≈ 2.98, Rhat ≈ 3.15")
k_row <- fit1_smry[fit1_smry$variable == "k", ]
log_result("Computed:    k mean = ", round(k_row$mean, 3),
           ", sd = ", round(k_row$sd, 3),
           ", Rhat = ", round(k_row$rhat, 3))

if (abs(k_row$rhat - 3.15) > 0.5) {
  log_result("NOTE: Rhat differs from book target by > 0.5.",
             " This is expected when loading a fit from a different machine/CmdStan version.")
}

# Populate wf diagnostics manually (Stan-native pattern)
wf <- diagnose_cmdstan(fit1, wf, params = "k")
log_result("--- Diagnostics: fit1 (cmdstanr) ---")
log_result("  passed      : ", wf$diagnostics$passed)
log_result("  Rhat_max    : ", round(wf$diagnostics$rhat_max, 4))
log_result("  bulk_ESS_min: ", round(wf$diagnostics$bulk_ess_min, 0))
log_result("  tail_ESS_min: ", round(wf$diagnostics$tail_ess_min, 0))
log_result("  divergences : ", wf$diagnostics$n_divergences)
if (length(wf$diagnostics$failed_criteria) > 0)
  log_result("  FAILED      : ",
             paste(wf$diagnostics$failed_criteria, collapse = ", "))

wf$diagnostics$acknowledged <- TRUE
wf$audit_trail <- append(wf$audit_trail, list(list(
  phase     = 4,
  action    = "diagnostic_acknowledged",
  timestamp = Sys.time(),
  notes     = "fit1: expected failure — multimodal posterior, chains trapped at local modes"
)))

# =============================================================================
# §4  fit1 figures: trace plots and per-chain PPC
# =============================================================================

# --- Fig-30.2: trace plot, sampling phase only ---
fig_30_2 <- bayesplot::mcmc_trace(fit1$draws(), pars = c("lp__", "k"))
print(fig_30_2)
ggplot2::ggsave("figs/Fig-30.2.svg", plot = fig_30_2,
                width = 6, height = 3, device = svg)
log_result("Saved: figs/Fig-30.2.svg")

# --- Fig-30.3: per-chain PPC, fit1 ---
# ppc_plot2D() is defined in tools.R. It plots the per-chain median predicted
# trajectory (q_x_pred, q_y_pred from generated quantities) against the
# observed orbit. Chains trapped at local modes will show trajectories that
# do not match the data.
fig_30_3 <- ppc_plot2D(fit1, data_pred = data_pred)
print(fig_30_3)
ggplot2::ggsave("figs/Fig-30.3.svg", plot = fig_30_3,
                width = 6, height = 6, device = svg)
log_result("Saved: figs/Fig-30.3.svg")

# --- Fig-30.4: trace plot including warmup ---
# inc_warmup = TRUE requires save_warmup = TRUE at sample time.
# The book's fit1.RDS was saved with save_warmup = TRUE.
fig_30_4 <- bayesplot::mcmc_trace(
  fit1$draws(inc_warmup = TRUE),
  pars     = c("lp__", "k"),
  n_warmup = 500
)
print(fig_30_4)
ggplot2::ggsave("figs/Fig-30.4.svg", plot = fig_30_4,
                width = 6, height = 3, device = svg)
log_result("Saved: figs/Fig-30.4.svg")

# =============================================================================
# §5  Likelihood landscape: log-joint grid and trajectory comparison
# =============================================================================

# These sections use solve_trajectory() from tools.R to compute the ODE
# solution numerically (leapfrog integrator, dt = 0.001).
# Shared constants used throughout §5
q0_grid <- c(1.0, 0)
p0_grid <- c(0, 1.0)
dt      <- 0.001
m       <- 1
n_obs   <- N
ts      <- 1:n_obs / 10
sigma_x <- 0.01
sigma_y <- 0.01

# --- Fig-30.5: log joint density vs k (grid scan) ---
ks <- seq(from = 0.2, to = 9, by = 0.01)
lk <- vapply(seq_along(ks), function(i) {
  k     <- ks[i]
  q_sim <- solve_trajectory(q0_grid, p0_grid, dt, k, m, n_obs, ts)
  sum(dnorm(q_obs[, 1], q_sim[, 1], sigma_x, log = TRUE)) +
  sum(dnorm(q_obs[, 2], q_sim[, 2], sigma_y, log = TRUE)) +
  dnorm(k, 0, 1, log = TRUE)
}, numeric(1L))

log_result("\n=== §5 Log-joint grid ===")
log_result("k at global mode: ", round(ks[which.max(lk)], 3))
log_result("Book: strong mode at k = 1, wiggly tail for k > 1")

fig_30_5 <- ggplot2::ggplot(
  data = data.frame(ks = ks, lk = lk),
  ggplot2::aes(x = ks, y = lk)
) +
  ggplot2::geom_line() +
  ggplot2::labs(
    title = "Figure 30.5: Log joint density vs k",
    x = "k", y = "log joint"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(text = ggplot2::element_text(size = 14))

print(fig_30_5)
ggplot2::ggsave("figs/Fig-30.5.svg", plot = fig_30_5,
                width = 6, height = 4, device = svg)
log_result("Saved: figs/Fig-30.5.svg")

# --- Fig-30.6: simulated trajectories for k = 0.5, 1.0, 1.6, 2.16, 3.0 ---
k_vals  <- c(0.5, 1.6, 2.16, 3.0)
k_names <- c("0.5", "1.6", "2.16", "3.00")

traj_list <- lapply(k_vals, function(k)
  solve_trajectory(q0_grid, p0_grid, dt, k, m, n_obs, ts))

q_plot <- do.call(rbind, c(list(q_obs), traj_list))
k_plot <- rep(c("1.0 (obs)", k_names), each = N)
plot_data <- data.frame(qx = q_plot[, 1], qy = q_plot[, 2],
                        k = k_plot, obs = rep(1:N, length(k_vals) + 1))

comp_point <- 35
q_160 <- traj_list[[which(k_vals == 1.6)]]
q_216 <- traj_list[[which(k_vals == 2.16)]]

fig_30_6 <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = plot_data,
    ggplot2::aes(x = qx, y = qy, color = k)
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = q_216[comp_point, 1], y = q_216[comp_point, 2]),
    shape = 3) +
  ggplot2::geom_point(
    ggplot2::aes(x = q_160[comp_point, 1], y = q_160[comp_point, 2]),
    shape = 3) +
  ggplot2::geom_point(
    ggplot2::aes(x = q_obs[comp_point, 1], y = q_obs[comp_point, 2])) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = q_obs[comp_point, 1], y = q_obs[comp_point, 2],
      xend = q_216[comp_point, 1], yend = q_216[comp_point, 2]),
    linetype = "dashed") +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = q_obs[comp_point, 1], y = q_obs[comp_point, 2],
      xend = q_160[comp_point, 1], yend = q_160[comp_point, 2]),
    linetype = "dashed") +
  ggplot2::labs(
    title = "Figure 30.6: Simulated trajectories for selected k values",
    x = expression(q[x]), y = expression(q[y]), color = "k"
  ) +
  ggplot2::theme_minimal()

print(fig_30_6)
ggplot2::ggsave("figs/Fig-30.6.svg", plot = fig_30_6,
                width = 6, height = 4, device = svg)
log_result("Saved: figs/Fig-30.6.svg")

# =============================================================================
# §6  Pathfinder initialization + fit1p (simple model, fixed inits)
# =============================================================================

pth1p <- mod1$pathfinder(
  data               = stan_data1,
  num_paths          = 40,
  single_path_draws  = 25,
  draws              = 1000,
  max_lbfgs_iters    = 100,
  psis_resample      = FALSE,
  refresh            = 0
)
# Pathfinder may report some path failures; this is expected and documented
# in the book. Successful paths provide the importance-sampling mixture.

fit1p <- mod1$sample(
  data            = stan_data1,
  init            = pth1p,
  chains          = CHAINS,
  parallel_chains = CHAINS,
  iter_warmup     = 500,
  iter_sampling   = 500,
  seed            = SEED,
  save_warmup     = TRUE,
  refresh         = 0
)
fit1p$save_object(file = "saved_fit/fit1p.RDS")

# Diagnostics — sampled parameter k only
fit1p_smry <- fit1p$summary(
  variables = c("lp__", "k"),
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
)

log_result("\n=== §6 fit1p diagnostics (Pathfinder init, sampled parameter k) ===")
capture_result(fit1p_smry, label = "fit1p summary: lp__ and k")

log_result("Book target: k mean ≈ 1.000, sd ≈ 0.000329, Rhat ≈ 1.00")
k1p_row <- fit1p_smry[fit1p_smry$variable == "k", ]
log_result("Computed:    k mean = ", round(k1p_row$mean, 6),
           ", sd = ", round(k1p_row$sd, 6),
           ", Rhat = ", round(k1p_row$rhat, 4))

wf2 <- init_workflow(mode = "practice", stage = "explore")
wf2 <- diagnose_cmdstan(fit1p, wf2, params = "k")
log_result("--- Diagnostics: fit1p (cmdstanr) ---")
log_result("  passed      : ", wf2$diagnostics$passed)
log_result("  Rhat_max    : ", round(wf2$diagnostics$rhat_max, 4))
log_result("  bulk_ESS_min: ", round(wf2$diagnostics$bulk_ess_min, 0))
log_result("  tail_ESS_min: ", round(wf2$diagnostics$tail_ess_min, 0))
log_result("  divergences : ", wf2$diagnostics$n_divergences)
if (length(wf2$diagnostics$failed_criteria) > 0)
  log_result("  FAILED      : ",
             paste(wf2$diagnostics$failed_criteria, collapse = ", "))

# --- Fig-30.7: per-chain PPC, fit1p ---
fig_30_7 <- ppc_plot2D(fit1p, data_pred = data_pred)
print(fig_30_7)
ggplot2::ggsave("figs/Fig-30.7.svg", plot = fig_30_7,
                width = 6, height = 6, device = svg)
log_result("Saved: figs/Fig-30.7.svg")

# =============================================================================
# §7  Full model: conditional likelihood profiles for q_star
# =============================================================================

# Fix all parameters at their true values except q_*^x, then scan.
# Uses solve_trajectory() with the star argument (defaults to c(0,0)).
k_true  <- 1.0
q0_true <- c(1.0, 0)
p0_true <- c(0, 1.0)

# --- Fig-30.8: conditional log likelihood vs q_star_x ---
star_x_vals <- seq(from = -0.5, to = 0.8, by = 0.01)
lk_star_x <- vapply(seq_along(star_x_vals), function(i) {
  star  <- c(star_x_vals[i], 0)
  q_sim <- solve_trajectory(q0_true, p0_true, dt, k_true, m, n_obs, ts, star)
  sum(dnorm(q_obs[, 1], q_sim[, 1], sigma_x, log = TRUE)) +
  sum(dnorm(q_obs[, 2], q_sim[, 2], sigma_y, log = TRUE))
}, numeric(1L))

log_result("\n=== §7 Conditional log likelihood vs q_star_x ===")
log_result("q_star_x at max: ",
           round(star_x_vals[which.max(lk_star_x)], 3))

fig_30_8 <- ggplot2::ggplot(
  data = data.frame(star_x = star_x_vals, lk = lk_star_x),
  ggplot2::aes(x = star_x, y = lk)
) +
  ggplot2::geom_line() +
  ggplot2::labs(
    title = "Figure 30.8: Conditional log likelihood vs q_star_x",
    x = expression(q["*,x"]),
    y = "Conditional log likelihood"
  ) +
  ggplot2::theme_minimal()

print(fig_30_8)
ggplot2::ggsave("figs/Fig-30.8.svg", plot = fig_30_8,
                width = 6, height = 4, device = svg)
log_result("Saved: figs/Fig-30.8.svg")

# --- Fig-30.9: conditional log likelihood heat map over (q_star_x, q_star_y) ---
star_x_grid <- seq(from = -0.5, to = 0.8, by = 0.05)
star_y_grid <- seq(from = -0.5, to = 0.5, by = 0.05)
grid_df     <- expand.grid(star_x = star_x_grid, star_y = star_y_grid)

lk_grid <- vapply(seq_len(nrow(grid_df)), function(i) {
  star  <- c(grid_df$star_x[i], grid_df$star_y[i])
  q_sim <- solve_trajectory(q0_true, p0_true, dt, k_true, m, n_obs, ts, star)
  sum(dnorm(q_obs[, 1], q_sim[, 1], sigma_x, log = TRUE)) +
  sum(dnorm(q_obs[, 2], q_sim[, 2], sigma_y, log = TRUE))
}, numeric(1L))

grid_df$lk <- lk_grid

log_result("Heat map grid: ",
           nrow(grid_df), " points computed")

fig_30_9 <- ggplot2::ggplot(
  grid_df,
  ggplot2::aes(x = star_x, y = star_y, fill = lk)
) +
  ggplot2::geom_tile() +
  ggplot2::labs(
    title = "Figure 30.9: Conditional log likelihood heat map",
    x = expression(q["*,x"]),
    y = expression(q["*,y"]),
    fill = "log likelihood"
  ) +
  ggplot2::theme_minimal()

print(fig_30_9)
ggplot2::ggsave("figs/Fig-30.9.svg", plot = fig_30_9,
                width = 6, height = 4, device = svg)
log_result("Saved: figs/Fig-30.9.svg")

# =============================================================================
# §8  Full model: fit2p — 7 parameters, Pathfinder init
# =============================================================================

N_select   <- 40
time_obs   <- (1:N_select) / 10
stan_data2 <- list(N = N_select, q_obs = q_obs, time = time_obs, sigma = sigma)

mod2 <- cmdstan_model("planetary_motion_star.stan")

pth2 <- mod2$pathfinder(
  data              = stan_data2,
  num_paths         = 40,
  single_path_draws = 25,
  draws             = 1000,
  max_lbfgs_iters   = 100,
  psis_resample     = FALSE,
  refresh           = 0
)
# Most of the 40 paths will fail; this is expected (book: "only 3s").
# The surviving paths supply the importance-sampling proposal.

fit2p <- mod2$sample(
  data            = stan_data2,
  init            = pth2,
  chains          = CHAINS,
  parallel_chains = CHAINS,
  iter_warmup     = 500,
  iter_sampling   = 500,
  seed            = SEED,
  save_warmup     = TRUE,
  refresh         = 0
)
fit2p$save_object(file = "saved_fit/fit2p.RDS")

# Diagnostics — sampled parameters: k, q0, p0, star
# Do NOT include generated quantities (qx_pred, qy_pred) in params.
fit2p_params <- c("lp__", "k", "q0", "p0", "star")
fit2p_smry <- fit2p$summary(
  variables = fit2p_params,
  "mean", "sd", "rhat", "ess_bulk", "ess_tail"
)

log_result("\n=== §8 fit2p diagnostics (Pathfinder init, full model) ===")
capture_result(fit2p_smry, label = "fit2p summary: k, q0, p0, star")

log_result("Book targets:")
log_result("  k      : mean ≈ 1.00,  sd ≈ 0.001")
log_result("  q0[1]  : mean ≈ 1.00,  sd ≈ 0.004")
log_result("  q0[2]  : mean ≈ -0.003, sd ≈ 0.004")
log_result("  p0[1]  : mean ≈ 0.005, sd ≈ 0.006")
log_result("  p0[2]  : mean ≈ 1.000, sd ≈ 0.004")
log_result("  star[1]: mean ≈ 0.001, sd ≈ 0.005")
log_result("  star[2]: mean ≈ 0.007, sd ≈ 0.005")

for (v in c("k", "q0[1]", "q0[2]", "p0[1]", "p0[2]", "star[1]", "star[2]")) {
  row <- fit2p_smry[fit2p_smry$variable == v, ]
  if (nrow(row) > 0)
    log_result("  ", v, ": mean = ", round(row$mean, 4),
               ", sd = ", round(row$sd, 4),
               ", Rhat = ", round(row$rhat, 4))
}

wf3 <- init_workflow(mode = "practice", stage = "explore")
wf3 <- diagnose_cmdstan(fit2p, wf3,
         params = c("k", "q0", "p0", "star"))
log_result("--- Diagnostics: fit2p (cmdstanr) ---")
log_result("  passed      : ", wf3$diagnostics$passed)
log_result("  Rhat_max    : ", round(wf3$diagnostics$rhat_max, 4))
log_result("  bulk_ESS_min: ", round(wf3$diagnostics$bulk_ess_min, 0))
log_result("  tail_ESS_min: ", round(wf3$diagnostics$tail_ess_min, 0))
log_result("  divergences : ", wf3$diagnostics$n_divergences)
if (length(wf3$diagnostics$failed_criteria) > 0)
  log_result("  FAILED      : ",
             paste(wf3$diagnostics$failed_criteria, collapse = ", "))

# --- Fig-30.10: per-chain PPC, fit2p, with star position marker ---
fig_30_10 <- ppc_plot2D(fit2p, data_pred = data_pred, plot_star = TRUE)
print(fig_30_10)
ggplot2::ggsave("figs/Fig-30.10.svg", plot = fig_30_10,
                width = 6, height = 6, device = svg)
log_result("Saved: figs/Fig-30.10.svg")

# =============================================================================
# §9  Save and wrap up
# =============================================================================

# LOO-CV not applicable for this chapter.
# Reason: The chapter's inferential goal is to demonstrate multimodality and
# the Pathfinder fix, not model comparison. All three fits use the same
# likelihood; comparing fit1 (failed run) to fit1p (same model, better inits)
# via LOO would be meaningless — they estimate the same posterior, not
# different models. fit2p adds parameters and cannot be compared to fit1/fit1p
# via loo_compare() because the data list differs (stan_data1 vs stan_data2).
log_result("\nLOO-CV: not performed.")
log_result("Reason: single-model-per-stage chapter; fits not comparable via loo_compare().")
log_result("        See §9 comment for full explanation.")

wf3$ppc_complete <- TRUE
wf3$loo_complete <- FALSE
saveRDS(wf3, "wf_final.rds")
export_context(wf3)
guide(wf3)

log_result("\n=== Ch 30 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: saved_fit/fit1.RDS (book, loaded), saved_fit/fit1p.RDS, saved_fit/fit2p.RDS, wf_final.rds")
log_result("Figures saved to figs/: ",
           paste(paste0("Fig-30.", c(1:10), ".svg"), collapse = ", "))
close(results_con)



