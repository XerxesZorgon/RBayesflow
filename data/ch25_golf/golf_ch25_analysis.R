# data/ch25_golf/golf_ch25_analysis.R
#
# Ch 25 — "Model building and expansion: Golf putting"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# This script fits a sequence of eleven Stan models to golf putting data,
# demonstrating iterative model expansion: logistic regression → geometry-based
# angle model → angle + distance model → normal approximation → binomial with
# proportional errors (vector and constant) → fully parameterised variants.
# The script follows the book's progression exactly, including Pathfinder
# initialisations for models with convergence difficulties. LOO comparison
# uses integrated PSIS-LOO via standalone generated quantities.
#
# This is a Stan-native chapter; run_phase3() is bypassed throughout.
# All models are fitted directly via cmdstanr model objects compiled from
# the .stan files in this folder.
#
# Deliberate omissions:
#   - Prior predictive checks (the book does not show them for this chapter)
#   - brms-level wf_state diagnostics for models 3–11 (cmdstanr only; wf_state
#     is used for models 1 and 2 only, where the interface is straightforward)
#   - Exact Fig 3, 4, 8 (geometry sketches using base R symbols) are reproduced
#     as informational plots and saved, but not deeply commented
#   - K-fold CV is not run (integrated PSIS-LOO is used throughout)
#
# Data:
#   golf_data.txt     — Berry (1996), 19 bins, 2–20 ft; downloaded from book repo
#   golf_data_new.txt — Broadie (2018), 31 bins, 0.28–75 ft; downloaded from book repo
#   Both cached to data/ subfolder on first run; no re-download if file exists.
#
# Models covered:
#   model_1 / fit_1  — golf_logistic.stan         : logistic regression, y|a,b ~ binomial_logit
#   model_2 / fit_2  — golf_angle_binomial.stan    : geometry-based angle model, one param sigma
#   model_3 / fit_3  — golf_angle_distance_binomial.stan : angle + distance (Broadie), new data
#   model_4 / fit_4  — golf_angle_distance_normal.stan   : normal approx to binomial + fudge sigma_y
#   model_4r/fit_4r  — golf_angle_distance_normal_with_resids.stan : same + residuals GQ
#   model_5 / fit_5  — golf_angle_distance_binomial_with_logit_errors.stan : logit-scale errors
#   model_6 / fit_6  — golf_angle_distance_binomial_with_proportional_errors.stan : prop errors (vec)
#   model_7 / fit_7  — golf_angle_distance_binomial_with_proportional_errors_2.stan : + dist_tol param
#   model_8 / fit_8  — golf_angle_distance_binomial_with_proportional_errors_3.stan : + overshot param
#   model_9 / fit_9  — golf_angle_distance_binomial_with_constant_errors.stan      : scalar epsilon
#   model_10/fit_10  — golf_angle_distance_binomial_with_constant_errors_2.stan    : + dist_tol param
#   model_11/fit_11  — golf_angle_distance_binomial_with_constant_errors_3.stan    : + overshot param
#
# Figures produced:
#   Fig-25.1.svg  — Data: proportion of successful putts vs distance with SE bars (old data)
#   Fig-25.2.svg  — Fitted logistic regression with 10 posterior draws (green) + mean (black)
#   Fig-25.3.svg  — Geometry sketch: ball, hole, threshold angle
#   Fig-25.4.svg  — Normal distribution of shot angle with ±2σ labels
#   Fig-25.5.svg  — Modelled Pr(success) for σ = 0.5°, 1°, 2°, 5°, 20°
#   Fig-25.6.svg  — Two models fit to old data: logistic + geometry-based
#   Fig-25.7.svg  — New data (red) + old data (blue) + old model prediction
#   Fig-25.8.svg  — Geometry sketch: angle + distance gray zone
#   Fig-25.9.svg  — Model 3 fit to new data
#   Fig-25.10.svg — Model 4 (normal approx) fit to new data
#   Fig-25.11.svg — Residuals from model 4
#   Fig-25.12.svg — Model 5 (logit errors) fit to new data
#   Fig-25.13.svg — Model 6 (proportional errors) fit to new data
#   Fig-25.14.svg — Residuals from model 6
#   Fig-25.15.svg — Model 7 (dist_tol param) fit to new data
#   Fig-25.16.svg — Residuals from model 7
#   Fig-25.17.svg — Model 8 (overshot param) fit to new data
#   Fig-25.18.svg — Residuals from model 8
#   Fig-25.19.svg — Bivariate posterior of dist_tol vs overshot (model 8)
#   Fig-25.20.svg — Model 11 (constant error) fit to new data
#   Fig-25.21.svg — Residuals from model 11
#   Fig-25.22.svg — Pointwise elpd difference: model 8 vs model 11
#
# Book-target posterior summaries (used for success criteria):
#   Model 1 (logistic):
#     a: 2.23 ± 0.06
#     b: -0.26 ± 0.01
#   Model 2 (angle-only):
#     sigma_degrees: 1.53 ± 0.02
#   Model 4 (angle+distance, normal):
#     sigma_degrees: ~1.0
#     sigma_distance: ~0.08
#     sigma_y: ~0.003
#   Model 7 (dist_tol param):
#     distance_tolerance: mean ~3.9, median ~4.0, sd 0.22
#   Model 11 (constant error):
#     sigma_angle: ~0.015, sigma_distance: ~0.13
#     distance_tolerance: ~4.4, overshot: ~1.1, epsilon: ~0.00060
#   LOO comparison (models 6/7/8):
#     elpd_loo: model 6 ~-185.6, model 7 ~-173.4, model 8 ~-173.1
#     elpd_diff (model 7 vs 6): ~-12.6 ± 4.5 (model 7 significantly better)
#     elpd_diff (model 8 vs 7): ~-0.3 ± 0.3 (no significant difference)
#   LOO comparison (models 8 vs 11):
#     elpd_diff (model 11 vs 8): ~-32.5 ± 18.0 (constant error clearly worse)
#
# Stan-native chapter notes:
#   run_phase3() is bypassed; cmdstanr is used directly for all models.
#   Models 5–8 include a per-observation latent vector (epsilon[J], J=31);
#   rubric requires iter = 4000, warmup = 2000 for these models.
#   Models 1–4, 9–11 use cmdstanr defaults (iter = 2000, warmup = 1000).
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch25_golf/.

# =============================================================================
# § 0  Environment setup
# =============================================================================

source("../../R/source_all.R")

library(posterior)   # as_draws_rvars()
library(bayesplot)   # mcmc_scatter()
library(loo)         # loo(), loo_compare(), pointwise()

options(brms.backend = "cmdstanr", mc.cores = 4)

SEED <- 42

# ---------------------------------------------------------------------------
# Thin helper for cmdstanr compilation + sampling (Stan-native chapters).
# Uses 4 chains; caller overrides iter/warmup for latent-parameter models.
# ---------------------------------------------------------------------------
cstan <- function(stan_file, data = list(), seed = SEED, chains = 4,
                  iter = 2000, warmup = 1000, init = NULL) {
  model <- cmdstanr::cmdstan_model(stan_file)
  args  <- list(
    data            = data,
    seed            = seed,
    chains          = chains,
    parallel_chains = chains,
    iter_sampling   = iter - warmup,
    iter_warmup     = warmup,
    refresh         = 0
  )
  if (!is.null(init)) args$init <- init
  do.call(model$sample, args)
}

# ---------------------------------------------------------------------------
# Figure helper: draw a base-R plot to the RStudio plot pane, then copy to
# SVG.  Drawing goes to the active device (the RStudio pane) and dev.copy()
# writes the file without opening a separate window.
#
# Usage:
#   save_fig("figs/Fig-25.1.svg", width = 6, height = 4, {
#     plot(...)
#     points(...)
#   })
#
# Never call dev.new(), svg(), or any other device-opening function inside a
# figure block.  That creates external graphics windows instead of updating
# the RStudio plot pane.
#
# For ggplot objects: assign the plot to a variable, print() it, then ggsave().
# ---------------------------------------------------------------------------
save_fig <- function(path, width = 6, height = 4, expr) {
  expr <- substitute(expr)
  eval(expr, parent.frame())          # draw to the active RStudio pane
  dev.copy(svg, path, width = width, height = height)
  dev.off()                           # close the copy device only
  invisible(path)
}

# Initialise wf_state for audit trail
wf <- init_workflow(mode = "practice", stage = "explore")

dir.create("figs", showWarnings = FALSE)
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

log_result("=== Ch 25 results log — ", format(Sys.time()), " ===")

# Utility functions (from the book)
logit    <- qlogis
invlogit <- plogis
fround   <- function(x, digits) format(round(x, digits), nsmall = digits)

# =============================================================================
# § 1  Data acquisition — download from book repo if not already cached
# =============================================================================

log_result("\n=== § 1  Data acquisition ===")

BASE_URL <- "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/master/golf/data"

fetch_data <- function(filename) {
  local_path <- file.path("data", filename)
  if (!file.exists(local_path)) {
    message("Downloading ", filename, " ...")
    download.file(
      url      = paste0(BASE_URL, "/", filename),
      destfile = local_path,
      quiet    = TRUE
    )
    message("  -> saved to ", local_path)
  } else {
    message(filename, " already cached at ", local_path)
  }
  invisible(local_path)
}

dir.create("data", showWarnings = FALSE)
fetch_data("golf_data.txt")
fetch_data("golf_data_new.txt")

# Old data (Berry 1996) — 19 observations, distances 2–20 feet
golf <- read.table("data/golf_data.txt", header = TRUE, skip = 2)
x    <- golf$x
y    <- golf$y
n    <- golf$n
J    <- length(y)

# New data (Broadie 2018) — 31 observations, distances 0.28–75 feet
golf_new <- read.table("data/golf_data_new.txt", header = TRUE, skip = 2)

# Physical constants (inches → feet)
r <- (1.68 / 2) / 12   # ball radius
R <- (4.25 / 2) / 12   # hole radius

log_result("Old data: J = ", J, " bins, distances ", min(x), "–", max(x), " ft")
log_result("New data: J = ", nrow(golf_new), " bins, distances ",
           min(golf_new$x), "–", max(golf_new$x), " ft")
log_result("r (ball) = ", round(r, 4), " ft  |  R (hole) = ", round(R, 4), " ft")

# Standard errors for plotting
se <- sqrt((y / n) * (1 - y / n) / n)

# =============================================================================
# § 2  Fig 25.1 — Raw data
# =============================================================================

log_result("\n=== § 2  Fig 25.1: Raw data ===")

save_fig("figs/Fig-25.1.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Data on putts in pro golf", type = "n")
  points(x, y / n, pch = 20, col = "blue")
  segments(x, y / n + se, x, y / n - se, lwd = .5, col = "blue")
  text(x + .4, y / n + se + .02,
       paste(y, "/", n, sep = ""), cex = .6, col = "gray40")
})

# =============================================================================
# § 3  Model 1: Logistic regression (old data)
# =============================================================================

log_result("\n=== § 3  Model 1: Logistic regression ===")

golf_data <- list(x = x, y = y, n = n, J = J)

model_1 <- cmdstanr::cmdstan_model("golf_logistic.stan")
fit_1   <- cstan("golf_logistic.stan", data = golf_data)

saveRDS(fit_1, "fit_1.rds")

draws_1 <- fit_1$draws(format = "df")
a_sim   <- draws_1$a
b_sim   <- draws_1$b
a_hat   <- mean(a_sim)
b_hat   <- mean(b_sim)
n_sims  <- nrow(draws_1)

log_result("\n=== M1: Logistic regression ===")
log_result("Book targets: a = 2.23 ± 0.06 | b = -0.26 ± 0.01")
log_result("Fitted:       a = ", round(a_hat, 2), " ± ", round(sd(a_sim), 2),
           " | b = ", round(b_hat, 2), " ± ", round(sd(b_sim), 2))

diag_1 <- fit_1$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_1 ---")
log_result("  Rhat_max     : ", round(max(fit_1$summary()$rhat,      na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_1$summary()$ess_bulk,  na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_1$summary()$ess_tail,  na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_1$num_divergent))

# Fig 25.2: Fitted logistic regression
save_fig("figs/Fig-25.2.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Fitted logistic regression", type = "n")
  for (i in sample(n_sims, 10)) {
    curve(invlogit(a_sim[i] + b_sim[i] * x),
          from = 0, to = 1.1 * max(x), lwd = 0.5, add = TRUE, col = "green")
  }
  curve(invlogit(a_hat + b_hat * x), from = 0, to = 1.1 * max(x), add = TRUE)
  points(x, y / n, pch = 20, col = "blue")
  segments(x, y / n + se, x, y / n - se, lwd = .5, col = "blue")
  text(11, .57, paste("Logistic regression,\n a = ", fround(a_hat, 2),
                      ", b = ", fround(b_hat, 2), sep = ""))
})

# =============================================================================
# § 4  Geometry sketches (Figs 25.3, 25.4, 25.5)
# =============================================================================

log_result("\n=== § 4  Geometry sketches ===")

# Fig 25.3: Geometry of ball, hole, and threshold angle
save_fig("figs/Fig-25.3.svg", width = 7, height = 2, {
  par(mar = c(0, 0, 0, 0))
  dist   <- 2
  r_plot <- r
  R_plot <- R
  plot(0, 0, xlim = c(-R_plot, dist + 3 * R_plot),
       ylim = c(-2 * R_plot, 2 * R_plot),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "", type = "n", asp = 1)
  symbols(0, 0, circles = r_plot, inches = FALSE, add = TRUE)
  symbols(dist, 0, circles = R_plot - r_plot, inches = FALSE, lty = 2, add = TRUE)
  symbols(dist, 0, circles = R_plot, inches = FALSE, add = TRUE)
  curve(0 * x, from = 0, to = dist, add = TRUE)
  curve(((R_plot - r_plot) / dist) * x, from = 0, to = dist, lty = 2, add = TRUE)
  curve(-((R_plot - r_plot) / dist) * x, from = 0, to = dist, lty = 2, add = TRUE)
  text(0.5 * dist, -1.5 * R_plot, "x")
  arrows(0.5 * dist + 0.05, -1.5 * R_plot, dist, -1.5 * R_plot, 2, length = .1)
  arrows(0.5 * dist - 0.05, -1.5 * R_plot, 0,    -1.5 * R_plot, 2, length = .1)
  text(dist + 1.2 * R_plot, .5 * R_plot, "R")
  arrows(dist + 1.2 * R_plot, .7 * R_plot, dist + 1.2 * R_plot, R_plot, length = .05)
  arrows(dist + 1.2 * R_plot, .3 * R_plot, dist + 1.2 * R_plot, 0,      length = .05)
  text(0, r_plot / 2, "r")
})

# Fig 25.4: Normal distribution of shot angle
save_fig("figs/Fig-25.4.svg", width = 7, height = 3, {
  par(mar = c(3, 3, 0, 0), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(-4, 4), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
       xlab = "Angle of shot", ylab = "", type = "n")
  axis(1, seq(-4, 4),
       c("", "", expression(-2 * sigma), "", 0, "", expression(2 * sigma), "", ""))
  curve(dnorm(x) / dnorm(0), add = TRUE)
})

# Fig 25.5: Model-implied success curves for different sigma values
save_fig("figs/Fig-25.5.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = expression(paste("Modeled Pr(success) for different values of ", sigma)),
       type = "n")
  sigma_degrees_plot <- c(0.5, 1, 2, 5, 20)
  x_text <- c(15, 10, 6, 4, 2)
  for (i in seq_along(sigma_degrees_plot)) {
    sigma_i  <- (pi / 180) * sigma_degrees_plot[i]
    x_grid_i <- seq(R - r, 1.1 * max(x), .01)
    p_grid_i <- 2 * pnorm(asin((R - r) / x_grid_i) / sigma_i) - 1
    lines(c(0, R - r, x_grid_i), c(1, 1, p_grid_i))
    text(x_text[i] + 0.7,
         2 * pnorm(asin((R - r) / x_text[i]) / sigma_i) - 1,
         bquote(sigma == .(sigma_degrees_plot[i]) * degree), adj = 0)
  }
})

# =============================================================================
# § 5  Model 2: Geometry-based angle model (old data)
# =============================================================================

log_result("\n=== § 5  Model 2: Geometry-based angle model ===")

golf_data <- c(golf_data, r = r, R = R)

model_2 <- cmdstanr::cmdstan_model("golf_angle_binomial.stan")
fit_2   <- cstan("golf_angle_binomial.stan", data = golf_data)

saveRDS(fit_2, "fit_2.rds")

draws_2           <- fit_2$draws(format = "df")
sigma_sim         <- draws_2$sigma
sigma_degrees_sim <- draws_2$sigma_degrees
sigma_hat         <- mean(sigma_sim)

log_result("\n=== M2: Angle-only geometry model ===")
log_result("Book target: sigma_degrees = 1.53 ± 0.02")
log_result("Fitted:      sigma_degrees = ",
           round(mean(sigma_degrees_sim), 2), " ± ",
           round(sd(sigma_degrees_sim), 2))

diag_2 <- fit_2$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_2 ---")
log_result("  Rhat_max     : ", round(max(fit_2$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_2$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_2$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_2$num_divergent))

# Fig 25.6: Two models on old data
save_fig("figs/Fig-25.6.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Two models fit to the golf putting data", type = "n")
  segments(x, y / n + se, x, y / n - se, lwd = .5)
  curve(invlogit(a_hat + b_hat * x), from = 0, to = 1.1 * max(x), add = TRUE)
  x_grid <- seq(R - r, 1.1 * max(x), .01)
  p_grid <- 2 * pnorm(asin((R - r) / x_grid) / sigma_hat) - 1
  lines(c(0, R - r, x_grid), c(1, 1, p_grid), col = "blue")
  points(x, y / n, pch = 20, col = "blue")
  text(10.3, .58, "Logistic regression")
  text(18.5, .24, "Geometry-based model", col = "blue")
})

# =============================================================================
# § 6  Testing on new data (Fig 25.7)
# =============================================================================

log_result("\n=== § 6  New data vs old model ===")

# Grid over the full new-data range — used by all subsequent model plots
x_grid_new <- seq(R - r, 1.1 * max(golf_new$x), .01)

# Fig 25.7: New data (red) + old data (blue) + old angle model
save_fig("figs/Fig-25.7.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking already-fit model to new data", type = "n")
  p_grid_new <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_hat) - 1
  lines(c(0, R - r, x_grid_new), c(1, 1, p_grid_new), col = "blue")
  points(golf$x,     golf$y / golf$n,         pch = 20, col = "blue")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
  legend(60, 0.4, legend = c("Old data", "New data"),
         col = c("blue", "red"), pch = 20)
})

# Set up new data list (used by models 3–11)
overshot           <- 1
distance_tolerance <- 3
golf_new_data <- list(
  x                  = golf_new$x,
  y                  = golf_new$y,
  n                  = golf_new$n,
  J                  = nrow(golf_new),
  r                  = r,
  R                  = R,
  overshot           = overshot,
  distance_tolerance = distance_tolerance
)

# =============================================================================
# § 7  Fig 25.8 — Geometry sketch: angle + distance gray zone
# =============================================================================

log_result("\n=== § 7  Fig 25.8: angle + distance sketch ===")

save_fig("figs/Fig-25.8.svg", width = 7, height = 2, {
  par(mar = c(0, 0, 0, 0))
  dist                    <- 2
  r_plot                  <- r
  R_plot                  <- R
  distance_tolerance_plot <- 0.6
  plot(0, 0,
       xlim = c(-R_plot, dist + 3 * R_plot + 1.5 * distance_tolerance_plot),
       ylim = c(-2 * R_plot, 2 * R_plot),
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "", type = "n", asp = 1)
  polygon(
    c(dist, dist,
      dist + distance_tolerance_plot, dist + distance_tolerance_plot),
    c(R_plot - r_plot, -(R_plot - r_plot),
      -(R_plot - r_plot) * (dist + distance_tolerance_plot) / dist,
       (R_plot - r_plot) * (dist + distance_tolerance_plot) / dist),
    border = NA, col = "gray"
  )
  symbols(0, 0, circles = r_plot, inches = FALSE, add = TRUE)
  symbols(dist, 0, circles = R_plot, inches = FALSE, add = TRUE)
  symbols(dist, 0, circles = R_plot - r_plot, inches = FALSE,
          lty = 2, bg = "gray", add = TRUE)
  curve(((R_plot - r_plot) / dist) * x,
        from = 0, to = dist + 1.5 * distance_tolerance_plot, lty = 2, add = TRUE)
  curve(-((R_plot - r_plot) / dist) * x,
        from = 0, to = dist + 1.5 * distance_tolerance_plot, lty = 2, add = TRUE)
  text(0.5 * dist, -1.5 * R_plot, "x")
  arrows(0.5 * dist + 0.05, -1.5 * R_plot, dist, -1.5 * R_plot, 2, length = .1)
  arrows(0.5 * dist - 0.05, -1.5 * R_plot, 0,    -1.5 * R_plot, 2, length = .1)
})

# =============================================================================
# § 8  Model 3: Angle + distance, binomial (new data)
# =============================================================================

log_result("\n=== § 8  Model 3: Angle + distance binomial ===")

model_3 <- cmdstanr::cmdstan_model("golf_angle_distance_binomial.stan")

# Initial attempt — may have convergence problems (book notes high Rhat, low ESS)
fit_3 <- cstan("golf_angle_distance_binomial.stan", data = golf_new_data)
# <check fit_3 summary before proceeding — book reports multimodality here>

# Re-initialise with Pathfinder
pth_3 <- model_3$pathfinder(
  data = golf_new_data, refresh = 0,
  num_paths = 20, max_lbfgs_iters = 100
)
fit_3 <- cstan("golf_angle_distance_binomial.stan", data = golf_new_data,
               init = pth_3)

saveRDS(fit_3, "fit_3.rds")

draws_3            <- fit_3$draws(format = "df")
sigma_angle_hat    <- mean(draws_3$sigma_angle)
sigma_distance_hat <- mean(draws_3$sigma_distance)

log_result("\n=== M3: Angle + distance binomial (new data) ===")
log_result("Book target: sigma_degrees ~0.76")
log_result("Fitted:      sigma_degrees = ", round(mean(draws_3$sigma_degrees), 2))

diag_3 <- fit_3$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_3 (after Pathfinder init) ---")
log_result("  Rhat_max     : ", round(max(fit_3$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_3$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_3$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_3$num_divergent))

# Fig 25.9: Model 3 fit
save_fig("figs/Fig-25.9.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

# =============================================================================
# § 9  Model 4: Normal approximation (new data)
# =============================================================================

log_result("\n=== § 9  Model 4: Normal approximation ===")

model_4 <- cmdstanr::cmdstan_model("golf_angle_distance_normal.stan")
fit_4   <- cstan("golf_angle_distance_normal.stan", data = golf_new_data)

saveRDS(fit_4, "fit_4.rds")

draws_4            <- fit_4$draws(format = "df")
sigma_angle_hat    <- mean(draws_4$sigma_angle)
sigma_distance_hat <- mean(draws_4$sigma_distance)

log_result("\n=== M4: Normal approximation (new data) ===")
log_result("Book targets: sigma_degrees ~1.0, sigma_distance ~0.08, sigma_y ~0.003")
log_result("Fitted: sigma_degrees = ", round(mean(draws_4$sigma_degrees), 2),
           " | sigma_distance = ",     round(mean(draws_4$sigma_distance), 3),
           " | sigma_y = ",            round(mean(draws_4$sigma_y), 4))

diag_4 <- fit_4$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_4 ---")
log_result("  Rhat_max     : ", round(max(fit_4$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_4$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_4$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_4$num_divergent))

# Fig 25.10: Model 4 fit
save_fig("figs/Fig-25.10.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

# Model 4 with residuals
model_4r <- cmdstanr::cmdstan_model("golf_angle_distance_normal_with_resids.stan")
fit_4r   <- cstan("golf_angle_distance_normal_with_resids.stan", data = golf_new_data)

saveRDS(fit_4r, "fit_4r.rds")

posterior_mean_residual_4 <- mean(
  posterior::as_draws_rvars(fit_4r$draws())$residual
)

# Fig 25.11: Residuals from model 4
save_fig("figs/Fig-25.11.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, posterior_mean_residual_4,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "y/n - fitted E(y/n)",
       main = "Residuals from fitted model", type = "n")
  abline(0, 0, col = "gray", lty = 2)
  lines(golf_new$x, posterior_mean_residual_4)
})

# =============================================================================
# § 10  Model 5: Binomial with logit-scale errors (new data)
# =============================================================================

log_result("\n=== § 10  Model 5: Binomial + logit errors ===")

model_5 <- cmdstanr::cmdstan_model("golf_angle_distance_binomial_with_logit_errors.stan")

# eta[J] is a per-observation latent vector -> iter = 4000, warmup = 2000
# First attempt (may produce warnings / poor mixing); Pathfinder fixes init
fit_5 <- cstan("golf_angle_distance_binomial_with_logit_errors.stan",
               data = golf_new_data, iter = 4000, warmup = 2000)

pth_5 <- model_5$pathfinder(
  data = golf_new_data, refresh = 0,
  num_paths = 20, max_lbfgs_iters = 100
)
fit_5 <- cstan("golf_angle_distance_binomial_with_logit_errors.stan",
               data = golf_new_data, iter = 4000, warmup = 2000,
               init = pth_5)

saveRDS(fit_5, "fit_5.rds")

draws_5            <- fit_5$draws(format = "df")
sigma_angle_hat    <- mean(draws_5$sigma_angle)
sigma_distance_hat <- mean(draws_5$sigma_distance)

log_result("Fitted sigma_angle = ",    round(mean(draws_5$sigma_angle),    4),
           " | sigma_distance = ",     round(mean(draws_5$sigma_distance), 3),
           " | sigma_eta = ",          round(mean(draws_5$sigma_eta),      3))

diag_5 <- fit_5$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_5 (after Pathfinder init) ---")
log_result("  Rhat_max     : ", round(max(fit_5$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_5$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_5$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_5$num_divergent))

# Fig 25.12: Model 5 fit
save_fig("figs/Fig-25.12.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

# =============================================================================
# § 11  Model 6: Binomial + proportional errors vector (new data)
# =============================================================================

log_result("\n=== § 11  Model 6: Proportional errors (vector epsilon) ===")

model_6 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_proportional_errors.stan"
)
# epsilon[J] per-observation -> iter = 4000, warmup = 2000
fit_6 <- cstan("golf_angle_distance_binomial_with_proportional_errors.stan",
               data = golf_new_data, iter = 4000, warmup = 2000)

saveRDS(fit_6, "fit_6.rds")

draws_6            <- fit_6$draws(format = "df")
sigma_angle_hat    <- mean(draws_6$sigma_angle)
sigma_distance_hat <- mean(draws_6$sigma_distance)

log_result("Fitted sigma_angle = ",    round(mean(draws_6$sigma_angle),    4),
           " | sigma_distance = ",     round(mean(draws_6$sigma_distance), 3),
           " | sigma_epsilon = ",      round(mean(draws_6$sigma_epsilon),  4))

diag_6 <- fit_6$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_6 ---")
log_result("  Rhat_max     : ", round(max(fit_6$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_6$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_6$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_6$num_divergent))

# Fig 25.13: Model 6 fit
save_fig("figs/Fig-25.13.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

# Residuals from model 6 (generated quantities block)
posterior_mean_residual_6 <- mean(
  posterior::as_draws_rvars(draws_6)$residual
)

# Fig 25.14: Residuals from model 6
save_fig("figs/Fig-25.14.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, posterior_mean_residual_6,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "y/n - fitted E(y/n)",
       main = "Residuals from fitted model", type = "n")
  abline(0, 0, col = "gray", lty = 2)
  lines(golf_new$x, posterior_mean_residual_6)
})

# =============================================================================
# § 12  Model 7: + distance_tolerance as parameter (new data)
# =============================================================================

log_result("\n=== § 12  Model 7: distance_tolerance as parameter ===")

golf_new_data_m7 <- golf_new_data
golf_new_data_m7$distance_tolerance <- NULL   # now a parameter

model_7 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_proportional_errors_2.stan"
)
# epsilon[J] still present -> iter = 4000; init from fit_6
fit_7 <- cstan("golf_angle_distance_binomial_with_proportional_errors_2.stan",
               data = golf_new_data_m7, iter = 4000, warmup = 2000,
               init = fit_6)

saveRDS(fit_7, "fit_7.rds")

draws_7              <- fit_7$draws(format = "df")
sigma_angle_hat      <- mean(draws_7$sigma_angle)
sigma_distance_hat   <- mean(draws_7$sigma_distance)
distance_tolerance_7 <- mean(draws_7$distance_tolerance)

log_result("\n=== M7: distance_tolerance parameter ===")
log_result("Book target: distance_tolerance mean ~3.9, median ~4.0, sd 0.22")
capture_result(fit_7$summary(variables = "distance_tolerance"),
               label = "M7 distance_tolerance posterior:")

diag_7 <- fit_7$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_7 ---")
log_result("  Rhat_max     : ", round(max(fit_7$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_7$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_7$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_7$num_divergent))

# Fig 25.15: Model 7 fit — use posterior mean of distance_tolerance
distance_tolerance <- distance_tolerance_7
overshot           <- 1

save_fig("figs/Fig-25.15.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

posterior_mean_residual_7 <- mean(
  posterior::as_draws_rvars(draws_7)$residual
)

# Fig 25.16: Residuals from model 7
save_fig("figs/Fig-25.16.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, posterior_mean_residual_7,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "y/n - fitted E(y/n)",
       main = "Residuals from fitted model", type = "n")
  abline(0, 0, col = "gray", lty = 2)
  lines(golf_new$x, posterior_mean_residual_7)
})

# =============================================================================
# § 13  Model 8: + overshot as parameter (new data)
# =============================================================================

log_result("\n=== § 13  Model 8: overshot as parameter ===")

golf_new_data_m8 <- golf_new_data
golf_new_data_m8$distance_tolerance <- NULL   # both now parameters
golf_new_data_m8$overshot           <- NULL

model_8 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_proportional_errors_3.stan"
)
# epsilon[J] still present -> iter = 4000; Pathfinder from fit_7
pth_8 <- model_8$pathfinder(
  data = golf_new_data_m8, refresh = 0,
  num_paths = 40, max_lbfgs_iters = 100,
  init = fit_7
)
fit_8 <- cstan("golf_angle_distance_binomial_with_proportional_errors_3.stan",
               data = golf_new_data_m8, iter = 4000, warmup = 2000,
               init = pth_8)

saveRDS(fit_8, "fit_8.rds")

draws_8            <- fit_8$draws(format = "df")
sigma_angle_hat    <- mean(draws_8$sigma_angle)
sigma_distance_hat <- mean(draws_8$sigma_distance)
distance_tolerance <- mean(draws_8$distance_tolerance)
overshot           <- mean(draws_8$overshot)

log_result("\n=== M8: distance_tolerance + overshot parameters ===")
log_result("Book targets: distance_tolerance mean ~3.5, sd ~0.43 | overshot mean ~0.87, sd ~0.13")
capture_result(fit_8$summary(variables = c("distance_tolerance", "overshot")),
               label = "M8 distance_tolerance + overshot posterior:")

diag_8 <- fit_8$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_8 ---")
log_result("  Rhat_max     : ", round(max(fit_8$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_8$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_8$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_8$num_divergent))

# Fig 25.17: Model 8 fit
save_fig("figs/Fig-25.17.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat) - 1
  p_distance_grid <- pnorm((distance_tolerance - overshot) /
                             ((x_grid_new + overshot) * sigma_distance_hat)) -
                     pnorm(-overshot /
                             ((x_grid_new + overshot) * sigma_distance_hat))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

posterior_mean_residual_8 <- mean(
  posterior::as_draws_rvars(draws_8)$residual
)

# Fig 25.18: Residuals from model 8
save_fig("figs/Fig-25.18.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, posterior_mean_residual_8,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "y/n - fitted E(y/n)",
       main = "Residuals from fitted model", type = "n")
  abline(0, 0, col = "gray", lty = 2)
  lines(golf_new$x, posterior_mean_residual_8)
})

# Fig 25.19: Bivariate scatter of distance_tolerance vs overshot (model 8)
p19 <- fit_8$draws(variables = c("distance_tolerance", "overshot")) |>
  bayesplot::mcmc_scatter() +
  bayesplot::theme_default(base_family = "sans", base_size = 16)
print(p19)
ggplot2::ggsave("figs/Fig-25.19.svg", plot = p19, width = 6, height = 4)

# =============================================================================
# § 14  LOO comparison: models 6 / 7 / 8 (integrated PSIS-LOO)
# =============================================================================

log_result("\n=== § 14  Integrated PSIS-LOO: models 6 / 7 / 8 ===")

gq_ll <- cmdstanr::cmdstan_model("golf_log_lik.stan")

loo_6 <- gq_ll$generate_quantities(
  fit_6$draws(variables = c("sigma_epsilon", "p_angle", "p_distance")),
  data = golf_new_data
)$draws(variables = "log_lik") |> loo::loo()

loo_7 <- gq_ll$generate_quantities(
  fit_7$draws(variables = c("sigma_epsilon", "p_angle", "p_distance")),
  data = golf_new_data_m7
)$draws(variables = "log_lik") |> loo::loo()

loo_8 <- gq_ll$generate_quantities(
  fit_8$draws(variables = c("sigma_epsilon", "p_angle", "p_distance")),
  data = golf_new_data_m8
)$draws(variables = "log_lik") |> loo::loo()

log_result("\n=== LOO comparison (models 6/7/8) ===")
log_result("Book targets: elpd_loo:  M6 ~-185.6 | M7 ~-173.4 | M8 ~-173.1")
log_result("Book targets: elpd_diff: M7 vs M6 ~-12.6 ± 4.5 | M8 vs M7 ~-0.3 ± 0.3")

capture_result(loo_6, label = "LOO M6:")
capture_result(loo_7, label = "LOO M7:")
capture_result(loo_8, label = "LOO M8:")

loo_tab_678 <- loo::loo_compare(list(
  `Distance tolerance and overshot fixed`           = loo_6,
  `Distance tolerance parameter and overshot fixed` = loo_7,
  `Distance tolerance and overshot parameters`      = loo_8
))
capture_result(loo_tab_678, label = "LOO compare M6/M7/M8:")
saveRDS(loo_tab_678, "loo_tab_678.rds")

# =============================================================================
# § 15  Models 9 / 10 / 11: Constant error term
# =============================================================================

log_result("\n=== § 15  Models 9–11: Constant error term ===")

# Model 9: scalar epsilon, fixed overshot and distance_tolerance
model_9 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_constant_errors.stan"
)
fit_9 <- cstan("golf_angle_distance_binomial_with_constant_errors.stan",
               data = golf_new_data)

saveRDS(fit_9, "fit_9.rds")

log_result("\n=== M9 summary (sigma_angle, sigma_distance, epsilon) ===")
capture_result(
  fit_9$summary(variables = c("sigma_angle", "sigma_distance", "epsilon")),
  label = "M9:"
)

# Model 10: scalar epsilon + distance_tolerance as parameter
model_10 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_constant_errors_2.stan"
)
pth_10 <- model_10$pathfinder(
  data = golf_new_data_m7, refresh = 0,
  num_paths = 40, max_lbfgs_iters = 100,
  init = fit_9
)
fit_10 <- cstan("golf_angle_distance_binomial_with_constant_errors_2.stan",
                data = golf_new_data_m7, init = pth_10)

saveRDS(fit_10, "fit_10.rds")

log_result("\n=== M10 summary (sigma_angle, sigma_distance, epsilon) ===")
capture_result(
  fit_10$summary(variables = c("sigma_angle", "sigma_distance", "epsilon")),
  label = "M10:"
)

# Model 11: scalar epsilon + distance_tolerance + overshot as parameters
model_11 <- cmdstanr::cmdstan_model(
  "golf_angle_distance_binomial_with_constant_errors_3.stan"
)
pth_11 <- model_11$pathfinder(
  data = golf_new_data_m8, refresh = 0,
  num_paths = 40, max_lbfgs_iters = 100,
  init = fit_10
)
fit_11 <- cstan("golf_angle_distance_binomial_with_constant_errors_3.stan",
                data = golf_new_data_m8, init = pth_11)

saveRDS(fit_11, "fit_11.rds")

log_result("\n=== M11 summary ===")
log_result("Book targets: sigma_angle ~0.015, sigma_distance ~0.13,")
log_result("              distance_tolerance ~4.4, overshot ~1.1, epsilon ~0.00060")
capture_result(
  fit_11$summary(variables = c("sigma_angle", "sigma_distance",
                               "distance_tolerance", "overshot", "epsilon")),
  label = "M11:"
)

diag_11 <- fit_11$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_11 ---")
log_result("  Rhat_max     : ", round(max(fit_11$summary()$rhat,     na.rm = TRUE), 4))
log_result("  bulk_ESS_min : ", round(min(fit_11$summary()$ess_bulk, na.rm = TRUE), 0))
log_result("  tail_ESS_min : ", round(min(fit_11$summary()$ess_tail, na.rm = TRUE), 0))
log_result("  divergences  : ", sum(diag_11$num_divergent))

# Fig 25.20: Model 11 fit — book plots the curve using draws_8 parameters
draws_11             <- fit_11$draws(format = "df")
sigma_angle_hat_8    <- mean(draws_8$sigma_angle)
sigma_distance_hat_8 <- mean(draws_8$sigma_distance)
distance_tolerance_8 <- mean(draws_8$distance_tolerance)
overshot_8           <- mean(draws_8$overshot)

save_fig("figs/Fig-25.20.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(0, 0, xlim = c(0, 1.1 * max(golf_new$x)), ylim = c(0, 1.02),
       xaxs = "i", yaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "Probability of success",
       main = "Checking model fit", type = "n")
  p_angle_grid    <- 2 * pnorm(asin((R - r) / x_grid_new) / sigma_angle_hat_8) - 1
  p_distance_grid <- pnorm((distance_tolerance_8 - overshot_8) /
                             ((x_grid_new + overshot_8) * sigma_distance_hat_8)) -
                     pnorm(-overshot_8 /
                             ((x_grid_new + overshot_8) * sigma_distance_hat_8))
  lines(c(0, R - r, x_grid_new), c(1, 1, p_angle_grid * p_distance_grid), col = "red")
  points(golf_new$x, golf_new$y / golf_new$n, pch = 20, col = "red")
})

posterior_mean_residual_11 <- mean(
  posterior::as_draws_rvars(draws_11)$residual
)

# Fig 25.21: Residuals from model 11
save_fig("figs/Fig-25.21.svg", width = 4.5, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, posterior_mean_residual_11,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)", ylab = "y/n - fitted E(y/n)",
       main = "Residuals from fitted model", type = "n")
  abline(0, 0, col = "gray", lty = 2)
  lines(golf_new$x, posterior_mean_residual_11)
})

# =============================================================================
# § 16  LOO comparison: model 8 vs model 11
# =============================================================================

log_result("\n=== § 16  LOO: model 8 vs model 11 ===")

gq_ll_c <- cmdstanr::cmdstan_model("golf_log_lik_constant_errors.stan")

loo_11 <- gq_ll_c$generate_quantities(
  fit_11$draws(variables = "p"),
  data = golf_new_data_m8
)$draws(variables = "log_lik") |> loo::loo()

log_result("Book targets: elpd_diff (model 11 vs model 8) ~-32.5 ± 18.0 (constant worse)")
loo_tab_8_11 <- loo::loo_compare(list(
  `Model 8 with varying error terms`  = loo_8,
  `Model 11 with constant error term` = loo_11
))
capture_result(loo_tab_8_11, label = "LOO compare M8 vs M11:")
saveRDS(loo_tab_8_11, "loo_tab_8_11.rds")

# Pointwise elpd difference
pointwise_elpd_diff <- loo::pointwise(loo_8, "elpd_loo") -
                       loo::pointwise(loo_11, "elpd_loo")

# Fig 25.22: Pointwise elpd difference
save_fig("figs/Fig-25.22.svg", width = 6, height = 4, {
  par(mar = c(3, 3, 2, 1), mgp = c(1.7, .5, 0), tck = -.02)
  plot(golf_new$x, pointwise_elpd_diff,
       xlim = c(0, 1.1 * max(golf_new$x)), xaxs = "i", bty = "l",
       xlab = "Distance from hole (feet)",
       ylab = "pointwise elpd_loo difference",
       main = "Predictive performance difference Model 8 vs Model 11")
  abline(0, 0, col = "gray", lty = 2)
})

# =============================================================================
# § 17  Wrap-up
# =============================================================================

log_result("\n=== Ch 25 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved RDS files: fit_1.rds, fit_2.rds, fit_3.rds, fit_4.rds, fit_4r.rds,")
log_result("                 fit_5.rds, fit_6.rds, fit_7.rds, fit_8.rds,")
log_result("                 fit_9.rds, fit_10.rds, fit_11.rds,")
log_result("                 loo_tab_678.rds, loo_tab_8_11.rds")
log_result("Figures: Fig-25.1.svg through Fig-25.22.svg in figs/")

close(results_con)
