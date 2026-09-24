# data/ch28_student_grades/ch28_student_grades_analysis.R
#
# Ch 28 — "Models for regression coefficients and variable selection: Student grades"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Predicts median math (n=407) and Portuguese (n=657 rows, ~382 complete cases)
# exam grades for Portuguese high school students using 26 predictors.
# Demonstrates prior comparison (uniform, scaled normal, horseshoe, R2D2),
# projection predictive variable selection via projpred, and subsampling LOO.
# Deliberately omits: tinytable display tables, latex2exp axis labels
# (replaced with plain text), khroma palette (replaced with RColorBrewer Set1),
# progressr progress bars.
#
# Data:
#   student-merged-all.csv (data/ch28_student_grades/)
#   Source: ROS-Examples GitHub (Cortez & Silva 2008)
#   682 rows × 32 cols; math outcome n=407, Portuguese outcome n=657 rows
#   Key variables: G1mat–G3mat, G1por–G3por, 26 predictors (school through absences)
#
# Models covered:
#   fitm_u   — brm(Gmat ~ ., flat prior, math)
#   fitm_n1  — brm(Gmat ~ ., normal(0,2.5) prior, math)
#   fitm_n1p — fitm_n1 prior-only
#   fitm_n1pt— fitm_n1 prior-only, truncated sigma
#   fitm_n2  — brm(Gmat ~ ., scaled normal prior, math)
#   fitm_n2p — fitm_n2 prior-only
#   fitm_hs  — brm(Gmat ~ ., horseshoe prior, math)
#   fitm_hsp — fitm_hs prior-only
#   fitm     — brm(Gmat ~ ., R2D2 prior, math, iter=5000 then refit iter=2000)
#   fitmp    — fitm prior-only
#   fitp     — brm(Gpor ~ ., R2D2 prior, Portuguese)
#
# Figures produced:
#   Fig-28.1.svg  — Math and Portuguese score distributions (dot plots)
#   Fig-28.2.svg  — Marginal posteriors flat prior (mcmc_areas)
#   Fig-28.3.svg  — Implied prior on R^2 full range (panel a)
#   Fig-28.4.svg  — Implied prior on R^2 zoomed (panel b)
#   Fig-28.5.svg  — Posterior R^2 zoomed (panel c)
#   Fig-28.6.svg  — LOO-R^2 zoomed (panel d)
#   Fig-28.7.svg  — Marginal posteriors R2D2 prior (mcmc_areas)
#   Fig-28.8.svg  — Bivariate Fedu vs Medu (mcmc_scatter)
#   Fig-28.9.svg  — PPC histogram (pp_check hist, 5 draws)
#   Fig-28.10.svg — LOO-PIT-ECDF
#   Fig-28.11.svg — Math fast variable selection path (vselm_fast)
#   Fig-28.12.svg — Math validated variable selection path (vselm)
#   Fig-28.13.svg — Math projected posterior (mcmc_areas, 4 predictors)
#   Fig-28.14.svg — Math search stability (cv_proportions)
#   Fig-28.15.svg — Portuguese marginal posteriors R2D2 (mcmc_areas)
#   Fig-28.16.svg — Portuguese fast variable selection path (vselp_fast)
#   Fig-28.17.svg — Portuguese validated variable selection path (vselp)
#   Fig-28.18.svg — Portuguese projected posterior (mcmc_areas, 8 predictors)
#   Fig-28.19.svg — Portuguese search stability (cv_proportions)
#
# Book-target posterior summaries:
#   bayes_R2(fitm_u): Estimate ≈ 0.30
#   loo_R2(fitm_u):   Estimate ≈ 0.18
#   LOO elpd_diff R2D2 vs wide normal: ≈ -4.1 (se 2.6)
#   bayes_R2(fitp):   Estimate ≈ 0.33
#   loo_R2(fitp):     Estimate ≈ 0.28
#   suggest_size(vselm): 4
#   suggest_size(vselp): 8 (book text says 7; code output says 8)
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch28_student_grades/.

.libPaths(c("C:/Users/johnx/Documents/WildPeaches/Projects/RBayesflow/renv/library/windows/R-4.6/x86_64-w64-mingw32", .libPaths()))
source("../../R/source_all.R")
library(dplyr)
library(matrixStats)
library(ggdist)
library(patchwork)
library(projpred)
library(doFuture)
library(doRNG)
library(RColorBrewer)

options(brms.backend = "cmdstanr", mc.cores = 4)
SEED <- 42

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
log_result("=== Ch 28 results log — ", format(Sys.time()), " ===")

student <- read.csv("student-merged-all.csv")
predictors <- c("school","sex","age","address","famsize","Pstatus",
                "Medu","Fedu","traveltime","studytime","failures","schoolsup",
                "famsup","paid","activities","nursery","higher","internet",
                "romantic","famrel","freetime","goout","Dalc","Walc","health",
                "absences")
p <- length(predictors)

student <- student %>%
  dplyr::mutate(dplyr::across(dplyr::matches("G[1-3](mat|por)"), ~ifelse(. == 0, NA, .))) %>%
  dplyr::mutate(
    Gmat = matrixStats::rowMedians(as.matrix(dplyr::select(., dplyr::matches("G[123]mat"))), na.rm = TRUE),
    Gpor = matrixStats::rowMedians(as.matrix(dplyr::select(., dplyr::matches("G[123]por"))), na.rm = TRUE)
  )

student_Gmat <- subset(student, is.finite(Gmat), select = c("Gmat", predictors))
student_Gmat <- student_Gmat[is.finite(rowMeans(student_Gmat)), ]
student_Gpor <- subset(student, is.finite(Gpor), select = c("Gpor", predictors))

nmat <- nrow(student_Gmat)
npor <- nrow(student_Gpor)
npor_complete <- sum(complete.cases(student_Gpor))

log_result("nmat = ", nmat)
log_result("npor = ", npor)
log_result("npor_complete (complete cases) = ", npor_complete)
log_result("NOTE: brms will fit Portuguese model on complete cases only (", npor_complete, " rows)")

studentstd_Gmat <- student_Gmat
studentstd_Gmat[, predictors] <- scale(student_Gmat[, predictors])
studentstd_Gpor <- student_Gpor
studentstd_Gpor[, predictors] <- scale(student_Gpor[, predictors])

log_result("sd(Gmat) = ", round(sd(student_Gmat$Gmat), 2), " (book: 3.3)")
log_result("sd(Gpor) = ", round(sd(student_Gpor$Gpor, na.rm = TRUE), 2), " (book: ~2.7)")

fig_28_1 <- (
  ggplot2::ggplot(student_Gmat, ggplot2::aes(x = Gmat)) +
    ggdist::geom_dots() +
    ggplot2::labs(x = "Median math exam score") +
    ggplot2::scale_x_continuous(limits = c(0, 20)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
) + (
  ggplot2::ggplot(student_Gpor, ggplot2::aes(x = Gpor)) +
    ggdist::geom_dots() +
    ggplot2::labs(x = "Median Portuguese exam score") +
    ggplot2::scale_x_continuous(limits = c(0, 20)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
)

print(fig_28_1)
ggplot2::ggsave("figs/Fig-28.1.svg", plot = fig_28_1,
                width = 8, height = 4, device = svg)
log_result("Fig-28.1.svg saved")

# =============================================================================
# Phase 3: Flat-prior fit (fitm_u) — §4 of book
# =============================================================================

# Custom bayes_R2 matching Gelman et al. (2019) — plain function, not S3 override
bayes_R2_paper <- function(fit, summary = TRUE, probs = c(0.025, 0.975)) {
  mupred    <- brms::posterior_epred(fit)
  var_mu    <- apply(mupred, 1, var)
  sigma2    <- as.matrix(fit, variable = "sigma")^2
  R2        <- var_mu / (var_mu + sigma2)
  colnames(R2) <- "R2"
  if (summary) R2 <- brms::posterior_summary(R2, probs = probs)
  R2
}

if (file.exists("fitm_u.rds")) {
  fitm_u <- readRDS("fitm_u.rds")
  log_result("fitm_u loaded from cache")
} else {
  wf <- init_workflow(mode = "practice", stage = "explore")
  result_u <- run_phase3(
    wf      = wf,
    formula = Gmat ~ .,
    family  = gaussian(),
    priors  = brms::prior(normal(0, 1e6), class = b),
    data    = studentstd_Gmat,
    seed    = SEED,
    chains  = 4,
    iter    = 2000,
    warmup  = 1000
  )
  fitm_u <- result_u$fit
  wf     <- result_u$wf
  saveRDS(fitm_u, "fitm_u.rds")
  log_result("fitm_u fitted and saved")
}

log_result("\n=== fitm_u: flat prior R2 ===")
log_result("Book target bayes_R2: Estimate ≈ 0.30")
log_result("Book target loo_R2:   Estimate ≈ 0.18")
capture_result(round(as.data.frame(bayes_R2_paper(fitm_u)), 2),
               label = "bayes_R2_paper(fitm_u)")
capture_result(round(as.data.frame(brms::loo_R2(fitm_u)), 2),
               label = "loo_R2(fitm_u)")
br2 <- bayes_R2_paper(fitm_u)["R2", "Estimate"]
lr2 <- brms::loo_R2(fitm_u)["R2", "Estimate"]
if (abs(br2 - 0.30) > 0.05)
  log_result("NOTE: bayes_R2 ", round(br2, 3), " differs from book 0.30 — possible prior/seed difference")
if (abs(lr2 - 0.18) > 0.05)
  log_result("NOTE: loo_R2 ", round(lr2, 3), " differs from book 0.18 — possible prior/seed difference")

drawsmu <- posterior::as_draws_df(fitm_u,
             variable = paste0("b_", predictors)) |>
           posterior::set_variables(predictors)
# Confirm exactly 26 variables (excludes b_Intercept)
stopifnot(length(posterior::variables(drawsmu)) == 26)

fig_28_2 <- bayesplot::mcmc_areas(drawsmu, prob_outer = 0.98,
                                   area_method = "scaled height") +
            ggplot2::xlim(c(-1.5, 1.5)) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Figure 28.2: Marginal posteriors, flat prior")
fig_28_2 <- fig_28_2 +
            ggplot2::scale_y_discrete(limits = rev(levels(fig_28_2$data$parameter)))

print(fig_28_2)
ggplot2::ggsave("figs/Fig-28.2.svg", plot = fig_28_2,
                width = 8, height = 6, device = svg)
log_result("Fig-28.2.svg saved")

# =============================================================================
# Phase 3: Four proper-prior posterior fits — §6 of book
# =============================================================================

# --- fitm_n1: normal(0, 2.5) ---
if (file.exists("fitm_n1.rds")) {
  fitm_n1 <- readRDS("fitm_n1.rds")
  log_result("fitm_n1 loaded from cache")
} else {
  fitm_n1 <- brms::brm(Gmat ~ ., data = studentstd_Gmat,
    prior   = brms::prior(normal(0, 2.5), class = b),
    warmup  = 1000, iter = 5000, refresh = 0,
    seed    = SEED, save_pars = brms::save_pars(all = TRUE))
  saveRDS(fitm_n1, "fitm_n1.rds")
  log_result("fitm_n1 fitted and saved")
}
log_result("--- Diagnostics: fitm_n1 ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitm_n1),        na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitm_n1),  na.rm = TRUE), 4))

# --- fitm_n2: scaled normal ---
# Book uses sd(studentstd_Gmat$Gmat) which equals 1.0 after standardisation.
# Numerically identical to sqrt(0.3/26)*1, but we match the book's expression.
scale_b <- sqrt(0.3 / 26) * sd(studentstd_Gmat$Gmat)
log_result("scale_b (scaled normal sd) = ", round(scale_b, 4),
           "  (note: sd(studentstd_Gmat$Gmat) = 1.0 after standardisation)")
if (file.exists("fitm_n2.rds")) {
  fitm_n2 <- readRDS("fitm_n2.rds")
  log_result("fitm_n2 loaded from cache")
} else {
  fitm_n2 <- brms::brm(Gmat ~ ., data = studentstd_Gmat,
    prior    = brms::prior(normal(0, scale_b), class = b),
    warmup   = 1000, iter = 5000, refresh = 0,
    seed     = SEED, save_pars = brms::save_pars(all = TRUE),
    stanvars = brms::stanvar(scale_b, name = "scale_b"))
  saveRDS(fitm_n2, "fitm_n2.rds")
  log_result("fitm_n2 fitted and saved")
}
log_result("--- Diagnostics: fitm_n2 ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitm_n2),        na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitm_n2),  na.rm = TRUE), 4))

# --- fitm_hs: regularized horseshoe ---
p0 <- 6
scale_slab   <- sd(studentstd_Gmat$Gmat) / sqrt(p0) * sqrt(0.3)
scale_global <- p0 / (p - p0) / sqrt(nrow(studentstd_Gmat))
log_result("horseshoe scale_global = ", round(scale_global, 4),
           "  scale_slab = ", round(scale_slab, 4))
if (file.exists("fitm_hs.rds")) {
  fitm_hs <- readRDS("fitm_hs.rds")
  log_result("fitm_hs loaded from cache")
} else {
  fitm_hs <- brms::brm(Gmat ~ ., data = studentstd_Gmat,
    prior  = brms::prior(horseshoe(scale_global = scale_global,
                                   scale_slab   = scale_slab), class = b),
    warmup = 1000, iter = 5000, refresh = 0,
    seed   = SEED, save_pars = brms::save_pars(all = TRUE))
  saveRDS(fitm_hs, "fitm_hs.rds")
  log_result("fitm_hs fitted and saved")
}
log_result("--- Diagnostics: fitm_hs ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitm_hs),        na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitm_hs),  na.rm = TRUE), 4))

# --- fitm: R2D2 (iter=5000, then refit at iter=2000 in Task 008) ---
if (file.exists("fitm.rds")) {
  fitm <- readRDS("fitm.rds")
  log_result("fitm loaded from cache")
} else {
  fitm <- brms::brm(Gmat ~ ., data = studentstd_Gmat,
    prior  = brms::prior(R2D2(mean_R2 = 1/3, prec_R2 = 3,
                              cons_D2 = 1/2), class = b),
    warmup = 1000, iter = 5000, refresh = 0,
    seed   = SEED, save_pars = brms::save_pars(all = TRUE))
  saveRDS(fitm, "fitm.rds")
  log_result("fitm (R2D2, iter=5000) fitted and saved")
}
log_result("--- Diagnostics: fitm (R2D2) ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitm),           na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitm),     na.rm = TRUE), 4))

# =============================================================================
# Prior-only fits for R^2 prior comparison — §6 of book
# =============================================================================

# fitm_n1p: normal(0,2.5) prior only
if (file.exists("fitm_n1p.rds")) {
  fitm_n1p <- readRDS("fitm_n1p.rds")
  log_result("fitm_n1p loaded from cache")
} else {
  fitm_n1p <- stats::update(fitm_n1, sample_prior = "only", refresh = 0)
  saveRDS(fitm_n1p, "fitm_n1p.rds")
  log_result("fitm_n1p fitted and saved")
}
log_result("Diagnostics skipped: prior-only fit (fitm_n1p)")

# fitm_n1pt: normal(0,2.5) prior only, truncated sigma for zoomed plot
if (file.exists("fitm_n1pt.rds")) {
  fitm_n1pt <- readRDS("fitm_n1pt.rds")
  log_result("fitm_n1pt loaded from cache")
} else {
  fitm_n1pt <- stats::update(fitm_n1, sample_prior = "only", refresh = 0,
    prior = c(brms::prior(normal(0, 2.5), class = b),
              brms::prior(student_t(3, 0, 3), lb = 5, class = sigma)))
  saveRDS(fitm_n1pt, "fitm_n1pt.rds")
  log_result("fitm_n1pt fitted and saved")
}
log_result("Diagnostics skipped: prior-only fit (fitm_n1pt)")

# fitm_n2p: scaled normal prior only
if (file.exists("fitm_n2p.rds")) {
  fitm_n2p <- readRDS("fitm_n2p.rds")
  log_result("fitm_n2p loaded from cache")
} else {
  fitm_n2p <- stats::update(fitm_n2, sample_prior = "only", refresh = 0)
  saveRDS(fitm_n2p, "fitm_n2p.rds")
  log_result("fitm_n2p fitted and saved")
}
log_result("Diagnostics skipped: prior-only fit (fitm_n2p)")

# fitm_hsp: horseshoe prior only
if (file.exists("fitm_hsp.rds")) {
  fitm_hsp <- readRDS("fitm_hsp.rds")
  log_result("fitm_hsp loaded from cache")
} else {
  fitm_hsp <- stats::update(fitm_hs, sample_prior = "only", refresh = 0)
  saveRDS(fitm_hsp, "fitm_hsp.rds")
  log_result("fitm_hsp fitted and saved")
}
log_result("Diagnostics skipped: prior-only fit (fitm_hsp)")

# fitmp: R2D2 prior only
if (file.exists("fitmp.rds")) {
  fitmp <- readRDS("fitmp.rds")
  log_result("fitmp loaded from cache")
} else {
  fitmp <- stats::update(fitm, sample_prior = "only", refresh = 0)
  saveRDS(fitmp, "fitmp.rds")
  log_result("fitmp fitted and saved")
}
log_result("Diagnostics skipped: prior-only fit (fitmp)")

# =============================================================================
# R^2 prior comparison plots (Figs 28.3–28.6) and LOO comparison — §6 of book
# =============================================================================

# Build R2 draws data frame for all 8 fits (4 prior-only + 4 posterior)
types      <- factor(c("Prior", "Posterior"), levels = c("Prior", "Posterior"))
types      <- rep(types, times = 4)
priornames <- factor(c("Wide normal", "Scaled normal", "RHS", "R2D2"),
                     levels = c("Wide normal", "Scaled normal", "RHS", "R2D2"))
priornames <- rep(priornames, each = 2)

# Panel a: prior on R^2, full range (0,1) — uses fitm_n1p, fitm_n1, fitm_n2p,
#           fitm_n2, fitm_hsp, fitm_hs, fitmp, fitm
fits_full <- list(fitm_n1p, fitm_n1, fitm_n2p, fitm_n2,
                  fitm_hsp, fitm_hs, fitmp, fitm)

clr <- RColorBrewer::brewer.pal(4, "Set1")
names(clr) <- c("Wide normal", "Scaled normal", "RHS", "R2D2")

r2_df_full <- do.call(rbind, lapply(seq_along(fits_full), function(i) {
  data.frame(
    R2        = as.numeric(bayes_R2_paper(fits_full[[i]], summary = FALSE)),
    type      = types[i],
    priorname = priornames[i]
  )
}))

fig_28_3 <- r2_df_full |>
  dplyr::filter(type == "Prior") |>
  ggplot2::ggplot(ggplot2::aes(x = R2, color = priorname)) +
  ggdist::stat_slab(density = "bounded", expand = TRUE, trim = FALSE,
                    alpha = 0.6, fill = NA, adjust = 2) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::xlim(c(0, 1)) +
  ggplot2::scale_color_manual(values = clr) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "none") +
  ggplot2::labs(x = "R2", title = "a) Implied prior on R2") +
  ggplot2::annotate("text", x = 0.94, y = 0.33, label = "Wide normal",
                    hjust = 1, color = clr["Wide normal"], size = 4) +
  ggplot2::annotate("text", x = 0.85, y = 0.12, label = "Scaled normal",
                    hjust = 1, color = clr["Scaled normal"], size = 4) +
  ggplot2::annotate("text", x = 0.06, y = 0.33, label = "Horseshoe",
                    hjust = 0, color = clr["RHS"], size = 4) +
  ggplot2::annotate("text", x = 0.15, y = 0.12, label = "R2D2",
                    hjust = 0, color = clr["R2D2"], size = 4)

print(fig_28_3)
ggplot2::ggsave("figs/Fig-28.3.svg", plot = fig_28_3,
                width = 6, height = 4, device = svg)
log_result("Fig-28.3.svg saved")

# Panel b: prior on R^2, zoomed (0.041, 0.419) — uses fitm_n1pt instead of fitm_n1p
fits_zoom <- list(fitm_n1pt, fitm_n1, fitm_n2p, fitm_n2,
                  fitm_hsp, fitm_hs, fitmp, fitm)

r2_df_zoom <- do.call(rbind, lapply(seq_along(fits_zoom), function(i) {
  data.frame(
    R2        = as.numeric(bayes_R2_paper(fits_zoom[[i]], summary = FALSE)),
    type      = types[i],
    priorname = priornames[i]
  )
}))

fig_28_4 <- r2_df_zoom |>
  dplyr::filter(type == "Prior") |>
  ggplot2::ggplot(ggplot2::aes(x = R2, color = priorname)) +
  ggdist::stat_slab(density = "bounded", expand = TRUE, trim = FALSE,
                    alpha = 0.5, fill = NA, adjust = 2) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::xlim(c(0.041, 0.419)) +
  ggplot2::scale_color_manual(values = clr) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "none") +
  ggplot2::labs(x = "R2", title = "b) Implied prior on R2 (zoomed)") +
  ggplot2::annotate("text", x = 0.415, y = 0.61, label = "Wide normal",
                    hjust = 1, color = clr["Wide normal"], size = 4) +
  ggplot2::annotate("text", x = 0.415, y = 0.34, label = "Scaled normal",
                    hjust = 1, color = clr["Scaled normal"], size = 4) +
  ggplot2::annotate("text", x = 0.06,  y = 0.87, label = "Horseshoe",
                    hjust = 0, color = clr["RHS"], size = 4) +
  ggplot2::annotate("text", x = 0.06,  y = 0.28, label = "R2D2",
                    hjust = 0, color = clr["R2D2"], size = 4)

print(fig_28_4)
ggplot2::ggsave("figs/Fig-28.4.svg", plot = fig_28_4,
                width = 6, height = 4, device = svg)
log_result("Fig-28.4.svg saved")

# Panel c: posterior R^2 zoomed
fig_28_5 <- r2_df_zoom |>
  dplyr::filter(type == "Posterior") |>
  ggplot2::ggplot(ggplot2::aes(x = R2, color = priorname)) +
  ggdist::stat_slab(density = "unbounded", expand = TRUE, trim = FALSE,
                    alpha = 0.5, fill = NA, adjust = 2) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::xlim(c(0.04, 0.42)) +
  ggplot2::scale_color_manual(values = clr) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "none") +
  ggplot2::labs(x = "R2", title = "c) Posterior-R2") +
  ggplot2::annotate("text", x = 0.28, y = 0.95, label = "Wide normal",
                    hjust = 0, color = clr["Wide normal"], size = 4) +
  ggplot2::annotate("text", x = 0.26, y = 0.95, label = "Scaled normal",
                    hjust = 1, color = clr["Scaled normal"], size = 4) +
  ggplot2::annotate("text", x = 0.19, y = 0.73, label = "Horseshoe",
                    hjust = 1, color = clr["RHS"], size = 4) +
  ggplot2::annotate("text", x = 0.21, y = 0.84, label = "R2D2",
                    hjust = 1, color = clr["R2D2"], size = 4)

print(fig_28_5)
ggplot2::ggsave("figs/Fig-28.5.svg", plot = fig_28_5,
                width = 6, height = 4, device = svg)
log_result("Fig-28.5.svg saved")

# Panel d: LOO-R^2 zoomed
loo_r2_df <- do.call(rbind, lapply(seq_along(fits_zoom), function(i) {
  data.frame(
    R2        = as.numeric(brms::loo_R2(fits_zoom[[i]], summary = FALSE)),
    type      = types[i],
    priorname = priornames[i]
  )
}))

fig_28_6 <- loo_r2_df |>
  dplyr::filter(type == "Posterior") |>
  ggplot2::ggplot(ggplot2::aes(x = R2, color = priorname)) +
  ggdist::stat_slab(density = "unbounded", expand = TRUE, trim = FALSE,
                    alpha = 0.5, fill = NA, adjust = 2) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::xlim(c(0.04, 0.42)) +
  ggplot2::scale_color_manual(values = clr) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank(),
                 axis.line.y = ggplot2::element_blank(),
                 axis.title.y = ggplot2::element_blank(),
                 legend.position = "none") +
  ggplot2::labs(x = "R2", title = "d) LOO-R2") +
  ggplot2::annotate("text", x = 0.175, y = 0.72, label = "Wide normal",
                    hjust = 1, color = clr["Wide normal"], size = 4) +
  ggplot2::annotate("text", x = 0.235, y = 0.72, label = "Scaled normal",
                    hjust = 0, color = clr["Scaled normal"], size = 4) +
  ggplot2::annotate("text", x = 0.175, y = 0.90, label = "Horseshoe",
                    hjust = 1, color = clr["RHS"], size = 4) +
  ggplot2::annotate("text", x = 0.23,  y = 0.90, label = "R2D2",
                    hjust = 0, color = clr["R2D2"], size = 4)

print(fig_28_6)
ggplot2::ggsave("figs/Fig-28.6.svg", plot = fig_28_6,
                width = 6, height = 4, device = svg)
log_result("Fig-28.6.svg saved")

# LOO comparison table
log_result("\n=== LOO comparison: four priors ===")
log_result("Book target: R2D2 best; Wide normal elpd_diff ≈ -4.1 (se 2.6)")
loo_tab <- brms::loo_compare(list(
  `Wide normal`  = loo::loo(fitm_n1),
  `Scaled normal` = loo::loo(fitm_n2),
  `Horseshoe`    = loo::loo(fitm_hs),
  `R2D2`         = loo::loo(fitm)
))
capture_result(loo_tab, label = "loo_compare: four priors")
r2d2_row  <- which(rownames(loo_tab) == "R2D2")
wn_row    <- which(rownames(loo_tab) == "Wide normal")
if (length(r2d2_row) > 0 && loo_tab[r2d2_row, "elpd_diff"] != 0)
  log_result("NOTE: R2D2 is not top model — investigate")
if (length(wn_row) > 0 && abs(loo_tab[wn_row, "elpd_diff"] - (-4.1)) > 2)
  log_result("NOTE: Wide normal elpd_diff ", round(loo_tab[wn_row, "elpd_diff"], 1),
             " differs from book -4.1 — seed/version difference expected")

# =============================================================================
# Refit fitm at iter=2000 for marginal posterior plots — §7 of book
# =============================================================================

if (file.exists("fitm_2000.rds")) {
  fitm_2000 <- readRDS("fitm_2000.rds")
  log_result("fitm_2000 loaded from cache")
} else {
  fitm_2000 <- stats::update(fitm, iter = 2000, warmup = 1000, refresh = 0)
  saveRDS(fitm_2000, "fitm_2000.rds")
  log_result("fitm_2000 (R2D2, iter=2000 refit) fitted and saved")
}
log_result("--- Diagnostics: fitm_2000 ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitm_2000),       na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitm_2000), na.rm = TRUE), 4))

# Fig-28.7: marginal posteriors, R2D2 prior (iter=2000 fit)
drawsm <- posterior::as_draws_df(fitm_2000,
            variable = paste0("b_", predictors)) |>
          posterior::set_variables(predictors)
stopifnot(length(posterior::variables(drawsm)) == 26)

fig_28_7 <- bayesplot::mcmc_areas(drawsm, prob_outer = 0.98,
                                   area_method = "scaled height") +
            ggplot2::xlim(c(-1.5, 1.5)) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Figure 28.7: Marginal posteriors, R2D2 prior")
fig_28_7 <- fig_28_7 +
            ggplot2::scale_y_discrete(limits = rev(levels(fig_28_7$data$parameter)))

print(fig_28_7)
ggplot2::ggsave("figs/Fig-28.7.svg", plot = fig_28_7,
                width = 8, height = 6, device = svg)
log_result("Fig-28.7.svg saved")

# Fig-28.8: bivariate scatter Fedu vs Medu
fig_28_8 <- bayesplot::mcmc_scatter(drawsm, pars = c("Fedu", "Medu"),
                                     size = 1, alpha = 0.1) +
            bayesplot::vline_0(linetype = "dashed") +
            bayesplot::hline_0(linetype = "dashed") +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Figure 28.8: Fedu vs Medu bivariate posterior")

print(fig_28_8)
ggplot2::ggsave("figs/Fig-28.8.svg", plot = fig_28_8,
                width = 4, height = 4, device = svg)
log_result("Fig-28.8.svg saved")

# Fig-28.9: PPC histogram (5 draws)
# PPC METHOD NOTE
# Standard bayesplot::ppc_dens_overlay() would be appropriate for a normal model, but
# the book uses pp_check(type="hist", ndraws=5) to show the bounded outcome
# (scores 0-20) vs the normal model's unbounded predictions. This directly
# illustrates the model's main limitation. We follow the book exactly.
y_mat    <- student_Gmat$Gmat
y_pp_mat <- brms::posterior_predict(fitm_2000, ndraws = 5)

fig_28_9 <- bayesplot::ppc_hist(y_mat, y_pp_mat) +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Figure 28.9: PPC histogram (5 draws)")

print(fig_28_9)
ggplot2::ggsave("figs/Fig-28.9.svg", plot = fig_28_9,
                width = 6, height = 4, device = svg)
log_result("Fig-28.9.svg saved")

# Fig-28.10: LOO-PIT-ECDF
fig_28_10 <- bayesplot::ppc_loo_pit_ecdf(
               y    = y_mat,
               yrep = brms::posterior_predict(fitm_2000),
               lw   = stats::weights(loo::loo(fitm_2000, save_psis = TRUE)$psis_object),
               method = "correlated") +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.10: LOO-PIT-ECDF")

print(fig_28_10)
ggplot2::ggsave("figs/Fig-28.10.svg", plot = fig_28_10,
                width = 5, height = 5, device = svg)
log_result("Fig-28.10.svg saved")

# =============================================================================
# Projection predictive variable selection — Math (fast search) — §9.1 of book
# =============================================================================

if (file.exists("vselm_fast.rds")) {
  vselm_fast <- readRDS("vselm_fast.rds")
  log_result("vselm_fast loaded from cache")
} else {
  vselm_fast <- projpred::cv_varsel(
    fitm_2000,
    nterms_max      = 27,
    validate_search = FALSE
  )
  saveRDS(vselm_fast, "vselm_fast.rds")
  log_result("vselm_fast fitted and saved")
}

# Fig-28.11: fast search path
fig_28_11 <- plot(vselm_fast,
                  stats           = c("elpd", "R2"),
                  deltas          = "mixed",
                  text_angle      = 45,
                  alpha           = 0.1,
                  size_position   = "primary_x_top",
                  show_cv_proportions = FALSE) +
             ggplot2::geom_vline(xintercept = seq(0, 25, by = 5),
                                 colour = "black", alpha = 0.1) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.11: Math fast variable selection path")

print(fig_28_11)
ggplot2::ggsave("figs/Fig-28.11.svg", plot = fig_28_11,
                width = 8, height = 6, device = svg)
log_result("Fig-28.11.svg saved")
log_result("NOTE: fast search (validate_search=FALSE) may overfit — use vselm for reliable size")

# =============================================================================
# Math validated variable selection + projected posterior — §9.1 of book
# =============================================================================

if (file.exists("vselm.rds")) {
  vselm <- readRDS("vselm.rds")
  log_result("vselm loaded from cache")
} else {
  doFuture::registerDoFuture()
  future::plan(future::multisession,
               workers = getOption("mc.cores", default = 1))
  vselm <- projpred::cv_varsel(
    fitm_2000,
    nterms_max      = 10,
    validate_search = TRUE,
    refit_prj       = TRUE,
    nloo            = 50,
    parallel        = TRUE,
    verbose         = FALSE
  )
  future::plan(future::sequential)
  saveRDS(vselm, "vselm.rds")
  log_result("vselm fitted and saved")
}

nselm <- projpred::suggest_size(vselm)
log_result("suggest_size(vselm) = ", nselm, "  (book target: 4)")
if (nselm != 4)
  log_result("NOTE: suggest_size = ", nselm, " differs from book 4 — seed/nloo difference expected")

rankm <- projpred::ranking(vselm, nterms = nselm)
log_result("Top ", nselm, " math predictors: ",
           paste(rankm$fulldata[seq_len(nselm)], collapse = ", "))

# Fig-28.12: validated search path
fig_28_12 <- plot(vselm,
                  stats               = c("elpd", "R2"),
                  deltas              = "mixed",
                  text_angle          = 45,
                  alpha               = 0.1,
                  size_position       = "primary_x_top",
                  show_cv_proportions = FALSE) +
             ggplot2::geom_vline(xintercept = seq(0, 10, by = 5),
                                 colour = "black", alpha = 0.1) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.12: Math validated variable selection path")

print(fig_28_12)
ggplot2::ggsave("figs/Fig-28.12.svg", plot = fig_28_12,
                width = 8, height = 5, device = svg)
log_result("Fig-28.12.svg saved")

# Projected posterior for selected model
projm <- projpred::project(vselm, nterms = nselm)
drawsm_proj <- posterior::as_draws_df(projm) |>
               posterior::subset_draws(
                 variable = paste0("b_", rankm$fulldata[seq_len(nselm)])) |>
               posterior::set_variables(rankm$fulldata[seq_len(nselm)])

# Fig-28.13: projected posterior marginals
fig_28_13 <- bayesplot::mcmc_areas(drawsm_proj, prob_outer = 0.98,
                                    area_method = "scaled height") +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = paste0("Figure 28.13: Math projected posterior (",
                                          nselm, " predictors)"))

print(fig_28_13)
ggplot2::ggsave("figs/Fig-28.13.svg", plot = fig_28_13,
                width = 8, height = 3, device = svg)
log_result("Fig-28.13.svg saved")

# Fig-28.14: search stability (cv_proportions)
cv_prop_m  <- projpred::cv_proportions(rankm, cumulate = TRUE)
fig_28_14  <- plot(cv_prop_m) +
              ggplot2::theme_minimal() +
              ggplot2::labs(title = "Figure 28.14: Math search stability (cv_proportions)")

print(fig_28_14)
ggplot2::ggsave("figs/Fig-28.14.svg", plot = fig_28_14,
                width = 8, height = 5, device = svg)
log_result("Fig-28.14.svg saved")

# =============================================================================
# Portuguese model (fitp) — R2D2 prior — §9.2 of book
# =============================================================================

if (file.exists("fitp.rds")) {
  fitp <- readRDS("fitp.rds")
  log_result("fitp loaded from cache")
} else {
  fitp <- brms::brm(Gpor ~ ., data = studentstd_Gpor,
    prior     = brms::prior(R2D2(mean_R2 = 1/3, prec_R2 = 3,
                                 cons_D2 = 1/2), class = b),
    warmup    = 1000, iter = 2000, refresh = 0,
    seed      = SEED, save_pars = brms::save_pars(all = TRUE))
  saveRDS(fitp, "fitp.rds")
  log_result("fitp fitted and saved")
}

log_result("--- Diagnostics: fitp ---")
log_result("  Rhat_max    : ", round(max(brms::rhat(fitp),       na.rm = TRUE), 4))
log_result("  neff_min    : ", round(min(brms::neff_ratio(fitp), na.rm = TRUE), 4))
log_result("  nobs(fitp)  : ", stats::nobs(fitp),
           "  (book reports 657 rows; complete cases only)")
if (stats::nobs(fitp) != 657)
  log_result("# INVESTIGATE: fitp fitted on ", stats::nobs(fitp),
             " rows, not 657. Portuguese data has NAs — brms drops incomplete cases.")

# LOO for fitp
fitp <- brms::add_criterion(fitp, criterion = "loo",
                             save_psis = TRUE)
saveRDS(fitp, "fitp.rds")

log_result("\n=== fitp: Portuguese R2D2 model ===")
log_result("Book target bayes_R2: Estimate ≈ 0.33")
log_result("Book target loo_R2:   Estimate ≈ 0.28")
capture_result(round(as.data.frame(bayes_R2_paper(fitp)), 2),
               label = "bayes_R2_paper(fitp)")
capture_result(round(as.data.frame(brms::loo_R2(fitp)), 2),
               label = "loo_R2(fitp)")
br2p <- bayes_R2_paper(fitp)["R2", "Estimate"]
lr2p <- brms::loo_R2(fitp)["R2", "Estimate"]
if (abs(br2p - 0.33) > 0.05)
  log_result("NOTE: bayes_R2 ", round(br2p, 3), " differs from book 0.33")
if (abs(lr2p - 0.28) > 0.05)
  log_result("NOTE: loo_R2 ", round(lr2p, 3), " differs from book 0.28")

# Fig-28.15: marginal posteriors, Portuguese R2D2 model
drawsp <- posterior::as_draws_df(fitp,
            variable = paste0("b_", predictors)) |>
          posterior::set_variables(predictors)
stopifnot(length(posterior::variables(drawsp)) == 26)

fig_28_15 <- bayesplot::mcmc_areas(drawsp, prob_outer = 0.98,
                                    area_method = "scaled height") +
             ggplot2::xlim(c(-1.5, 1.5)) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.15: Portuguese marginal posteriors, R2D2 prior")
fig_28_15 <- fig_28_15 +
             ggplot2::scale_y_discrete(limits = rev(levels(fig_28_15$data$parameter)))

print(fig_28_15)
ggplot2::ggsave("figs/Fig-28.15.svg", plot = fig_28_15,
                width = 8, height = 6, device = svg)
log_result("Fig-28.15.svg saved")

# =============================================================================
# Portuguese variable selection + projected posterior — §9.2 of book
# =============================================================================

# Fast search
if (file.exists("vselp_fast.rds")) {
  vselp_fast <- readRDS("vselp_fast.rds")
  log_result("vselp_fast loaded from cache")
} else {
  vselp_fast <- projpred::cv_varsel(
    fitp,
    nterms_max      = 27,
    validate_search = FALSE
  )
  saveRDS(vselp_fast, "vselp_fast.rds")
  log_result("vselp_fast fitted and saved")
}

# Fig-28.16: Portuguese fast search path
fig_28_16 <- plot(vselp_fast,
                  stats               = c("elpd", "R2"),
                  deltas              = "mixed",
                  text_angle          = 45,
                  alpha               = 0.1,
                  size_position       = "primary_x_top",
                  show_cv_proportions = FALSE) +
             ggplot2::geom_vline(xintercept = seq(0, 25, by = 5),
                                 colour = "black", alpha = 0.1) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.16: Portuguese fast variable selection path")

print(fig_28_16)
ggplot2::ggsave("figs/Fig-28.16.svg", plot = fig_28_16,
                width = 8, height = 6, device = svg)
log_result("Fig-28.16.svg saved")

# Validated search
if (file.exists("vselp.rds")) {
  vselp <- readRDS("vselp.rds")
  log_result("vselp loaded from cache")
} else {
  doFuture::registerDoFuture()
  future::plan(future::multisession,
               workers = getOption("mc.cores", default = 1))
  vselp <- projpred::cv_varsel(
    fitp,
    nterms_max      = 10,
    validate_search = TRUE,
    refit_prj       = TRUE,
    nloo            = 50,
    parallel        = TRUE,
    verbose         = FALSE
  )
  future::plan(future::sequential)
  saveRDS(vselp, "vselp.rds")
  log_result("vselp fitted and saved")
}

nselp <- projpred::suggest_size(vselp)
log_result("suggest_size(vselp) = ", nselp,
           "  (book text says 7; code output says 8)")
if (!nselp %in% c(7L, 8L))
  log_result("NOTE: suggest_size = ", nselp,
             " differs from both book values (7 and 8) — seed/nloo difference")

rankp <- projpred::ranking(vselp, nterms = nselp)
log_result("Top ", nselp, " Portuguese predictors: ",
           paste(rankp$fulldata[seq_len(nselp)], collapse = ", "))

# Fig-28.17: Portuguese validated search path
fig_28_17 <- plot(vselp,
                  stats               = c("elpd", "R2"),
                  deltas              = "mixed",
                  text_angle          = 45,
                  alpha               = 0.1,
                  size_position       = "primary_x_top",
                  show_cv_proportions = FALSE) +
             ggplot2::geom_vline(xintercept = seq(0, 10, by = 5),
                                 colour = "black", alpha = 0.1) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.17: Portuguese validated variable selection path")

print(fig_28_17)
ggplot2::ggsave("figs/Fig-28.17.svg", plot = fig_28_17,
                width = 8, height = 5, device = svg)
log_result("Fig-28.17.svg saved")

# Projected posterior for selected Portuguese model
projp <- projpred::project(vselp, nterms = nselp)
drawsp_proj <- posterior::as_draws_df(projp) |>
               posterior::subset_draws(
                 variable = paste0("b_", rankp$fulldata[seq_len(nselp)])) |>
               posterior::set_variables(rankp$fulldata[seq_len(nselp)])

# Fig-28.18: Portuguese projected posterior marginals
fig_28_18 <- bayesplot::mcmc_areas(drawsp_proj, prob_outer = 0.98,
                                    area_method = "scaled height") +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = paste0("Figure 28.18: Portuguese projected posterior (",
                                          nselp, " predictors)"))

print(fig_28_18)
ggplot2::ggsave("figs/Fig-28.18.svg", plot = fig_28_18,
                width = 8, height = 3, device = svg)
log_result("Fig-28.18.svg saved")

# Fig-28.19: Portuguese search stability
cv_prop_p <- projpred::cv_proportions(rankp, cumulate = TRUE)
fig_28_19 <- plot(cv_prop_p) +
             ggplot2::theme_minimal() +
             ggplot2::labs(title = "Figure 28.19: Portuguese search stability (cv_proportions)")

print(fig_28_19)
ggplot2::ggsave("figs/Fig-28.19.svg", plot = fig_28_19,
                width = 8, height = 5, device = svg)
log_result("Fig-28.19.svg saved")

# =============================================================================
# Save and wrap up
# =============================================================================

wf <- init_workflow(mode = "practice", stage = "explore")
wf$ppc_complete <- TRUE
wf$loo_complete <- TRUE
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

log_result("\n=== Ch 28 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: fitm_u.rds, fitm_n1.rds, fitm_n1p.rds, fitm_n1pt.rds,")
log_result("       fitm_n2.rds, fitm_n2p.rds, fitm_hs.rds, fitm_hsp.rds,")
log_result("       fitm.rds, fitmp.rds, fitm_2000.rds, fitp.rds,")
log_result("       vselm_fast.rds, vselm.rds, vselp_fast.rds, vselp.rds,")
log_result("       wf_final.rds")
log_result("Figures saved to figs/: Fig-28.1.svg through Fig-28.19.svg")
close(results_con)
