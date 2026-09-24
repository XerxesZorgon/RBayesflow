# data/sleep_ch17/sleep_ch17_analysis.R
#
# Ch 17 — "Prior specification for regression models: Reanalysis of a sleep study"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Replicates the book's core three models using the sleepstudy data (lme4).
# The chapter's focus is PRIOR SPECIFICATION, not model comparison. This
# script covers three models of increasing complexity, showing how the
# priors scale with the model structure.
#
# Data (from lme4 package):
#   sleepstudy: 18 subjects, 10 days of sleep deprivation each (Days 0-9).
#   Book pre-filters to Days >= 2 and re-centres Day 0 = first day of real
#   deprivation (the first two days are "adaptation"). N = 144 observations.
#
# Models covered:
#   fit1  Reaction ~ 1 + Days                          (Gaussian, no pooling)
#   fit3  Reaction ~ 1 + Days + (1 | Subject)          (Gaussian, varying intercept)
#   fit4  Reaction ~ 1 + Days + (1 + Days | Subject)   (Gaussian, varying int + slope, LKJ prior)
#
# Skipped (see write-up for pointers):
#   priorsense sensitivity (Figs 17.3, 17.4)
#   Lognormal, Student-t, log-t variants (fit_ln2, fit5, fit4t, fit4lt)
#   LOO-PIT model checking (Figs 17.10-17.12)
#
# Figures reproduced:
#   Fig 17.1  Data overview: Reaction vs Days, faceted by Subject
#   Fig 17.2  Prior predictive check for fit1 (via RBayesflow Phase 2 overlay)
#   Fig 17.6  Person-specific regression lines from fit4
#   Fig 17.7  Partial pooling: independent lm() vs fit4 estimates
#
# Book-target posterior summaries (used for success criteria):
#   fit1  Intercept 268.17, Days 11.34, sigma 51.21
#   fit3  Intercept 267.56, Days 11.41, sd(Intercept) 43.56, sigma 30.46
#   fit4  Intercept 267.54, Days 11.26, sd(Intercept) 31.72,
#         sd(Days) 7.25, cor(Intercept, Days) 0.21, sigma 25.93
#
# Run interactively in RStudio (Ctrl+Enter line by line, or section by section).
# Working directory must be data/sleep_ch17/.

source("../../R/source_all.R")

pause <- function(msg = "Press Enter to continue...") {
  readline(prompt = paste0("\n[Pause] ", msg, " "))
  invisible(NULL)
}

# ── Load and prepare data ─────────────────────────────────────────────────────

data("sleepstudy", package = "lme4")

# Book pre-filter: drop Days 0 and 1 (adaptation period) and re-centre so
# Day 0 = first day of real sleep deprivation. This gives N = 144.
sleepstudy       <- subset(sleepstudy, Days >= 2)
sleepstudy$Days  <- sleepstudy$Days - 2

cat("sleepstudy: N =", nrow(sleepstudy),
    " Subjects =", length(unique(sleepstudy$Subject)),
    " Days range:", min(sleepstudy$Days), "-", max(sleepstudy$Days), "\n\n")

# ── Phase 1 ───────────────────────────────────────────────────────────────────

wf <- init_workflow(mode = "learn", stage = "explore")

# Real data (not simulated). When run_phase1() offers non-Bayesian off-ramps,
# select "full_stan" — the book commits to a Bayesian analysis throughout in
# order to teach prior specification.
wf <- run_phase1(wf, simulated = FALSE)

# --- Figure 17.1: data overview (Reaction vs Days per Subject) ---
cat("\n=== Fig 17.1: Sleep deprivation data by subject ===\n")

fig17_1 <- ggplot2::ggplot(sleepstudy,
    ggplot2::aes(x = Days, y = Reaction)) +
  ggplot2::geom_point(colour = "steelblue") +
  ggplot2::facet_wrap(~ Subject, ncol = 6) +
  ggplot2::labs(
    title = "Figure 17.1: Sleep deprivation data for 18 subjects",
    x     = "Days of sleep deprivation",
    y     = "Reaction time (ms)"
  ) +
  ggplot2::theme_minimal()
print(fig17_1)

guide(wf)
pause("Phase 1 complete. Continue to Model 1.")

# =============================================================================
# MODEL 1 — Simple linear regression (Section 17.1)
#
#   y_i ~ normal(mu_i, sigma)
#   mu_i = b0 + b1 * Days_i
#
# Priors (from the book, page 277):
#   b0 (intercept):      normal(250, 100)
#   b1 (Days slope):     normal(0, 20)
#   sigma (residual sd): exponential(1/50)   =  exponential(0.02)
#
# Reasoning:
#   b0 ~ normal(250, 100) puts ~68% prior probability on average reaction time
#     (around the mean day, day 3.5) in (150, 350) ms — plausible.
#   b1 ~ normal(0, 20) puts ~95% prior probability on per-day change in
#     (-40, 40) ms/day — generous, rules out absurd slopes.
#   sigma ~ exponential(0.02) has mean 50 ms — the expected residual scale.
#
# Note: brms centres predictors internally, so "Intercept" during sampling is
# the mean-centred intercept. brms then reports the intercept on the original
# (uncentred) scale, so the book's target of 268.17 is what appears in the
# standard summary output — no back-transformation needed.
# =============================================================================

cat("\n\n=== MODEL 1: Simple linear regression (Section 17.1) ===\n")

priors_m1 <- c(
  brms::prior(normal(250, 100),  class = Intercept),
  brms::prior(normal(0, 20),     class = b),
  brms::prior(exponential(0.02), class = sigma)
)

# Phase 2 — prior predictive check
wf <- run_phase2(
  wf      = wf,
  formula = Reaction ~ 1 + Days,
  family  = gaussian(),
  priors  = priors_m1,
  data    = sleepstudy,
  seed    = 42
)
pause("M1 Phase 2: Prior predictive plot shown. Draws should mostly fall between 0 and 500 ms (compare book Fig 17.2). Continue to fit.")

# Phase 3 — fit
result_m1 <- run_phase3(
  wf      = wf,
  formula = Reaction ~ 1 + Days,
  family  = gaussian(),
  priors  = priors_m1,
  data    = sleepstudy,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m1 <- result_m1$fit
wf     <- result_m1$wf
saveRDS(fit_m1, "fit_m1.rds")

# Phase 4 — diagnostics
wf <- run_diagnostics(fit_m1, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M1 Phase 4: Diagnostics reviewed. Continue to posterior summary.")

# Book result (page 278):
#   Intercept: 268.17, sd 7.80, 95% CI (253.03, 283.38)
#   Days:       11.34, sd 1.88, 95% CI (  7.62,  15.10)
#   sigma:      51.21, sd 3.14, 95% CI ( 45.71,  57.85)
cat("\n=== M1: Posterior summary (compare to book p.278) ===\n")
cat("Book: Intercept 268.17, 95% CI (253.03, 283.38)\n")
cat("Book: Days       11.34, 95% CI (  7.62,  15.10)\n")
cat("Book: sigma      51.21, 95% CI ( 45.71,  57.85)\n\n")
print(brms::fixef(fit_m1))
cat("\nsigma:\n")
print(brms::VarCorr(fit_m1))

wf$ppc_complete <- TRUE
export_context(wf)

pause("M1 complete. The simple regression ignores that each subject contributed 8 observations — its standard error on the slope is inflated by that unmodeled within-subject dependence. Continue to M3.")

# =============================================================================
# MODEL 3 — Varying intercept (Section 17.2)
#
#   y_i ~ normal(mu_i, sigma)
#   mu_i = b0_j[i] + b1 * Days_i
#   b0_j ~ normal(b0, tau0)                [hierarchical prior on intercept]
#
# Non-centred parameterisation (what brms actually samples):
#   b0_j = b0 + tau0 * z0_j
#   z0_j ~ normal(0, 1)
#
# Priors (from the book, page 281):
#   b0    (overall intercept):                  normal(250, 100)   [same as fit1]
#   b1    (Days slope):                         normal(0, 20)      [same as fit1]
#   tau0  (between-subject sd of intercepts):   exponential(1/25)  =  exponential(0.04)
#   sigma (within-subject residual sd):         exponential(1/25)  =  exponential(0.04)
#
# Reasoning:
#   The single exponential(1/50) prior on sigma from fit1 has been "split"
#   into two exponential(1/25) priors (each with mean 25 ms) covering the
#   between-subject and within-subject variation separately.
# =============================================================================

cat("\n\n=== MODEL 3: Varying intercept (Section 17.2) ===\n")

priors_m3 <- c(
  brms::prior(normal(250, 100),  class = Intercept),
  brms::prior(normal(0, 20),     class = b),
  brms::prior(exponential(0.04), class = sd,    group = Subject),
  brms::prior(exponential(0.04), class = sigma)
)

# Phase 2 — prior predictive check
wf <- run_phase2(
  wf      = wf,
  formula = Reaction ~ 1 + Days + (1 | Subject),
  family  = gaussian(),
  priors  = priors_m3,
  data    = sleepstudy,
  seed    = 42
)
pause("M3 Phase 2: Prior predictive plot shown. Continue to fit.")

# Phase 3 — fit
result_m3 <- run_phase3(
  wf      = wf,
  formula = Reaction ~ 1 + Days + (1 | Subject),
  family  = gaussian(),
  priors  = priors_m3,
  data    = sleepstudy,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m3 <- result_m3$fit
wf     <- result_m3$wf
saveRDS(fit_m3, "fit_m3.rds")

# Phase 4 — diagnostics
wf <- run_diagnostics(fit_m3, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M3 Phase 4: Diagnostics reviewed. Continue to posterior summary.")

# Book result (page 281):
#   Intercept:                267.56, sd 11.16, 95% CI (245.89, 288.50)
#   Days:                      11.41, sd  1.13, 95% CI (  9.14,  13.57)
#   sd(Intercept) Subject:     43.56, sd  8.23, 95% CI ( 30.90,  62.07)
#   sigma:                     30.46, sd  2.00, 95% CI ( 26.95,  34.77)
cat("\n=== M3: Posterior summary (compare to book p.281) ===\n")
cat("Book: Intercept              267.56, 95% CI (245.89, 288.50)\n")
cat("Book: Days                    11.41, 95% CI (  9.14,  13.57)\n")
cat("Book: sd(Intercept) Subject   43.56, 95% CI ( 30.90,  62.07)\n")
cat("Book: sigma                   30.46, 95% CI ( 26.95,  34.77)\n\n")
print(brms::fixef(fit_m3))
print(brms::VarCorr(fit_m3))

wf$ppc_complete <- TRUE
export_context(wf)

pause("M3 complete. Note the slope's uncertainty tightened (sd ~1.1 vs ~1.9 in M1) once within-subject dependence was modelled. Continue to M4.")

# =============================================================================
# MODEL 4 — Varying intercept + varying slope + LKJ correlation prior (Section 17.2)
#
#   y_i ~ normal(mu_i, sigma)
#   mu_i = b0_j[i] + b1_j[i] * Days_i
#   (b0_j, b1_j) ~ MVN((b0, b1), Sigma)
#
# The 2x2 covariance matrix Sigma is decomposed into standard deviations
# (tau0, tau1) and a correlation matrix C, with C given an LKJ prior:
#
#   Sigma = diag(tau0, tau1) * C * diag(tau0, tau1)
#   C = [[1, rho], [rho, 1]]
#   C ~ LKJ(eta = 1)
#
# LKJ(1) is UNIFORM over 2x2 correlation matrices, so the prior on rho is
# uniform on (-1, 1). For higher-dimensional correlation matrices, LKJ(1)
# concentrates around zero — a side effect of the marginal being non-uniform
# even when the joint is (see book Fig 17.5).
#
# Priors (from the book, page 283):
#   b0:                                 normal(250, 100)
#   b1:                                 normal(0, 20)
#   tau0 = sd(Intercept) Subject:       exponential(1/25)  =  exponential(0.04)
#   tau1 = sd(Days) Subject:            exponential(1/10)  =  exponential(0.1)
#   cor(Intercept, Days) Subject:       LKJ(1)
#   sigma:                              exponential(1/25)  =  exponential(0.04)
# =============================================================================

cat("\n\n=== MODEL 4: Varying intercept + slope with LKJ prior (Section 17.2) ===\n")

priors_m4 <- c(
  brms::prior(normal(250, 100),  class = Intercept),
  brms::prior(normal(0, 20),     class = b),
  brms::prior(exponential(0.04), class = sd,    group = Subject, coef = Intercept),
  brms::prior(exponential(0.1),  class = sd,    group = Subject, coef = Days),
  brms::prior(lkj(1),            class = cor,   group = Subject),
  brms::prior(exponential(0.04), class = sigma)
)

# Phase 2 — prior predictive check
wf <- run_phase2(
  wf      = wf,
  formula = Reaction ~ 1 + Days + (1 + Days | Subject),
  family  = gaussian(),
  priors  = priors_m4,
  data    = sleepstudy,
  seed    = 42
)
pause("M4 Phase 2: Prior predictive plot shown. Continue to fit (this will take a bit longer than M3).")

# Phase 3 — fit
result_m4 <- run_phase3(
  wf      = wf,
  formula = Reaction ~ 1 + Days + (1 + Days | Subject),
  family  = gaussian(),
  priors  = priors_m4,
  data    = sleepstudy,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m4 <- result_m4$fit
wf     <- result_m4$wf
saveRDS(fit_m4, "fit_m4.rds")

# Phase 4 — diagnostics
wf <- run_diagnostics(fit_m4, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M4 Phase 4: Diagnostics reviewed. Continue to posterior summary.")

# Book result (page 283):
#   Intercept:                          267.54, sd  8.61, 95% CI (250.72, 284.84)
#   Days:                                11.26, sd  2.03, 95% CI (  7.18,  15.21)
#   sd(Intercept) Subject:               31.72, sd  7.34, 95% CI ( 19.96,  48.80)
#   sd(Days) Subject:                     7.25, sd  1.82, 95% CI (  4.23,  11.44)
#   cor(Intercept, Days) Subject:         0.21, sd  0.30, 95% CI ( -0.38,   0.78)
#   sigma:                               25.93, sd  1.81, 95% CI ( 22.70,  29.73)
cat("\n=== M4: Posterior summary (compare to book p.283) ===\n")
cat("Book: Intercept                        267.54, 95% CI (250.72, 284.84)\n")
cat("Book: Days                               11.26, 95% CI (  7.18,  15.21)\n")
cat("Book: sd(Intercept) Subject              31.72, 95% CI ( 19.96,  48.80)\n")
cat("Book: sd(Days) Subject                    7.25, 95% CI (  4.23,  11.44)\n")
cat("Book: cor(Intercept, Days) Subject        0.21, 95% CI ( -0.38,   0.78)\n")
cat("Book: sigma                              25.93, 95% CI ( 22.70,  29.73)\n\n")
print(brms::fixef(fit_m4))
print(brms::VarCorr(fit_m4))

wf$ppc_complete <- TRUE
export_context(wf)

pause("M4 complete. The overall slope's uncertainty widened (sd ~2.0 vs ~1.1 in M3) because between-subject variation in slopes is now modelled explicitly. Continue to Figure 17.6.")

# --- Figure 17.6: person-specific regression lines from fit4 ---
cat("\n=== Fig 17.6: Person-specific regression lines from M4 ===\n")

# For each subject, extract the posterior-median intercept and slope.
# brms::coef() returns a 3D array: [Subject, Estimate/lo/hi, Parameter].
# We grab the "Estimate" (posterior median) slice.
subj_coefs <- brms::coef(fit_m4)$Subject[, "Estimate", ]

subj_df <- data.frame(
  Subject   = rownames(subj_coefs),
  Intercept = subj_coefs[, "Intercept"],
  Slope     = subj_coefs[, "Days"]
)

# Build fitted lines by joining subject coefficients with the data grid.
line_df        <- merge(subj_df, sleepstudy, by = "Subject")
line_df$fitted <- line_df$Intercept + line_df$Slope * line_df$Days

fig17_6 <- ggplot2::ggplot(sleepstudy,
    ggplot2::aes(x = Days, y = Reaction)) +
  ggplot2::geom_point(colour = "steelblue") +
  ggplot2::geom_line(data = line_df,
                     ggplot2::aes(y = fitted),
                     colour = "steelblue", linewidth = 0.5) +
  ggplot2::facet_wrap(~ Subject, ncol = 6) +
  ggplot2::labs(
    title = "Figure 17.6: M4 person-specific regression lines",
    x     = "Days of sleep deprivation",
    y     = "Reaction time (ms)"
  ) +
  ggplot2::theme_minimal()
print(fig17_6)

pause("Figure 17.6 shown. Each subject's line is pulled toward the population mean by the multilevel prior. Continue to Figure 17.7.")

# --- Figure 17.7: partial pooling — independent lm() vs M4 estimates ---
cat("\n=== Fig 17.7: Partial pooling — no-pooling vs multilevel estimates ===\n")

# No-pooling estimates: fit a separate lm() per subject.
subjects <- levels(sleepstudy$Subject)
no_pool_list <- lapply(subjects, function(s) {
  d     <- sleepstudy[sleepstudy$Subject == s, ]
  fit_s <- lm(Reaction ~ Days, data = d)
  data.frame(
    Subject   = s,
    Intercept = unname(coef(fit_s)["(Intercept)"]),
    Slope     = unname(coef(fit_s)["Days"])
  )
})
no_pool_df       <- do.call(rbind, no_pool_list)
no_pool_df$model <- "Independent models"

# Multilevel (M4) estimates from subj_df above.
multi_df <- data.frame(
  Subject   = subj_df$Subject,
  Intercept = subj_df$Intercept,
  Slope     = subj_df$Slope,
  model     = "Joint multilevel model"
)

fig17_7_df <- rbind(no_pool_df, multi_df)

fig17_7 <- ggplot2::ggplot(fig17_7_df,
    ggplot2::aes(x = Intercept, y = Slope,
                 colour = model, shape = model)) +
  # Grey segment connecting each subject's two estimates so the pooling
  # motion is visible.
  ggplot2::geom_line(ggplot2::aes(group = Subject),
                     colour = "grey60", linewidth = 0.3,
                     show.legend = FALSE) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    data = subset(fig17_7_df, model == "Independent models"),
    ggplot2::aes(label = Subject),
    nudge_y = 0.7, size = 3, show.legend = FALSE
  ) +
  ggplot2::labs(
    title = "Figure 17.7: Partial pooling of subject-specific coefficients",
    x     = "Intercept estimate",
    y     = "Slope estimate",
    colour = NULL, shape = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "top")
print(fig17_7)

pause("Figure 17.7 shown. Grey segments show how each subject's estimate is pulled toward the group centre; extreme estimates are pulled the most.")

# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete <- TRUE
wf$loo_complete <- FALSE   # LOO not run — this script covers priors, not comparison
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

cat("\n=== Ch 17 analysis complete ===\n")
cat("Saved: fit_m1.rds, fit_m3.rds, fit_m4.rds, wf_final.rds\n")
cat("Figures reproduced: 17.1, 17.2 (via Phase 2 overlay), 17.6, 17.7\n")
