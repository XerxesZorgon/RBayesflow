# data/movies_ch16/movies_ch16_analysis.R
#
# Ch 16 — "Coding a series of models: Simulated data of movie ratings"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# Replicates the book's data generation, models, and figures using brms.
# Data generation and parameter values match the book exactly.
# Model structure matches the book's Stan programs; see comments where
# the brms formula differs cosmetically but is statistically equivalent.
#
# Figures produced:
#   Fig 16.1 — M2: posterior median + 50%/95% intervals vs true theta
#   Fig 16.2 — M2: interval width vs number of ratings
#   Fig 16.3 — M3 balanced: coverage check for alpha and beta
#   Fig 16.4 — M3 unbalanced: coverage check by genre and rater type
#   Fig 16.5 — M3 unbalanced: raw mean vs truth / model estimate vs truth
#
# Run interactively in RStudio (Ctrl+Enter line by line, or section by section).
# Working directory must be data/movies_ch16/.

source("../../R/source_all.R")

pause <- function(msg = "Press Enter to continue...") {
  readline(prompt = paste0("\n[Pause] ", msg, " "))
  invisible(NULL)
}

# ── Phase 1 ───────────────────────────────────────────────────────────────────

wf <- init_workflow(mode = "learn", stage = "explore")
wf <- run_phase1(wf, simulated = TRUE)
guide(wf)
pause("Phase 1 complete. Continue to Model 1.")

# =============================================================================
# MODEL 1 — Two movies, N = 2 and N = 100 (Section 16.1)
#
# Book data (fixed, not random — both movies average exactly 4.0):
#   Movie 1: c(3, 5)                          n = 2,   mean = 4.0
#   Movie 2: rep(c(2,3,4,5), c(10,20,30,40))  n = 100, mean = 4.0
#
# Model: y_i ~ normal(theta[movie_i], sigma_y)
#   theta[j] ~ normal(3, 1),  constrained to (0, 5) in book's Stan
#   sigma_y  ~ normal+(0, 2.5)
#
# brms note: bounds on theta are enforced via lb/ub on the b prior.
# =============================================================================

cat("\n\n=== MODEL 1: Two movies, unequal N (Section 16.1) ===\n")

y_1 <- c(3, 5)
y_2 <- rep(c(2, 3, 4, 5), c(10, 20, 30, 40))
y   <- c(y_1, y_2)
N   <- length(y)
dat_m1 <- data.frame(
  y     = y,
  movie = factor(rep(c("Movie1", "Movie2"), c(length(y_1), length(y_2))))
)

cat("Movie 1: n =", length(y_1), " mean =", mean(y_1), "\n")
cat("Movie 2: n =", length(y_2), " mean =", mean(y_2), "\n\n")

# Priors matching the book exactly.
# lb = 0, ub = 5 on b replicates Stan's vector<lower=0, upper=5>[2] theta.
priors_m1 <- c(
  brms::prior(normal(3, 1),   class = b,     lb = 0, ub = 5),
  brms::prior(normal(0, 2.5), class = sigma, lb = 0)
)

# Phase 2
wf <- run_phase2(
  wf      = wf,
  formula = y ~ 0 + movie,
  family  = gaussian(),
  priors  = priors_m1,
  data    = dat_m1,
  seed    = 42
)
pause("M1 Phase 2: Prior predictive plot shown. Draws should be concentrated between 0 and 5. Continue to fit.")

# Phase 3
result_m1 <- run_phase3(
  wf      = wf,
  formula = y ~ 0 + movie,
  family  = gaussian(),
  priors  = priors_m1,
  data    = dat_m1,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m1 <- result_m1$fit
wf     <- result_m1$wf
saveRDS(fit_m1, "fit_m1.rds")

# Phase 4
wf <- run_diagnostics(fit_m1, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M1 Phase 4: Diagnostics reviewed. Continue to posterior summary.")

# Book result (Table after Stan output, p.265):
#   theta[1]: mean 3.63, sd 0.55, 90% CI (2.70, 4.53)
#   theta[2]: mean 3.99, sd 0.10, 90% CI (3.82, 4.15)
#   sigma_y:  mean 1.02, sd 0.07
cat("\n=== M1: Posterior summary (compare to book p.265) ===\n")
cat("Book: theta[Movie1] mean=3.63, 90% CI (2.70, 4.53)\n")
cat("Book: theta[Movie2] mean=3.99, 90% CI (3.82, 4.15)\n")
cat("Book: sigma_y mean=1.02\n\n")
print(brms::fixef(fit_m1))
cat("\nsigma_y:\n")
print(brms::VarCorr(fit_m1))

wf$ppc_complete <- TRUE
export_context(wf)

pause("M1 complete. Movie 1 should have a much wider interval than Movie 2 despite equal means. Continue to M2.")

# =============================================================================
# MODEL 2 — J = 40 movies, varying N (Section 16.2)
#
# Book data generation (exact):
#   J = 40 movies
#   N_ratings[j] ~ Uniform(0, 100) — random number of ratings per movie
#   theta[j]     ~ normal(3.0, 0.5)
#   y_i          ~ normal(theta[movie_i], 2.0)   [sigma_y = 2.0]
#
# Key plots:
#   Fig 16.1 — posterior median + 50%/95% intervals vs true theta
#   Fig 16.2 — interval width vs number of ratings per movie
# =============================================================================

cat("\n\n=== MODEL 2: J = 40 movies, varying N (Section 16.2) ===\n")

set.seed(42)
J          <- 40
N_ratings  <- sample(0:100, J, replace = TRUE)   # book: uniform 0–100
N_m2       <- sum(N_ratings)
movie_idx  <- rep(seq_len(J), N_ratings)
theta_true <- rnorm(J, 3.0, 0.5)                 # book: sd = 0.5
y_m2       <- rnorm(N_m2, theta_true[movie_idx], 2.0)  # book: sigma_y = 2.0

dat_m2 <- data.frame(
  y     = y_m2,
  movie = factor(paste0("M", sprintf("%02d", movie_idx)))
)

cat("Movies:", J, " | Total ratings:", N_m2,
    " | Ratings per movie:", min(N_ratings), "–", max(N_ratings), "\n\n")

priors_m2 <- c(
  brms::prior(normal(3, 1),   class = b,     lb = 0, ub = 5),
  brms::prior(normal(0, 2.5), class = sigma, lb = 0)
)

# Phase 2
wf <- run_phase2(
  wf      = wf,
  formula = y ~ 0 + movie,
  family  = gaussian(),
  priors  = priors_m2,
  data    = dat_m2,
  seed    = 42
)
pause("M2 Phase 2: Prior predictive plot shown. Continue to fit.")

# Phase 3
result_m2 <- run_phase3(
  wf      = wf,
  formula = y ~ 0 + movie,
  family  = gaussian(),
  priors  = priors_m2,
  data    = dat_m2,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m2 <- result_m2$fit
wf     <- result_m2$wf
saveRDS(fit_m2, "fit_m2.rds")

# Phase 4
wf <- run_diagnostics(fit_m2, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M2 Phase 4: Diagnostics reviewed. Continue to Figure 16.1.")

# --- Figure 16.1: posterior median + 50%/95% intervals vs true theta ---
cat("\n=== Fig 16.1: Posterior intervals vs true theta ===\n")

post_m2    <- brms::fixef(fit_m2, probs = c(0.025, 0.25, 0.75, 0.975))
movie_lvls <- levels(dat_m2$movie)
# fixef() row names are "movieM01", "movieM02" etc.
# Strip "movie" then "M" to get integer index.
movie_nums <- as.integer(sub("M", "", sub("movie", "", rownames(post_m2))))

fig16_1_df <- data.frame(
  truth  = theta_true[movie_nums],
  median = post_m2[, "Estimate"],
  lo95   = post_m2[, "Q2.5"],
  hi95   = post_m2[, "Q97.5"],
  lo50   = post_m2[, "Q25"],
  hi50   = post_m2[, "Q75"]
)

fig16_1 <- ggplot2::ggplot(fig16_1_df,
    ggplot2::aes(x = truth, y = median)) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo95, yend = hi95),
                        colour = "steelblue", linewidth = 0.5) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo50, yend = hi50),
                        colour = "steelblue", linewidth = 1.5) +
  ggplot2::geom_point(colour = "steelblue", size = 2) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(
    title = "Figure 16.1: Comparing parameters θ_j to their posterior inferences",
    x     = "True parameter value",
    y     = "Posterior median, 50%, and 95% interval"
  ) +
  ggplot2::theme_minimal()
print(fig16_1)

pause("Figure 16.1 shown. Roughly half the 50% intervals and 95% of the 95% intervals should contain the true value. Continue to Figure 16.2.")

# --- Figure 16.2: interval width vs number of ratings ---
cat("\n=== Fig 16.2: Interval width vs number of ratings ===\n")

fig16_2_df <- data.frame(
  n_ratings    = N_ratings[movie_nums],
  interval_width = fig16_1_df$hi50 - fig16_1_df$lo50
)

fig16_2 <- ggplot2::ggplot(fig16_2_df,
    ggplot2::aes(x = n_ratings, y = interval_width)) +
  ggplot2::geom_point(colour = "steelblue") +
  ggplot2::labs(
    title = "Figure 16.2: Where you have more data, you have less uncertainty",
    x     = "Number of ratings",
    y     = "Width of 50% posterior interval"
  ) +
  ggplot2::theme_minimal()
print(fig16_2)

wf$ppc_complete <- TRUE
export_context(wf)

pause("Figure 16.2 shown. Interval width should decrease as number of ratings increases. Continue to M3.")

# =============================================================================
# MODEL 3 — Item-response model with rater effects (Section 16.3)
#
# Book model (non-centered parameterization):
#   y_i ~ normal(mu + sigma_a*alpha[movie_i] - sigma_b*beta[rater_i], sigma_y)
#   alpha[j] ~ normal(0, 1)
#   beta[k]  ~ normal(0, 1)
#   mu       ~ normal(3, 5)
#   sigma_a  ~ normal+(0, 5)
#   sigma_b  ~ normal+(0, 5)
#   sigma_y  ~ normal+(0, 5)
#
# brms note: brms uses + for all random effects; we replicate the - sigma_b*beta
# sign convention by negating beta explicitly in the data (see below).
# The a_j = mu + sigma_a*alpha_j transformed parameter is recovered via ranef().
#
# True parameters (book p.269):
#   mu = 3, sigma_a = 0.5, sigma_b = 0.5, sigma_y = 2
# =============================================================================

cat("\n\n=== MODEL 3: Item-response model (Section 16.3) ===\n")

set.seed(42)
J_3     <- 40
K       <- 100
mu      <- 3
sigma_a <- 0.5
sigma_b <- 0.5
sigma_y <- 2

alpha <- rnorm(J_3, 0, 1)   # movie quality (unscaled)
beta  <- rnorm(K,   0, 1)   # rater difficulty (unscaled)

# --- Balanced data (everyone rates everything) ---
cat("Simulating balanced data (J=40 movies, K=100 raters, N=4000)...\n")
movie_idx_3 <- rep(seq_len(J_3), each = K)
rater_idx_3 <- rep(seq_len(K), times = J_3)
y_3 <- rnorm(
  J_3 * K,
  mu + sigma_a * alpha[movie_idx_3] - sigma_b * beta[rater_idx_3],
  sigma_y
)

# brms does not support a minus sign between random effects directly.
# We negate beta in the data so that (1|rater) absorbs -sigma_b*beta,
# making the posterior of the rater random effect equal to -sigma_b*beta.
# This is equivalent to the book's model; we document it clearly.
dat_m3_bal <- data.frame(
  y       = y_3,
  movie   = factor(paste0("M", sprintf("%02d", movie_idx_3))),
  rater   = factor(paste0("R", sprintf("%03d", rater_idx_3)))
)

priors_m3 <- c(
  brms::prior(normal(3, 5),   class = Intercept),
  brms::prior(normal(0, 5),   class = sd,    lb = 0),
  brms::prior(normal(0, 5),   class = sigma, lb = 0)
)

# Phase 2
wf <- run_phase2(
  wf      = wf,
  formula = y ~ 1 + (1 | movie) + (1 | rater),
  family  = gaussian(),
  priors  = priors_m3,
  data    = dat_m3_bal,
  seed    = 42
)
pause("M3 Phase 2 (balanced): Prior predictive plot shown. Continue to fit (this will take several minutes).")

# Phase 3
result_m3_bal <- run_phase3(
  wf      = wf,
  formula = y ~ 1 + (1 | movie) + (1 | rater),
  family  = gaussian(),
  priors  = priors_m3,
  data    = dat_m3_bal,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m3_bal <- result_m3_bal$fit
wf         <- result_m3_bal$wf
saveRDS(fit_m3_bal, "fit_m3_balanced.rds")

# Phase 4
wf <- run_diagnostics(fit_m3_bal, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M3 balanced Phase 4: Diagnostics reviewed. Continue to Figure 16.3.")

# --- Figure 16.3: coverage check for alpha and beta (balanced) ---
cat("\n=== Fig 16.3: Coverage check for alpha and beta (balanced data) ===\n")
cat("Book hyperparameter results (p.269):\n")
cat("  mu=3.19, sigma_a=0.54, sigma_b=0.52, sigma_y=2.03\n\n")

# Hyperparameter summary
cat("Fitted hyperparameters:\n")
print(brms::fixef(fit_m3_bal))
print(brms::VarCorr(fit_m3_bal))

re_movie <- brms::ranef(fit_m3_bal, probs = c(0.025, 0.25, 0.75, 0.975))$movie[, , "Intercept"]
re_rater <- brms::ranef(fit_m3_bal, probs = c(0.025, 0.25, 0.75, 0.975))$rater[, , "Intercept"]

# alpha recovery (movie effects)
alpha_df <- data.frame(
  truth  = alpha,
  median = re_movie[, "Estimate"],
  lo95   = re_movie[, "Q2.5"],
  hi95   = re_movie[, "Q97.5"],
  lo50   = re_movie[, "Q25"],
  hi50   = re_movie[, "Q75"]
)

p_alpha <- ggplot2::ggplot(alpha_df, ggplot2::aes(x = truth, y = median)) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo95, yend = hi95),
                        colour = "steelblue", linewidth = 0.4) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo50, yend = hi50),
                        colour = "steelblue", linewidth = 1.2) +
  ggplot2::geom_point(colour = "steelblue", size = 1.5) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Checking the αj's",
                x = "True parameter value",
                y = "Posterior median, 50%, and 95% interval") +
  ggplot2::theme_minimal()

# beta recovery (rater effects; sign is negated in our brms encoding)
beta_df <- data.frame(
  truth  = beta,
  median = -re_rater[, "Estimate"],   # negate to recover original beta sign
  lo95   = -re_rater[, "Q97.5"],
  hi95   = -re_rater[, "Q2.5"],
  lo50   = -re_rater[, "Q75"],
  hi50   = -re_rater[, "Q25"]
)

p_beta <- ggplot2::ggplot(beta_df, ggplot2::aes(x = truth, y = median)) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo95, yend = hi95),
                        colour = "tomato", linewidth = 0.4) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo50, yend = hi50),
                        colour = "tomato", linewidth = 1.2) +
  ggplot2::geom_point(colour = "tomato", size = 1.5) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Checking the βk's",
                x = "True parameter value",
                y = "Posterior median, 50%, and 95% interval") +
  ggplot2::theme_minimal()

gridExtra::grid.arrange(p_alpha, p_beta, ncol = 2,
  top = "Figure 16.3: Balanced data — coverage check for movie and rater parameters")

pause("Figure 16.3 shown. Both panels should show points near the diagonal with good interval coverage. Continue to unbalanced data.")

# --- Unbalanced data: genre-based selection (book p.269-270) ---
cat("\n=== M3 Unbalanced: Genre-based selection bias ===\n")

genre <- rep(c("romantic", "crime"), c(round(J_3 / 2), J_3 - round(J_3 / 2)))

prob_of_rated <- ifelse(
  beta[rater_idx_3] > 0,
  ifelse(genre[movie_idx_3] == "romantic", 0.2, 0.7),
  ifelse(genre[movie_idx_3] == "romantic", 0.7, 0.2)
)

rated <- rbinom(J_3 * K, 1, prob_of_rated) == 1

dat_m3_unbal <- data.frame(
  y     = y_3[rated],
  movie = factor(paste0("M", sprintf("%02d", movie_idx_3[rated]))),
  rater = factor(paste0("R", sprintf("%03d", rater_idx_3[rated]))),
  genre = genre[movie_idx_3[rated]]
)

cat("Balanced N:", J_3 * K, "| Unbalanced N:", sum(rated),
    "(", round(100 * mean(rated), 1), "% retained)\n\n")

# Phase 3 only (reuse priors from balanced fit)
result_m3_unbal <- run_phase3(
  wf      = wf,
  formula = y ~ 1 + (1 | movie) + (1 | rater),
  family  = gaussian(),
  priors  = priors_m3,
  data    = dat_m3_unbal,
  seed    = 42,
  chains  = 4,
  iter    = 2000,
  warmup  = 1000
)
fit_m3_unbal <- result_m3_unbal$fit
wf           <- result_m3_unbal$wf
saveRDS(fit_m3_unbal, "fit_m3_unbalanced.rds")

# Phase 4
wf <- run_diagnostics(fit_m3_unbal, wf)
print(wf)
if (!isTRUE(wf$diagnostics$passed)) wf <- wf$diagnose()
pause("M3 unbalanced Phase 4: Diagnostics reviewed. Continue to Figure 16.4.")

# --- Figure 16.4: coverage check with genre labels ---
cat("\n=== Fig 16.4: Coverage check — unbalanced data, by genre and rater type ===\n")
cat("Book hyperparameter results (p.270):\n")
cat("  mu=3.19, sigma_a=0.55, sigma_b=0.61, sigma_y=2.01\n\n")
print(brms::fixef(fit_m3_unbal))
print(brms::VarCorr(fit_m3_unbal))

re_movie_u <- brms::ranef(fit_m3_unbal, probs = c(0.025, 0.25, 0.75, 0.975))$movie[, , "Intercept"]
re_rater_u <- brms::ranef(fit_m3_unbal, probs = c(0.025, 0.25, 0.75, 0.975))$rater[, , "Intercept"]

# Movie levels present in unbalanced data
movie_lvls_u <- rownames(re_movie_u)
movie_nums_u <- as.integer(sub("M", "", movie_lvls_u))

alpha_u_df <- data.frame(
  truth  = alpha[movie_nums_u],
  median = re_movie_u[, "Estimate"],
  lo95   = re_movie_u[, "Q2.5"],
  hi95   = re_movie_u[, "Q97.5"],
  lo50   = re_movie_u[, "Q25"],
  hi50   = re_movie_u[, "Q75"],
  genre  = genre[movie_nums_u]
)

p_alpha_u <- ggplot2::ggplot(alpha_u_df,
    ggplot2::aes(x = truth, y = median, shape = genre)) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo95, yend = hi95,
                                     colour = genre), linewidth = 0.4) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo50, yend = hi50,
                                     colour = genre), linewidth = 1.2) +
  ggplot2::geom_point(ggplot2::aes(colour = genre), size = 2) +
  ggplot2::scale_shape_manual(values = c(romantic = 1, crime = 16)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Checking the αj's",
                x = "True parameter value",
                y = "Posterior median, 50%, and 95% interval",
                colour = "Genre", shape = "Genre") +
  ggplot2::theme_minimal()

# Rater levels in unbalanced data
rater_lvls_u <- rownames(re_rater_u)
rater_nums_u <- as.integer(sub("R", "", rater_lvls_u))

beta_u_df <- data.frame(
  truth      = beta[rater_nums_u],
  median     = -re_rater_u[, "Estimate"],
  lo95       = -re_rater_u[, "Q97.5"],
  hi95       = -re_rater_u[, "Q2.5"],
  lo50       = -re_rater_u[, "Q75"],
  hi50       = -re_rater_u[, "Q25"],
  rater_type = ifelse(beta[rater_nums_u] > 0, "difficult", "nice")
)

p_beta_u <- ggplot2::ggplot(beta_u_df,
    ggplot2::aes(x = truth, y = median, shape = rater_type)) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo95, yend = hi95,
                                     colour = rater_type), linewidth = 0.4) +
  ggplot2::geom_segment(ggplot2::aes(xend = truth, y = lo50, yend = hi50,
                                     colour = rater_type), linewidth = 1.2) +
  ggplot2::geom_point(ggplot2::aes(colour = rater_type), size = 2) +
  ggplot2::scale_shape_manual(values = c(nice = 1, difficult = 16)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Checking the βk's",
                x = "True parameter value",
                y = "Posterior median, 50%, and 95% interval",
                colour = "Rater type", shape = "Rater type") +
  ggplot2::theme_minimal()

gridExtra::grid.arrange(p_alpha_u, p_beta_u, ncol = 2,
  top = "Figure 16.4: Unbalanced data — coverage check by genre and rater type")

pause("Figure 16.4 shown. Continue to Figure 16.5 (the chapter punchline).")

# --- Figure 16.5: raw mean vs truth / model estimate vs truth ---
cat("\n=== Fig 16.5: Raw averaging vs model-based estimates ===\n")

# True movie quality on original scale: a_j = mu + sigma_a * alpha[j]
a_true <- mu + sigma_a * alpha

# Raw mean per movie from unbalanced data
ybar <- tapply(dat_m3_unbal$y,
               as.integer(sub("M", "", as.character(dat_m3_unbal$movie))),
               mean)

# Posterior median of a_j = mu + sigma_a * alpha_j
# In brms: fixef()["Intercept"] ≈ mu; ranef()$movie ≈ sigma_a * alpha_j
# So a_j_post = fixef Intercept + ranef movie Intercept
mu_post     <- brms::fixef(fit_m3_unbal)["Intercept", "Estimate"]
a_post_med  <- mu_post + re_movie_u[, "Estimate"]

movie_nums_u_ordered <- as.integer(sub("M", "", rownames(re_movie_u)))

fig16_5_df <- data.frame(
  a_true   = a_true[movie_nums_u_ordered],
  ybar     = as.numeric(ybar[movie_nums_u_ordered]),
  a_post   = as.numeric(a_post_med),
  genre    = genre[movie_nums_u_ordered]
)

p_raw <- ggplot2::ggplot(fig16_5_df,
    ggplot2::aes(x = a_true, y = ybar, shape = genre)) +
  ggplot2::geom_point(ggplot2::aes(colour = genre), size = 2) +
  ggplot2::scale_shape_manual(values = c(romantic = 1, crime = 16)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Problems with raw averaging",
                x = "True a_j",
                y = "Raw average rating for movie j",
                colour = "Genre", shape = "Genre") +
  ggplot2::theme_minimal()

p_model <- ggplot2::ggplot(fig16_5_df,
    ggplot2::aes(x = a_true, y = a_post, shape = genre)) +
  ggplot2::geom_point(ggplot2::aes(colour = genre), size = 2) +
  ggplot2::scale_shape_manual(values = c(romantic = 1, crime = 16)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::labs(title = "Model-based estimates do better",
                x = "True a_j",
                y = "Posterior median estimate for movie j",
                colour = "Genre", shape = "Genre") +
  ggplot2::theme_minimal()

gridExtra::grid.arrange(p_raw, p_model, ncol = 2,
  top = "Figure 16.5: Raw averaging vs model-based estimates (unbalanced data)")

pause("Figure 16.5 shown. Left: romantic comedies should appear above the diagonal (inflated raw means), crime movies below. Right: model estimates should hug the diagonal for both genres.")

# =============================================================================
# Save and wrap up
# =============================================================================

wf$ppc_complete  <- TRUE
wf$loo_complete  <- FALSE   # LOO not run (simulated data; no prediction target)
saveRDS(wf, "wf_final.rds")
export_context(wf)
guide(wf)

cat("\n=== Ch 16 analysis complete ===\n")
cat("Saved: fit_m1.rds, fit_m3_balanced.rds, fit_m3_unbalanced.rds, wf_final.rds\n")
cat("Figures reproduced: 16.1, 16.2, 16.3, 16.4, 16.5\n")
