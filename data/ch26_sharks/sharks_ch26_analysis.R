# data/ch26_sharks/sharks_ch26_analysis.R
#
# Ch 26 — "Model building with latent variables: Markov models for animal movement"
# Gelman, Vehtari et al., Bayesian Workflow (2026)
#
# This script fits a sequence of 2-state hidden Markov models (HMMs) to
# white-shark GPS tracking data from Gansbaii, South Africa, replicating
# the analysis in Ch 26 of the Bayesian Workflow book.
#
# All three models use cmdstanr directly; brms/run_phase3() is NOT used
# because HMMs marginalise over discrete latent states and require a
# hand-coded Stan forward algorithm that brms cannot generate.
# RBayesflow's wf_state is still used for diagnostic tracking and
# export_context(), but phase 3 is bypassed per the rubric.
#
# v2 changes (after first run):
#   - step_turn_hmm.stan: state-specific priors on mu, tighter sigma prior.
#   - R script: Pathfinder initialisation (Lesson 3) for Models 1 and 2.
#   - Rhat 1.73 / ESS 6 in v1 was label-switching between chains; fixed
#     by separating mu priors so states cannot swap roles at initialisation.
#
# Data:
#   whiteshark_trackdata.RData — positional GPS tracks of white sharks
#   Source: Towner et al. (2016), "Sex-Specific and Individual Preferences
#           for Hunting Strategies in White Sharks", Functional Ecology 30:1397.
#   N = 4584 observations across 76 tracks (after exclusions)
#   Key variables: Long, Lat, SharksexTrackNo, CDB (chum indicator)
#   Acquired: download.file() from raw.githubusercontent.com in §Data acquisition
#
# Models covered:
#   fit_2stateHMM          : step_turn_hmm.stan — baseline 2-state HMM
#   fit_2stateHMM_covariates: step_turn_hmm_covariates.stan — tpm with covariates
#   fit_2stateHMM_tpmcov_crencp: step_turn_hmm_covariates_cre_ncp.stan
#                              — covariates + non-centred individual random effects
#
# Figures produced:
#   Fig 26.1  — WSF1 T1 and T10 tracks (longitude / latitude)
#   Fig 26.2  — MCMC marginal histograms by chain for mu (fit_2stateHMM)
#   Fig 26.3  — Step-length and turning-angle state-dependent distributions
#   Fig 26.4  — State-1 probability overlaid on WSF1 tracks
#   Fig 26.5  — State-1 probability time series with 95% CI for WSF1
#   Fig 26.6  — Pseudo-residual Q-Q plot (step lengths)
#   Fig 26.7  — State-decoding histograms: PPD, state decodings, chum
#   Fig 26.8  — (base-R) tpm entries vs time for male, no chum (model 3)
#   Fig 26.9  — (base-R) tpm entries vs time for female, no chum (model 3)
#   Fig 26.10 — Transition probabilities: male white shark, no chum (ggplot)
#   Fig 26.11 — Transition probabilities: female white shark, no chum (ggplot)
#   Fig 26.12 — Transition probabilities: WSF1 T1 individual (ggplot)
#
# Book-target posterior summaries (used for success criteria):
#   No specific book targets printed in the chapter; visual agreement with
#   Figures 1-10 in the book is the acceptance criterion.
#   mu[1] expected ~0.05-0.15 (area-restricted search / slow state)
#   mu[2] expected ~0.20-0.50 (directed travel / fast state)
#   tpm diagonal entries expected > 0.8 (persistence)
#
# Run interactively in RStudio (line by line or section by section).
# Working directory must be data/ch26_sharks/.

# ============================================================
## §3 Environment setup
# ============================================================
source("../../R/source_all.R")

library(tidyverse)
library(CircStats)   # dvm(), rvm()
library(patchwork)
library(lubridate)
library(moveHMM)
library(cmdstanr)

options(mc.cores = 4)

SEED <- 42

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

log_result("=== Ch 26 results log v2 — ", format(Sys.time()), " ===")
log_result("Stan-native chapter: brms/run_phase3() bypassed.")
log_result("Stan files: step_turn_hmm.stan, step_turn_hmm_covariates.stan,")
log_result("            step_turn_hmm_covariates_cre_ncp.stan")
log_result("v2: state-specific mu priors + Pathfinder init to fix Rhat 1.73 / ESS 6.")

# ============================================================
## §4 wf_state initialisation
# ============================================================
wf <- init_workflow(mode = "practice", stage = "explore")
log_result("RBayesflow wf_state initialised: mode=practice, stage=explore")

# ============================================================
## §5 Data acquisition
# ============================================================
data_path <- file.path("data", "whiteshark_trackdata.RData")

if (!file.exists(data_path)) {
  data_url <- paste0(
    "https://raw.githubusercontent.com/avehtari/Bayesian-Workflow/",
    "master/sharks/data/whiteshark_trackdata.RData"
  )
  message("Downloading whiteshark_trackdata.RData ...")
  download.file(data_url, destfile = data_path, mode = "wb")
  message("Download complete.")
} else {
  message("Using cached whiteshark_trackdata.RData")
}

load(data_path)
log_result("Data loaded: whiteshark_trackdata.RData")
log_result("  dim(sharks.HMMtracks.df): ",
           paste(dim(sharks.HMMtracks.df), collapse = " x "))

# ============================================================
## §6 Data preparation (Phase 1)
# ============================================================
sharks.HMMtracks.df$ID <- sharks.HMMtracks.df$SharksexTrackNo

moveHMM_wsdata <- moveHMM::prepData(
  sharks.HMMtracks.df[, c("ID", "Long", "Lat")],
  type       = c("LL"),
  coordNames = c("Long", "Lat")
)

sharks.HMMtracks.df$steplength <- moveHMM_wsdata$step
sharks.HMMtracks.df$turnang    <- moveHMM_wsdata$angle

ws_HMM_full <- dplyr::filter(
  sharks.HMMtracks.df,
  !SharksexTrackNo %in% c("WSF9 T4", "WSF9 T3 B ")
)

ws_HMM <- ws_HMM_full[, c("dateTime", "SharkName", "SharksexTrackNo",
                            "steplength", "turnang", "year", "month", "CDB")]

log_result("ws_HMM rows after exclusions: ", nrow(ws_HMM))
log_result("Unique tracks: ", length(unique(ws_HMM$SharksexTrackNo)))

ws_HMM$steplength[is.na(ws_HMM$steplength)]       <- -100
ws_HMM$steplength[which(ws_HMM$steplength > 1.5)] <- -100
ws_HMM$turnang[is.na(ws_HMM$turnang)]              <- -100

log_result("Missing step lengths encoded as -100: ",
           sum(ws_HMM$steplength < 0))

# ============================================================
## §7 Figure 26.1 — Raw tracks for WSF1 T1 and WSF1 T10
# ============================================================
fig26_1 <- ggplot2::ggplot(
  data = dplyr::filter(sharks.HMMtracks.df,
                       SharksexTrackNo %in% c("WSF1 T1", "WSF1 T10")),
  ggplot2::aes(Long, Lat)
) +
  ggplot2::geom_path(ggplot2::aes(group = SharksexTrackNo),
                     alpha = 0.5, color = "grey") +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::facet_wrap(~SharksexTrackNo) +
  ggplot2::theme_minimal() +
  ggplot2::ylab("Latitude") +
  ggplot2::xlab("Longitude") +
  ggplot2::ggtitle("Figure 26.1: White shark tracks — WSF1 T1 and WSF1 T10")

print(fig26_1)
ggplot2::ggsave("figs/Fig-26.1.svg", fig26_1, width = 8, height = 4)

# ============================================================
## §8 Helper functions
# ============================================================

# Evaluate posterior Gamma state-dependent densities over a grid
eval_gamma_sdd <- function(post_draws, no_samples) {
  xval        <- seq(0.001, 1, len = 200)
  state1_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  state2_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  for (j in seq_len(no_samples)) {
    state1_dens[, j] <- dgamma(xval,
                                shape = as.numeric(post_draws[j, "shape[1]"]),
                                rate  = as.numeric(post_draws[j, "rate[1]"]))
    state2_dens[, j] <- dgamma(xval,
                                shape = as.numeric(post_draws[j, "shape[2]"]),
                                rate  = as.numeric(post_draws[j, "rate[2]"]))
  }
  data.frame(
    xval        = rep(xval, times = no_samples),
    state1_dens = c(state1_dens),
    state2_dens = c(state2_dens),
    sample      = rep(seq_len(no_samples), each = 200)
  )
}

# Evaluate posterior von Mises state-dependent densities over a grid
eval_vonMises_sdd <- function(post_draws, no_samples) {
  xval        <- seq(-pi, pi, len = 200)
  state1_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  state2_dens <- matrix(NA, nrow = 200, ncol = no_samples)
  for (j in seq_len(no_samples)) {
    state1_dens[, j] <- CircStats::dvm(xval,
                                        mu    = as.numeric(post_draws[j, "loc[1]"]),
                                        kappa = as.numeric(post_draws[j, "kappa[1]"]))
    state2_dens[, j] <- CircStats::dvm(xval,
                                        mu    = as.numeric(post_draws[j, "loc[2]"]),
                                        kappa = as.numeric(post_draws[j, "kappa[2]"]))
  }
  data.frame(
    xval        = rep(xval, times = no_samples),
    state1_dens = c(state1_dens),
    state2_dens = c(state2_dens),
    sample      = rep(seq_len(no_samples), each = 200)
  )
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ============================================================
## §9 Model 1: Baseline 2-state HMM
# ============================================================

log_result("\n=== Model 1: baseline 2-state HMM (v2 — Pathfinder init) ===")

stanHMM_2states <- list(
  Nstates     = 2,
  Tlen        = nrow(ws_HMM),
  Ntracks     = length(unique(ws_HMM$SharksexTrackNo)),
  track_index = as.integer(as.factor(ws_HMM$SharksexTrackNo)),
  steplength  = ws_HMM$steplength,
  angle       = ws_HMM$turnang
)

log_result("Stan data dimensions: Tlen = ", stanHMM_2states$Tlen,
           ", Ntracks = ", stanHMM_2states$Ntracks)

model_2stateHMM <- cmdstanr::cmdstan_model("step_turn_hmm.stan")

# --- Pathfinder initialisation (Lesson 3 from session log) ---
# Run Pathfinder once on each chain starting point; use the resulting
# draws as init for the full NUTS sampler.  This reliably finds the
# posterior basin before NUTS warmup begins, preventing the state-swap
# trap that produced Rhat 1.73 in v1.
#
# We provide 4 Pathfinder runs (one per chain) with informative starts:
#   mu = (0.08, 0.35), sigma = (0.1, 0.15), tpm diagonal near 0.95
# This anchors each run in the correct basin while still sampling
# the uncertainty in those parameters.

make_hmm_init <- function(mu1 = 0.08, mu2 = 0.35,
                           s1  = 0.10, s2  = 0.15) {
  list(
    mu          = c(mu1, mu2),
    sigma       = c(s1,  s2),
    mixp        = c(0.02, 0.02),
    xangle      = c(-1.0,  2.0),  # state 1 ~resident (mean angle ~pi), state 2 ~directed
    yangle      = c( 0.0,  0.0),
    tpm_raw     = list(c(0.95, 0.05), c(0.05, 0.95)),
    initial_dist = c(0.5, 0.5)
  )
}

set.seed(SEED)
pf_2state <- tryCatch(
  model_2stateHMM$pathfinder(
    data      = stanHMM_2states,
    init      = make_hmm_init,   # single zero-arg function; cmdstanr calls it per path
    num_paths = 4,
    seed      = SEED,
    refresh   = 0
  ),
  error = function(e) {
    message("Pathfinder failed (", conditionMessage(e), "); falling back to fixed init.")
    NULL
  }
)

# Build per-chain inits from Pathfinder draws
if (!is.null(pf_2state)) {
  pf_draws   <- pf_2state$draws(format = "df")
  pf_n_draws <- nrow(pf_draws)

  # Sample 4 rows (one per chain) spread across the Pathfinder output
  pf_idx <- round(seq(1, pf_n_draws, length.out = 4))

  build_init_from_pf <- function(i) {
    row <- pf_draws[pf_idx[i], ]
    list(
      mu          = c(as.numeric(row[["mu[1]"]]),   as.numeric(row[["mu[2]"]])),
      sigma       = c(as.numeric(row[["sigma[1]"]]),as.numeric(row[["sigma[2]"]])),
      mixp        = c(as.numeric(row[["mixp[1]"]]), as.numeric(row[["mixp[2]"]])),
      xangle      = c(as.numeric(row[["xangle[1]"]]),as.numeric(row[["xangle[2]"]])),
      yangle      = c(as.numeric(row[["yangle[1]"]]),as.numeric(row[["yangle[2]"]])),
      tpm_raw     = list(
        c(as.numeric(row[["tpm_raw[1,1]"]]), as.numeric(row[["tpm_raw[1,2]"]])),
        c(as.numeric(row[["tpm_raw[2,1]"]]), as.numeric(row[["tpm_raw[2,2]"]]))
      ),
      initial_dist = c(as.numeric(row[["initial_dist[1]"]]),
                       as.numeric(row[["initial_dist[2]"]]))
    )
  }

  hmm_inits <- lapply(1:4, build_init_from_pf)
  log_result("Pathfinder complete. mu[1] range from PF: ",
             round(range(pf_draws[["mu[1]"]]), 3))
  log_result("Pathfinder complete. mu[2] range from PF: ",
             round(range(pf_draws[["mu[2]"]]), 3))
} else {
  log_result("Pathfinder skipped. Using fixed informative init for all 4 chains.")
  hmm_inits <- replicate(4, make_hmm_init(), simplify = FALSE)
}

# Full NUTS sampler with Pathfinder inits
set.seed(SEED)
fit_2stateHMM <- model_2stateHMM$sample(
  data          = stanHMM_2states,
  init          = hmm_inits,
  chains        = 4,
  iter_warmup   = 2000,
  iter_sampling = 2000,
  refresh       = 200,
  seed          = SEED
)

saveRDS(fit_2stateHMM, "fit_2stateHMM.rds")
log_result("fit_2stateHMM saved.")

# --- Phase 4 diagnostics ---
wf$fit_timestamp <- Sys.time()
wf$stan_backend  <- "cmdstanr"
wf$formula       <- "2-state HMM: Gamma(shape,rate) + ZeroMass + vonMises"
wf$family        <- "HMM / step-turn"

fit_2stateHMM_summ <- fit_2stateHMM$summary(
  variables = c("mu", "sigma", "mixp", "shape", "rate", "loc", "kappa",
                "tpm", "initial_dist", "lp__")
)

log_result("--- Diagnostics: fit_2stateHMM ---")
diag_df <- fit_2stateHMM$diagnostic_summary(quiet = TRUE)
log_result("  divergences (total): ",
           sum(diag_df$num_divergent, na.rm = TRUE))
log_result("  max_treedepth hits : ",
           sum(diag_df$num_max_treedepth, na.rm = TRUE))
rhat_max <- max(fit_2stateHMM_summ$rhat, na.rm = TRUE)
ess_min  <- min(fit_2stateHMM_summ$ess_bulk, na.rm = TRUE)
log_result("  Rhat_max    : ", round(rhat_max, 4))
log_result("  ESS_bulk_min: ", round(ess_min, 0))

wf$diagnostics$rhat_max     <- rhat_max
wf$diagnostics$bulk_ess_min  <- ess_min
wf$diagnostics$n_divergences <- sum(diag_df$num_divergent, na.rm = TRUE)
wf$diagnostics$passed <- (rhat_max < 1.01 &&
                           ess_min > 400 &&
                           wf$diagnostics$n_divergences == 0)
log_result("  passed      : ", wf$diagnostics$passed)

if (!wf$diagnostics$passed)
  log_result("  NOTE: if Rhat still elevated, inspect mu chains with",
             " mcmc_trace(fit_2stateHMM$draws('mu')) for residual",
             " label-switching and consider tightening mu[2] prior further.")

# --- Parameter summary ---
log_result("\n=== Model 1 parameter summary (key params) ===")
key_vars <- dplyr::filter(fit_2stateHMM_summ,
                           variable %in% c("mu[1]", "mu[2]",
                                           "shape[1]", "shape[2]",
                                           "rate[1]",  "rate[2]",
                                           "tpm[1,1]", "tpm[1,2]",
                                           "tpm[2,1]", "tpm[2,2]"))
capture_result(key_vars, label = "Key parameter posterior summaries:")

# ============================================================
## §10 Figure 26.2 — MCMC marginal histograms by chain for mu
# ============================================================
fit_2stateHMM_draws <- fit_2stateHMM$draws(
  format    = "df",
  variables = c("mu", "sigma", "mixp", "shape", "rate",
                "xangle", "yangle", "kappa", "loc", "tpm",
                "initial_dist", "lp__")
)

# <inspect mu[1] and mu[2] chains for bimodality before continuing>

fig26_2 <- bayesplot::mcmc_hist_by_chain(fit_2stateHMM_draws,
                                          regex_pars = "^mu",
                                          pars       = "lp__")
print(fig26_2)
ggplot2::ggsave("figs/Fig-26.2.svg", fig26_2, width = 8, height = 5)
log_result("Fig 26.2 saved.")

# ============================================================
## §11 Figure 26.3 — State-dependent distributions
# ============================================================
sdd_steplength <- eval_gamma_sdd(fit_2stateHMM_draws, no_samples = 1000)
sdd_angle      <- eval_vonMises_sdd(fit_2stateHMM_draws, no_samples = 1000)

sl_sdd <- ggplot2::ggplot(dplyr::filter(ws_HMM, steplength >= 0)) +
  ggplot2::geom_histogram(ggplot2::aes(steplength,
                                       y = ggplot2::after_stat(density)),
                          bins = 100, alpha = 0.2) +
  ggplot2::xlim(-0.01, 1) +
  ggplot2::geom_line(data = sdd_steplength,
                     ggplot2::aes(xval, 0.5 * state1_dens, group = sample),
                     color = cbPalette[2], alpha = 0.2) +
  ggplot2::geom_line(data = sdd_steplength,
                     ggplot2::aes(xval, 0.5 * state2_dens, group = sample),
                     color = cbPalette[3], alpha = 0.2) +
  ggplot2::geom_line(data = sdd_steplength,
                     ggplot2::aes(xval,
                                  0.5 * state1_dens + 0.5 * state2_dens,
                                  group = sample),
                     color = "darkgrey", alpha = 0.02) +
  ggplot2::theme_minimal() +
  ggplot2::xlab("") +
  ggplot2::ylab("density") +
  ggplot2::annotate("text", x = 0.12, y = 6,   label = "state 1",
                    color = cbPalette[2]) +
  ggplot2::annotate("text", x = 0.45, y = 1.5, label = "state 2",
                    color = cbPalette[3]) +
  ggplot2::ggtitle("step length state-dependent distributions")

angle_sdd <- ggplot2::ggplot(dplyr::filter(ws_HMM, turnang >= -pi)) +
  ggplot2::geom_histogram(ggplot2::aes(turnang,
                                       y = ggplot2::after_stat(density)),
                          bins = 100, alpha = 0.2) +
  ggplot2::xlim(-pi, pi) +
  ggplot2::geom_line(data = sdd_angle,
                     ggplot2::aes(xval, 0.5 * state1_dens, group = sample),
                     color = cbPalette[2], alpha = 0.2) +
  ggplot2::geom_line(data = sdd_angle,
                     ggplot2::aes(xval, 0.5 * state2_dens, group = sample),
                     color = cbPalette[3], alpha = 0.2) +
  ggplot2::geom_line(data = sdd_angle,
                     ggplot2::aes(xval,
                                  0.5 * state1_dens + 0.5 * state2_dens,
                                  group = sample),
                     color = "darkgrey", alpha = 0.02) +
  ggplot2::theme_minimal() +
  ggplot2::xlab("") +
  ggplot2::ylab("density") +
  ggplot2::annotate("text", x = -2.5, y = 0.25, label = "state 1",
                    color = cbPalette[2]) +
  ggplot2::annotate("text", x = 1.5,  y = 0.25, label = "state 2",
                    color = cbPalette[3]) +
  ggplot2::ggtitle("turning angle state-dependent distributions")

fig26_3 <- sl_sdd + angle_sdd
print(fig26_3)
ggplot2::ggsave("figs/Fig-26.3.svg", fig26_3, width = 10, height = 4)
log_result("Fig 26.3 saved.")

# ============================================================
## §12 Figures 26.4 & 26.5 — Local state decoding via forward-backward
# ============================================================
Tlen_actual <- nrow(ws_HMM)
log_result("Tlen_actual = ", Tlen_actual)

state_probs_draws <- fit_2stateHMM$draws(
  variables = c("state_probs"),
  format    = "draws_matrix"
)

# state_probs[1,1..T] then state_probs[2,1..T] in draws_matrix columns
half <- Tlen_actual
state1_cols <- seq_len(half)
state2_cols <- half + seq_len(half)

state_probs_means <- data.frame(
  state1prob = colMeans(state_probs_draws[seq_len(1000), state2_cols]),
  state2prob = colMeans(state_probs_draws[seq_len(1000), state1_cols])
)
state1_probs_quants <- data.frame(
  state1prob025 = apply(state_probs_draws[seq_len(1000), state2_cols],
                        2, quantile, probs = 0.025),
  state1prob975 = apply(state_probs_draws[seq_len(1000), state2_cols],
                        2, quantile, probs = 0.975)
)
state2_probs_quants <- data.frame(
  state2prob025 = apply(state_probs_draws[seq_len(1000), state1_cols],
                        2, quantile, probs = 0.025),
  state2prob975 = apply(state_probs_draws[seq_len(1000), state1_cols],
                        2, quantile, probs = 0.975)
)

ws_HMM_rep <- ws_HMM_full[, c("dateTime", "SharkName", "SharksexTrackNo",
                                "steplength", "turnang", "year", "month",
                                "CDB", "Lat", "Long")]
ws_HMM_rep$state1prob    <- state_probs_means[, 1]
ws_HMM_rep$state2prob    <- state_probs_means[, 2]
ws_HMM_rep$state1prob025 <- state1_probs_quants[, 1]
ws_HMM_rep$state1prob975 <- state1_probs_quants[, 2]
ws_HMM_rep$state2prob025 <- state2_probs_quants[, 1]
ws_HMM_rep$state2prob975 <- state2_probs_quants[, 2]

fig26_4 <- ggplot2::ggplot(
  data = dplyr::filter(ws_HMM_rep,
                       SharksexTrackNo %in% c("WSF1 T1", "WSF1 T10")),
  ggplot2::aes(Long, Lat)
) +
  ggplot2::geom_path(ggplot2::aes(group = SharksexTrackNo),
                     alpha = 0.5, color = "grey") +
  ggplot2::geom_point(ggplot2::aes(color = state1prob)) +
  ggplot2::labs(color = "State 1\nProbability") +
  ggplot2::scale_color_viridis_c(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  ggplot2::facet_wrap(~SharksexTrackNo) +
  ggplot2::theme_minimal() +
  ggplot2::ylab("Latitude") +
  ggplot2::xlab("Longitude") +
  ggplot2::ggtitle("Figure 26.4: State-1 probability on WSF1 tracks")

print(fig26_4)
ggplot2::ggsave("figs/Fig-26.4.svg", fig26_4, width = 9, height = 4)
log_result("Fig 26.4 saved.")

fig26_5 <- ggplot2::ggplot(
  data = dplyr::filter(ws_HMM_rep,
                       SharksexTrackNo %in% c("WSF1 T1", "WSF1 T10")),
  ggplot2::aes(dateTime, state1prob)
) +
  ggplot2::theme_minimal() +
  ggplot2::geom_point() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = state1prob025,
                                    ymax = state1prob975),
                       alpha = 0.5) +
  ggplot2::facet_wrap(~SharksexTrackNo, nrow = 2, scales = "free_x") +
  ggplot2::ylab("State 1 Probability") +
  ggplot2::xlab("Time") +
  ggplot2::ggtitle("Figure 26.5: State-1 probability time series with 95% CI")

print(fig26_5)
ggplot2::ggsave("figs/Fig-26.5.svg", fig26_5, width = 9, height = 6)
log_result("Fig 26.5 saved.")

# ============================================================
## §13 Figure 26.6 — Pseudo-residual Q-Q plot
# ============================================================
pseudo_residuals <- fit_2stateHMM$draws(
  variables = "pseudo_residuals",
  format    = "draws_matrix"
)

pr_mean          <- colMeans(pseudo_residuals[seq_len(1000), ])
pr_missing_index <- which(pr_mean == 0)
df_resid         <- data.frame(y = pr_mean[-pr_missing_index])

log_result("Pseudo-residuals: n non-missing = ", nrow(df_resid))

fig26_6 <- ggplot2::ggplot(df_resid, ggplot2::aes(sample = y)) +
  ggplot2::stat_qq() +
  ggplot2::stat_qq_line() +
  ggplot2::theme_minimal() +
  ggplot2::xlab("Theoretical Quantiles") +
  ggplot2::ylab("Sample Quantiles") +
  ggplot2::ggtitle("Figure 26.6: Pseudo-residual Q-Q plot (step lengths)")

print(fig26_6)
ggplot2::ggsave("figs/Fig-26.6.svg", fig26_6, width = 5, height = 5)
log_result("Fig 26.6 saved.")

# ============================================================
## §14 Figure 26.7 — State-decoding histograms
# ============================================================
post_state_samples  <- fit_2stateHMM$draws("state_sequence", format = "draws_array")
init_dist_draws     <- fit_2stateHMM$draws("initial_dist",   format = "draws_array")
tpm_draws           <- fit_2stateHMM$draws("tpm",            format = "draws_array")

n_ppd_sims  <- 1000
state_samples <- matrix(NA_integer_, nrow = Tlen_actual, ncol = n_ppd_sims)

set.seed(SEED)
for (j in seq_len(n_ppd_sims)) {
  tpm_iter <- matrix(
    data = c(tpm_draws[j, 1, 1], tpm_draws[j, 1, 2],
             tpm_draws[j, 1, 3], tpm_draws[j, 1, 4]),
    nrow = 2, byrow = TRUE
  )
  init_p <- as.numeric(init_dist_draws[j, 1, ])
  state_samples[1, j] <- sample(2:1, size = 1, prob = init_p)
  for (t in 2:Tlen_actual) {
    if (ws_HMM$SharksexTrackNo[t] == ws_HMM$SharksexTrackNo[t - 1]) {
      state_samples[t, j] <- sample(
        2:1, size = 1,
        prob = tpm_iter[state_samples[t - 1, j], ]
      )
    } else {
      state_samples[t, j] <- sample(2:1, size = 1, prob = init_p)
    }
  }
}

track1_idx <- 1:91

sim_state_counts          <- apply(state_samples[track1_idx, ], 2, table)
sim_state_props           <- sim_state_counts / length(track1_idx)
post_state_samples_counts <- apply(post_state_samples[, 1, track1_idx], 1, table)
post_state_samples_props  <- post_state_samples_counts / length(track1_idx)

ws_HMM_rep$chum <- ifelse(ws_HMM_full$CDB == "x", 1L, 0L)
ws_HMM_rep$chum[which(is.na(ws_HMM_rep$CDB))] <- 0L

chum_ws1_tr1 <- which(ws_HMM_rep$chum[track1_idx] == 1)

post_state_samples_chumcounts <- apply(
  post_state_samples[, 1, chum_ws1_tr1], 1, table
)
n_chum <- length(chum_ws1_tr1)
post_state_samples_st2chumprops <- vapply(
  seq_len(n_ppd_sims),
  function(j) {
    ct <- post_state_samples_chumcounts[[j]]
    if (length(ct) == 0) return(0)
    max(ct) / n_chum
  },
  numeric(1L)
)

log_result("Chum observations in WSF1 T1: ", n_chum)

fig26_7 <- ggplot2::ggplot(
  data = data.frame(y = sim_state_props[1, ]),
  ggplot2::aes(y)
) +
  ggplot2::geom_histogram(bins = 100, fill = "darkgrey") +
  ggplot2::theme_minimal() +
  ggplot2::xlim(-0.05, 1) +
  ggplot2::annotate("text", label = "Posterior\nPredictive\nSimulations",
                    x = 0.63, y = 200, col = "darkgrey") +
  ggplot2::geom_histogram(
    data = data.frame(x = post_state_samples_props[1, ]),
    ggplot2::aes(x), bins = 100, fill = "black"
  ) +
  ggplot2::annotate("text", label = "State\nDecodings",
                    x = 0.35, y = 200, col = "black") +
  ggplot2::geom_histogram(
    data = data.frame(p = 1 - post_state_samples_st2chumprops),
    ggplot2::aes(p), bins = 100, fill = "blue"
  ) +
  ggplot2::annotate("text", label = "Chum",
                    x = 0.07, y = 200, col = "blue") +
  ggplot2::xlab("Proportion of State 1 Observations") +
  ggplot2::ylab("Count") +
  ggplot2::ggtitle("Figure 26.7: State decoding histograms")

print(fig26_7)
ggplot2::ggsave("figs/Fig-26.7.svg", fig26_7, width = 7, height = 5)
log_result("Fig 26.7 saved.")

# ============================================================
## §15 Model 2: 2-state HMM with covariates in tpm
# ============================================================
log_result("\n=== Model 2: 2-state HMM with covariates in tpm ===")

ws_HMM$tod_cos  <- cos((2 * pi * (lubridate::hour(ws_HMM$dateTime) * 60 +
                                    lubridate::minute(ws_HMM$dateTime))) / 1440)
ws_HMM$tod_sin  <- sin((2 * pi * (lubridate::hour(ws_HMM$dateTime) * 60 +
                                    lubridate::minute(ws_HMM$dateTime))) / 1440)
ws_HMM$chum     <- ifelse(ws_HMM$CDB == "x", 1L, 0L)
ws_HMM$chum[which(is.na(ws_HMM$chum))] <- 0L
ws_HMM$sex_char <- substring(ws_HMM$SharksexTrackNo, 3, 3)
ws_HMM$sex      <- ifelse(ws_HMM$sex_char == "F", 0L, 1L)

HMM_covar <- cbind(
  1,
  ws_HMM$chum,
  ws_HMM$sex,
  ws_HMM$tod_cos,
  ws_HMM$tod_sin
)
log_result("HMM_covar dim: ", paste(dim(HMM_covar), collapse = " x "))

stanHMM_2states_covariates <- list(
  Nstates     = 2,
  Tlen        = nrow(ws_HMM),
  Ntracks     = length(unique(ws_HMM$SharksexTrackNo)),
  track_index = as.integer(as.factor(ws_HMM$SharksexTrackNo)),
  steplength  = ws_HMM$steplength,
  angle       = ws_HMM$turnang,
  nCovs       = 4,
  covs        = HMM_covar
)

model_2stateHMM_covariates <- cmdstanr::cmdstan_model(
  "step_turn_hmm_covariates.stan"
)

# Pathfinder init for Model 2
# The covariates model has more parameters (beta matrix) but the
# state-dependent SDD parameters have the same label-switching risk.
# Seed the Pathfinder with the same informative starting point.

make_cov_init <- function(mu1 = 0.08, mu2 = 0.35) {
  list(
    mu          = c(mu1, mu2),
    sigma       = c(0.10, 0.15),
    mixp        = c(0.02, 0.02),
    xangle      = c(-1.0,  2.0),
    yangle      = c( 0.0,  0.0),
    beta        = matrix(c(-2, -2, 0, 0, 0,
                            -2, -2, 0, 0, 0),
                         nrow = 2, ncol = 5, byrow = TRUE),
    initial_dist = c(0.5, 0.5)
  )
}

set.seed(SEED)
pf_cov <- tryCatch(
  model_2stateHMM_covariates$pathfinder(
    data      = stanHMM_2states_covariates,
    init      = make_cov_init,   # zero-arg function; cmdstanr calls it per path
    num_paths = 4,
    seed      = SEED,
    refresh   = 0
  ),
  error = function(e) {
    message("Pathfinder (Model 2) failed (", conditionMessage(e), "); falling back to fixed init.")
    NULL
  }
)

if (!is.null(pf_cov)) {
  pf_cov_draws <- pf_cov$draws(format = "df")
  n_pf_cov     <- nrow(pf_cov_draws)
  pf_cov_idx   <- round(seq(1, n_pf_cov, length.out = 4))

  extract_cov_init <- function(i) {
    row <- pf_cov_draws[pf_cov_idx[i], ]
    # Pull all beta[r,c] entries from the Pathfinder draw
    beta_mat <- matrix(NA_real_, nrow = 2, ncol = 5)
    for (r in 1:2)
      for (cc in 1:5)
        beta_mat[r, cc] <- as.numeric(row[[paste0("beta[", r, ",", cc, "]")]])
    list(
      mu          = c(as.numeric(row[["mu[1]"]]),    as.numeric(row[["mu[2]"]])),
      sigma       = c(as.numeric(row[["sigma[1]"]]), as.numeric(row[["sigma[2]"]])),
      mixp        = c(as.numeric(row[["mixp[1]"]]),  as.numeric(row[["mixp[2]"]])),
      xangle      = c(as.numeric(row[["xangle[1]"]]),as.numeric(row[["xangle[2]"]])),
      yangle      = c(as.numeric(row[["yangle[1]"]]),as.numeric(row[["yangle[2]"]])),
      beta        = beta_mat,
      initial_dist = c(as.numeric(row[["initial_dist[1]"]]),
                       as.numeric(row[["initial_dist[2]"]]))
    )
  }

  cov_inits <- lapply(1:4, extract_cov_init)
  log_result("Pathfinder (Model 2) complete.")
} else {
  log_result("Pathfinder (Model 2) skipped. Using fixed informative init for all 4 chains.")
  cov_inits <- replicate(4, make_cov_init(), simplify = FALSE)
}

set.seed(SEED)
fit_2stateHMM_covariates <- model_2stateHMM_covariates$sample(
  data          = stanHMM_2states_covariates,
  init          = cov_inits,
  chains        = 4,
  iter_warmup   = 2000,
  iter_sampling = 2000,
  refresh       = 200,
  seed          = SEED
)

saveRDS(fit_2stateHMM_covariates, "fit_2stateHMM_covariates.rds")
log_result("fit_2stateHMM_covariates saved.")

fit_2stateHMM_cov_summ <- fit_2stateHMM_covariates$summary(
  variables = c("mu", "shape", "rate", "lp__")
)
diag_df2 <- fit_2stateHMM_covariates$diagnostic_summary(quiet = TRUE)
log_result("--- Diagnostics: fit_2stateHMM_covariates ---")
log_result("  divergences : ", sum(diag_df2$num_divergent, na.rm = TRUE))
log_result("  Rhat_max    : ",
           round(max(fit_2stateHMM_cov_summ$rhat, na.rm = TRUE), 4))
log_result("  ESS_bulk_min: ",
           round(min(fit_2stateHMM_cov_summ$ess_bulk, na.rm = TRUE), 0))

# ============================================================
## §16 Model 3: Covariates + non-centred individual random effects
# ============================================================
log_result("\n=== Model 3: Covariates + NCP random effects ===")

stanHMM_2states_tpmcov_crencp <- list(
  Nstates     = 2,
  Tlen        = nrow(ws_HMM),
  Ntracks     = length(unique(ws_HMM$SharksexTrackNo)),
  track_index = as.integer(as.factor(ws_HMM$SharksexTrackNo)),
  steplength  = ws_HMM$steplength,
  angle       = ws_HMM$turnang,
  nCovs       = 4,
  covs        = HMM_covar[, -1],
  Nsharks     = length(unique(ws_HMM$SharksexTrackNo)),
  shark_index = as.integer(as.factor(ws_HMM$SharksexTrackNo))
)

model_2stateHMM_tpmcov_crencp <- cmdstanr::cmdstan_model(
  "step_turn_hmm_covariates_cre_ncp.stan"
)

# Model 3: 1 chain as in book; NCP parameterisation converges
# well without Pathfinder (no label-switching in sigma_re / z_re)
set.seed(SEED)
fit_2stateHMM_tpmcov_crencp <- model_2stateHMM_tpmcov_crencp$sample(
  data          = stanHMM_2states_tpmcov_crencp,
  init          = list(make_hmm_init()),  # reuse SDD init; beta/RE default to 0
  chains        = 1,
  iter_warmup   = 2000,
  iter_sampling = 2000,
  refresh       = 200,
  seed          = SEED
)

saveRDS(fit_2stateHMM_tpmcov_crencp, "fit_2stateHMM_tpmcov_crencp.rds")
log_result("fit_2stateHMM_tpmcov_crencp saved (1 chain).")

diag_df3 <- fit_2stateHMM_tpmcov_crencp$diagnostic_summary(quiet = TRUE)
fit3_summ <- fit_2stateHMM_tpmcov_crencp$summary(
  variables = c("mu", "mu_tpm", "sigma_re", "lp__")
)
log_result("--- Diagnostics: fit_2stateHMM_tpmcov_crencp ---")
log_result("  divergences : ", sum(diag_df3$num_divergent, na.rm = TRUE))
log_result("  Rhat_max    : ",
           round(max(fit3_summ$rhat, na.rm = TRUE), 4))
log_result("  ESS_bulk_min: ",
           round(min(fit3_summ$ess_bulk, na.rm = TRUE), 0))

# ============================================================
## §17 Figures 26.8 & 26.9 — tpm entries vs time (base-R, model 3)
# ============================================================
beta    <- fit_2stateHMM_tpmcov_crencp$draws("beta",        format = "draws_array")
randeff <- fit_2stateHMM_tpmcov_crencp$draws("randeff_tpm", format = "draws_array")
mu_tpm  <- fit_2stateHMM_tpmcov_crencp$draws("mu_tpm",      format = "draws_array")

grid_tod <- 480:1200

DM_male_nochum <- cbind(
  1,
  rep(0, length(grid_tod)),
  rep(1, length(grid_tod)),
  cos(2 * pi * grid_tod / 1440),
  sin(2 * pi * grid_tod / 1440)
)

DM_female_nochum <- cbind(
  1,
  rep(0, length(grid_tod)),
  rep(0, length(grid_tod)),
  cos(2 * pi * grid_tod / 1440),
  sin(2 * pi * grid_tod / 1440)
)

build_tpm_grid <- function(beta_mat_rbind, DM) {
  moveHMM:::trMatrix_rcpp(nbStates = 2, beta = t(beta_mat_rbind), covs = DM)
}

# Fig 26.8 — base-R: male, no chum
dev.new()
par(mfrow = c(2, 2))
for (i in 1:2) {
  for (j in 1:2) {
    bm <- rbind(c(mu_tpm[1, 1, 1], beta[1, 1, c(1, 3, 5, 7)]),
                c(mu_tpm[1, 1, 2], beta[1, 1, c(2, 4, 6, 8)]))
    tpm_base <- build_tpm_grid(bm, DM_male_nochum)
    plot(grid_tod, tpm_base[i, j, ], type = "l", ylim = c(0, 1),
         col = "grey", lwd = 0.5,
         xlab = "minute of the day",
         ylab = paste0("Pr(", i, " -> ", j, ")"))
    for (k in 2:100) {
      bm <- rbind(c(mu_tpm[k, 1, 1], beta[k, 1, c(1, 3, 5, 7)]),
                  c(mu_tpm[k, 1, 2], beta[k, 1, c(2, 4, 6, 8)]))
      tpm_k <- build_tpm_grid(bm, DM_male_nochum)
      lines(grid_tod, tpm_k[i, j, ], col = "grey", lwd = 0.1)
    }
  }
}
dev.copy(svg, "figs/Fig-26.8.svg", width = 7, height = 7)
dev.off()
par(mfrow = c(1, 1))
log_result("Fig 26.8 saved (base-R male/no-chum tpm grid).")

# Fig 26.9 — base-R: female, no chum
dev.new()
par(mfrow = c(2, 2))
for (i in 1:2) {
  for (j in 1:2) {
    bm <- rbind(c(mu_tpm[1, 1, 1], beta[1, 1, c(1, 3, 5, 7)]),
                c(mu_tpm[1, 1, 2], beta[1, 1, c(2, 4, 6, 8)]))
    tpm_base <- build_tpm_grid(bm, DM_female_nochum)
    plot(grid_tod, tpm_base[i, j, ], type = "l", ylim = c(0, 1),
         col = "grey", lwd = 0.5,
         xlab = "minute of the day",
         ylab = paste0("Pr(", i, " -> ", j, ")"))
    for (k in 2:100) {
      bm <- rbind(c(mu_tpm[k, 1, 1], beta[k, 1, c(1, 3, 5, 7)]),
                  c(mu_tpm[k, 1, 2], beta[k, 1, c(2, 4, 6, 8)]))
      tpm_k <- build_tpm_grid(bm, DM_female_nochum)
      lines(grid_tod, tpm_k[i, j, ], col = "grey", lwd = 0.1)
    }
  }
}
dev.copy(svg, "figs/Fig-26.9.svg", width = 7, height = 7)
dev.off()
par(mfrow = c(1, 1))
log_result("Fig 26.9 saved (base-R female/no-chum tpm grid).")

# ============================================================
## §18 Figures 26.10-26.12 — ggplot tpm ribbons
# ============================================================

build_tpm_df <- function(mu_tpm, beta, DM, n_draws = 500) {
  bm    <- rbind(c(mu_tpm[1, 1, 1], beta[1, 1, c(1, 3, 5, 7)]),
                 c(mu_tpm[1, 1, 2], beta[1, 1, c(2, 4, 6, 8)]))
  tpm_0 <- build_tpm_grid(bm, DM)
  df <- data.frame(
    omega11 = tpm_0[1, 1, ], omega12 = tpm_0[1, 2, ],
    omega21 = tpm_0[2, 1, ], omega22 = tpm_0[2, 2, ],
    dateTime = grid_tod, draw = 1L
  )
  for (k in 2:n_draws) {
    bm    <- rbind(c(mu_tpm[k, 1, 1], beta[k, 1, c(1, 3, 5, 7)]),
                   c(mu_tpm[k, 1, 2], beta[k, 1, c(2, 4, 6, 8)]))
    tpm_k <- build_tpm_grid(bm, DM)
    df <- rbind(df, data.frame(
      omega11 = tpm_k[1, 1, ], omega12 = tpm_k[1, 2, ],
      omega21 = tpm_k[2, 1, ], omega22 = tpm_k[2, 2, ],
      dateTime = grid_tod, draw = k
    ))
  }
  colnames(df)[1:4] <- c("omega[11](t)", "omega[12](t)",
                          "omega[21](t)", "omega[22](t)")
  tidyr::pivot_longer(df, !c(dateTime, draw), names_to = "omega")
}

tpm_gg <- function(long_df, title_str) {
  ggplot2::ggplot(long_df, ggplot2::aes(dateTime, value)) +
    ggplot2::geom_line(ggplot2::aes(group = draw),
                       alpha = 0.1, col = "darkgrey") +
    ggplot2::facet_wrap(~omega, nrow = 2, ncol = 2,
                        labeller = ggplot2::label_parsed) +
    ggplot2::theme_minimal() +
    ggplot2::ylab("Probability") +
    ggplot2::xlab("Time") +
    ggplot2::scale_x_continuous(
      breaks = c(540, 720, 900, 1080),
      labels = c("09:00", "12:00", "15:00", "18:00")
    ) +
    ggplot2::ggtitle(title_str)
}

p1 <- tpm_gg(build_tpm_df(mu_tpm, beta, DM_male_nochum),
              "Figure 26.10: Male white shark, no chum")
print(p1)
ggplot2::ggsave("figs/Fig-26.10.svg", p1, width = 8, height = 6)
log_result("Fig 26.10 saved.")

p2 <- tpm_gg(build_tpm_df(mu_tpm, beta, DM_female_nochum),
              "Figure 26.11: Female white shark, no chum")
print(p2)
ggplot2::ggsave("figs/Fig-26.11.svg", p2, width = 8, height = 6)
log_result("Fig 26.11 saved.")

# WSF1 T1 individual (uses randeff_tpm for that shark)
wsf1_t1_idx   <- which(ws_HMM$SharksexTrackNo == "WSF1 T1")
wsf1_t1_covar <- cbind(1, HMM_covar[wsf1_t1_idx, -1])
wsf1_dateTime <- 60 * lubridate::hour(ws_HMM$dateTime[wsf1_t1_idx]) +
                 lubridate::minute(ws_HMM$dateTime[wsf1_t1_idx])

wsf1_t1_df <- {
  bm    <- rbind(c(randeff[1, 1, 1], beta[1, 1, c(1, 3, 5, 7)]),
                 c(randeff[1, 1, 2], beta[1, 1, c(2, 4, 6, 8)]))
  tpm_0 <- moveHMM:::trMatrix_rcpp(nbStates = 2, beta = t(bm),
                                    covs = wsf1_t1_covar)
  df <- data.frame(
    omega11 = tpm_0[1, 1, ], omega12 = tpm_0[1, 2, ],
    omega21 = tpm_0[2, 1, ], omega22 = tpm_0[2, 2, ],
    dateTime = wsf1_dateTime, draw = 1L
  )
  for (k in 2:500) {
    bm    <- rbind(c(mu_tpm[k, 1, 1], beta[k, 1, c(1, 3, 5, 7)]),
                   c(mu_tpm[k, 1, 2], beta[k, 1, c(2, 4, 6, 8)]))
    tpm_k <- moveHMM:::trMatrix_rcpp(nbStates = 2, beta = t(bm),
                                      covs = wsf1_t1_covar)
    df <- rbind(df, data.frame(
      omega11 = tpm_k[1, 1, ], omega12 = tpm_k[1, 2, ],
      omega21 = tpm_k[2, 1, ], omega22 = tpm_k[2, 2, ],
      dateTime = wsf1_dateTime, draw = k
    ))
  }
  colnames(df)[1:4] <- c("omega[11](t)", "omega[12](t)",
                          "omega[21](t)", "omega[22](t)")
  tidyr::pivot_longer(df, !c(dateTime, draw), names_to = "omega")
}

wsf1_chum_times <- wsf1_dateTime[which(HMM_covar[wsf1_t1_idx, 2] == 1)]

p3 <- ggplot2::ggplot(wsf1_t1_df, ggplot2::aes(dateTime, value)) +
  ggplot2::geom_line(ggplot2::aes(group = draw),
                     alpha = 0.1, col = "darkgrey") +
  ggplot2::geom_vline(xintercept = wsf1_chum_times,
                      linetype = 1, col = "lightgrey", alpha = 0.15) +
  ggplot2::facet_wrap(~omega, nrow = 2, ncol = 2,
                      labeller = ggplot2::label_parsed) +
  ggplot2::theme_minimal() +
  ggplot2::ylab("Probability") +
  ggplot2::xlab("Time") +
  ggplot2::scale_x_continuous(
    breaks = c(540, 720, 900, 1080),
    labels = c("09:00", "12:00", "15:00", "18:00")
  ) +
  ggplot2::ggtitle("Figure 26.12: White shark female 1, track 1")

print(p3)
ggplot2::ggsave("figs/Fig-26.12.svg", p3, width = 8, height = 6)
log_result("Fig 26.12 saved.")

fig26_combined <- (p1 + p2) / p3
print(fig26_combined)
ggplot2::ggsave("figs/Fig-26.combined.svg", fig26_combined,
                width = 10, height = 12)
log_result("Fig 26.combined saved.")

# ============================================================
## §19 LOO — not applicable
# ============================================================
log_result("\n=== LOO ===")
log_result("Reason: LOO is not computed for these HMM models. The latent-state")
log_result("marginalisation in the forward algorithm does not yield a pointwise")
log_result("log-likelihood decomposition compatible with PSIS-LOO without a")
log_result("leave-one-observation-out re-marginalisation, which is O(T^2) in")
log_result("cost. Model comparison in this chapter relies on visual and")
log_result("pseudo-residual diagnostics rather than LOO.")

wf$loo_complete <- FALSE

# ============================================================
## §20 Export wf_state and wrap-up
# ============================================================
export_context(wf)
log_result("\nwf_context.json written.")

log_result("\n=== Ch 26 analysis complete — ", format(Sys.time()), " ===")
log_result("Saved: fit_2stateHMM.rds, fit_2stateHMM_covariates.rds,")
log_result("       fit_2stateHMM_tpmcov_crencp.rds")
log_result("Figures: ",
  paste(
    paste0("Fig-26.", c(1:12, "combined"), ".svg"),
    collapse = ", "
  )
)

close(results_con)
