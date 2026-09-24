// step_turn_hmm.stan
//
// 2-state HMM for animal movement data (white sharks, Gansbaii, South Africa).
// State-dependent distributions:
//   step lengths : zero-inflated Gamma  (mixp = prob of zero step)
//   turning angles: von Mises
//
// Multiple independent tracks handled via track_index; missing / out-of-range
// observations are encoded as -100 and skipped in the likelihood.
//
// Generated quantities:
//   state_probs[n, t]  : local state probabilities (forward-backward)
//   state_sequence[t]  : one FFBS sample of the latent state path
//   pseudo_residuals[t]: forecast pseudo-residuals for step lengths
//
// v2: tighter priors on mu (separated by state) and sigma to prevent
//     label-switching across chains (Rhat 1.73 / ESS 6 in v1).
//
// Source: Ch 26 of Gelman, Vehtari et al. (2026) Bayesian Workflow.
//         Adapted from Leos-Barajas & Michelot (2018), arXiv:1806.10639.

data {
  int<lower=1> Nstates;
  int<lower=1> Tlen;
  int<lower=1> Ntracks;
  array[Tlen] int<lower=1, upper=Ntracks> track_index;
  vector[Tlen] steplength;   // -100 = missing/outlier
  vector[Tlen] angle;        // -100 = missing
}

parameters {
  positive_ordered[Nstates] mu;           // Gamma mean (ordered: no label-swap)
  vector<lower=0>[Nstates] sigma;         // Gamma SD
  vector<lower=0, upper=1>[Nstates] mixp; // zero-mass probability
  vector[Nstates] xangle;                 // Cartesian von Mises x-component
  vector[Nstates] yangle;                 // Cartesian von Mises y-component
  array[Nstates] simplex[Nstates] tpm_raw; // rows of the tpm
  simplex[Nstates] initial_dist;
}

transformed parameters {
  vector<lower=0>[Nstates] shape;
  vector<lower=0>[Nstates] rate;
  vector<lower=-pi(), upper=pi()>[Nstates] loc;
  vector<lower=0>[Nstates] kappa;
  matrix[Nstates, Nstates] tpm;

  for (n in 1:Nstates) {
    shape[n] = mu[n] * mu[n] / (sigma[n] * sigma[n]);
    rate[n]  = mu[n]       / (sigma[n] * sigma[n]);
    loc[n]   = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]^2 + yangle[n]^2);
  }
  for (i in 1:Nstates)
    for (j in 1:Nstates)
      tpm[i, j] = tpm_raw[i][j];
}

model {
  // --- Priors ---
  // State-specific priors on mu keep the two states well-separated
  // and prevent chains from swapping roles across warmup.
  // State 1 (area-restricted search / slow): mean ~0.08 km
  // State 2 (directed travel / fast):        mean ~0.35 km
  // These values are based on visual inspection of the step-length histogram
  // and are consistent with the moveHMM literature for white shark data.
  mu[1] ~ normal(0.08, 0.05);
  mu[2] ~ normal(0.35, 0.12);

  // Tighter sigma prior prevents overlap-induced label-switching.
  // exponential(5) has mean 0.2, 95% mass below 0.6 — reasonable for
  // step-length SD in km at this spatial scale.
  sigma ~ exponential(5.0);

  mixp   ~ beta(1, 8);          // most steps are non-zero
  xangle ~ normal(0, 2);
  yangle ~ normal(0, 1);
  for (i in 1:Nstates)
    tpm_raw[i] ~ dirichlet(rep_vector(2.0, Nstates));
  initial_dist ~ dirichlet(rep_vector(1.0, Nstates));

  // --- Log-transition matrix (transposed: log_tpm_tr[to, from]) ---
  matrix[Nstates, Nstates] log_tpm_tr;
  for (i in 1:Nstates)
    for (j in 1:Nstates)
      log_tpm_tr[j, i] = log(tpm[i, j]);

  // --- Forward algorithm (log-scale) ---
  vector[Nstates] lp;
  vector[Nstates] lp_p1;

  for (t in 1:Tlen) {
    if (t == 1 || track_index[t] != track_index[t - 1]) {
      for (n in 1:Nstates)
        lp[n] = log(initial_dist[n]);
    }
    for (n in 1:Nstates) {
      lp_p1[n] = log_sum_exp(to_row_vector(lp) + log_tpm_tr[n, 1:Nstates]);
      if (steplength[t] >= 0) {
        if (steplength[t] == 0)
          lp_p1[n] += log(mixp[n]);
        else
          lp_p1[n] += log1m(mixp[n]) + gamma_lpdf(steplength[t] | shape[n], rate[n]);
      }
      if (angle[t] >= -pi())
        lp_p1[n] += von_mises_lpdf(angle[t] | loc[n], kappa[n]);
    }
    lp = lp_p1;
    if (t == Tlen || track_index[t + 1] != track_index[t])
      target += log_sum_exp(lp);
  }
}

generated quantities {
  array[Nstates, Tlen] real state_probs;
  array[Tlen] int state_sequence;
  array[Tlen] real pseudo_residuals;

  {
    matrix[Nstates, Nstates] log_tpm_tr;
    for (i in 1:Nstates)
      for (j in 1:Nstates)
        log_tpm_tr[j, i] = log(tpm[i, j]);

    // Forward pass
    array[Nstates, Tlen] real logalpha;
    vector[Nstates] lp;
    vector[Nstates] lp_p1;

    for (t in 1:Tlen) {
      if (t == 1 || track_index[t] != track_index[t - 1]) {
        for (n in 1:Nstates)
          lp[n] = log(initial_dist[n]);
      }
      for (n in 1:Nstates) {
        lp_p1[n] = log_sum_exp(to_row_vector(lp) + log_tpm_tr[n, 1:Nstates]);
        if (steplength[t] >= 0) {
          if (steplength[t] == 0)
            lp_p1[n] += log(mixp[n]);
          else
            lp_p1[n] += log1m(mixp[n]) + gamma_lpdf(steplength[t] | shape[n], rate[n]);
        }
        if (angle[t] >= -pi())
          lp_p1[n] += von_mises_lpdf(angle[t] | loc[n], kappa[n]);
        logalpha[n, t] = lp_p1[n];
      }
      lp = lp_p1;
    }

    // Backward pass
    array[Nstates, Tlen] real logbeta;
    for (t0 in 1:Tlen) {
      int tt = Tlen - t0 + 1;
      if (tt == Tlen || track_index[tt + 1] != track_index[tt]) {
        for (n in 1:Nstates) lp_p1[n] = 0;
      } else {
        for (n in 1:Nstates) {
          lp_p1[n] = log_sum_exp(to_row_vector(lp) + log_tpm_tr[n, 1:Nstates]);
          if (steplength[tt + 1] >= 0) {
            if (steplength[tt + 1] == 0)
              lp_p1[n] += log(mixp[n]);
            else
              lp_p1[n] += log1m(mixp[n]) + gamma_lpdf(steplength[tt + 1] | shape[n], rate[n]);
          }
          if (angle[tt + 1] >= -pi())
            lp_p1[n] += von_mises_lpdf(angle[tt + 1] | loc[n], kappa[n]);
        }
      }
      lp = lp_p1;
      for (n in 1:Nstates) logbeta[n, tt] = lp[n];
    }

    // State probabilities
    for (t in 1:Tlen) {
      real llk = log_sum_exp(to_vector(logalpha[, t]) + to_vector(logbeta[, t]));
      for (n in 1:Nstates)
        state_probs[n, t] = exp(logalpha[n, t] + logbeta[n, t] - llk);
    }

    // FFBS state-sequence sample (backward sampling)
    for (t0 in 1:Tlen) {
      int tt = Tlen - t0 + 1;
      if (tt == Tlen || track_index[tt + 1] != track_index[tt]) {
        real llk_end = log_sum_exp(to_vector(logalpha[, tt]));
        vector[Nstates] pr;
        for (n in 1:Nstates) pr[n] = exp(logalpha[n, tt] - llk_end);
        state_sequence[tt] = categorical_rng(pr);
      } else {
        int ns = state_sequence[tt + 1];
        vector[Nstates] pr;
        for (n in 1:Nstates) pr[n] = exp(logalpha[n, tt]) * tpm[n, ns];
        real s = sum(pr);
        pr = (s > 0) ? pr / s : rep_vector(1.0 / Nstates, Nstates);
        state_sequence[tt] = categorical_rng(pr);
      }
    }

    // Pseudo-residuals (step lengths, forecast/one-step-ahead)
    for (t in 1:Tlen) {
      if (steplength[t] < 0) {
        pseudo_residuals[t] = 0;  // missing; flagged by caller
      } else {
        vector[Nstates] alpha_pred;
        if (t == 1 || track_index[t] != track_index[t - 1]) {
          alpha_pred = initial_dist;
        } else {
          real nm = 0;
          vector[Nstates] ap;
          for (n in 1:Nstates) nm += exp(logalpha[n, t - 1]);
          for (n in 1:Nstates) ap[n] = exp(logalpha[n, t - 1]) / nm;
          alpha_pred = tpm' * ap;
        }
        real cdf_val = 0;
        for (n in 1:Nstates) {
          real cn = (steplength[t] == 0)
            ? mixp[n]
            : mixp[n] + (1 - mixp[n]) * gamma_cdf(steplength[t] | shape[n], rate[n]);
          cdf_val += alpha_pred[n] * cn;
        }
        cdf_val = fmin(fmax(cdf_val, 1e-6), 1.0 - 1e-6);
        pseudo_residuals[t] = inv_Phi(cdf_val);
      }
    }
  }
}
