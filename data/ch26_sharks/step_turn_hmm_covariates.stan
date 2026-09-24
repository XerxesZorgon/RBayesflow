// step_turn_hmm_covariates.stan
//
// Time-inhomogeneous 2-state HMM: covariates in the transition probability matrix.
// Covariates (chum, sex, tod_cos, tod_sin) enter via a multinomial-logit link:
//
//   eta[i,j,t] = covs[t,] %*% beta[row_ij,]
//   tpm[i,j,t] = softmax(eta[i,,t])[j]
//
// One beta row per off-diagonal transition (for N=2: 2 rows).
// State-dependent distributions identical to step_turn_hmm.stan.
// The covariates matrix includes a leading intercept column (column 1 = 1).
//
// Source: Ch 26 of Gelman, Vehtari et al. (2026) Bayesian Workflow.

data {
  int<lower=1> Nstates;
  int<lower=1> Tlen;
  int<lower=1> Ntracks;
  array[Tlen] int<lower=1, upper=Ntracks> track_index;
  vector[Tlen] steplength;
  vector[Tlen] angle;
  int<lower=1> nCovs;                       // number of covariates (excl. intercept)
  matrix[Tlen, nCovs + 1] covs;             // design matrix: intercept in col 1
}

parameters {
  positive_ordered[Nstates] mu;
  vector<lower=0>[Nstates] sigma;
  vector<lower=0, upper=1>[Nstates] mixp;
  vector[Nstates] xangle;
  vector[Nstates] yangle;
  matrix[Nstates * (Nstates - 1), nCovs + 1] beta; // regression coefficients
  simplex[Nstates] initial_dist;
}

transformed parameters {
  vector<lower=0>[Nstates] shape;
  vector<lower=0>[Nstates] rate;
  vector<lower=-pi(), upper=pi()>[Nstates] loc;
  vector<lower=0>[Nstates] kappa;

  for (n in 1:Nstates) {
    shape[n] = mu[n] * mu[n] / (sigma[n] * sigma[n]);
    rate[n]  = mu[n]       / (sigma[n] * sigma[n]);
    loc[n]   = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]^2 + yangle[n]^2);
  }
}

model {
  // --- Priors ---
  mu     ~ normal(0.2, 0.3);
  sigma  ~ student_t(3, 0, 0.5);
  mixp   ~ beta(1, 5);
  xangle ~ normal(0, 2);
  yangle ~ normal(0, 1);
  for (r in 1:(Nstates * (Nstates - 1)))
    beta[r] ~ normal(0, 2);
  initial_dist ~ dirichlet(rep_vector(1.0, Nstates));

  // --- Forward algorithm with time-varying tpm ---
  vector[Nstates] lp;
  vector[Nstates] lp_p1;

  for (t in 1:Tlen) {
    // Compute log-tpm at time t
    matrix[Nstates, Nstates] tpm_t;
    matrix[Nstates, Nstates] log_tpm_tr_t;
    {
      int beta_row = 1;
      for (i in 1:Nstates) {
        for (j in 1:Nstates) {
          if (i == j)
            tpm_t[i, j] = 1.0;
          else {
            tpm_t[i, j] = exp(covs[t] * beta[beta_row]');
            beta_row += 1;
          }
        }
        // Normalise row
        real row_sum = sum(tpm_t[i]);
        for (j in 1:Nstates)
          log_tpm_tr_t[j, i] = log(tpm_t[i, j] / row_sum);
      }
    }

    if (t == 1 || track_index[t] != track_index[t - 1]) {
      for (n in 1:Nstates) lp[n] = log(initial_dist[n]);
    }
    for (n in 1:Nstates) {
      lp_p1[n] = log_sum_exp(to_row_vector(lp) + log_tpm_tr_t[n]);
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
