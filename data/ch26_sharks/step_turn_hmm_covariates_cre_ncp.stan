// step_turn_hmm_covariates_cre_ncp.stan
//
// 2-state HMM with:
//   - Covariates in the transition probability matrix (time-inhomogeneous)
//   - Individual-level random effects in the tpm (non-centered parameterisation)
//
// Random effects model for the tpm intercept of each shark:
//   randeff_tpm[shark, k] = mu_tpm[k] + sigma_re[k] * z_re[shark, k]
// where k indexes off-diagonal transition rows (k = 1: 1->2, k = 2: 2->1).
//
// The covariates matrix excludes the intercept column (covs has nCovs columns),
// because the individual random effect replaces the population intercept.
//
// Source: Ch 26 of Gelman, Vehtari et al. (2026) Bayesian Workflow.

data {
  int<lower=1> Nstates;
  int<lower=1> Tlen;
  int<lower=1> Ntracks;
  array[Tlen] int<lower=1, upper=Ntracks> track_index;
  vector[Tlen] steplength;
  vector[Tlen] angle;
  int<lower=1> nCovs;                          // covariates (NO intercept column)
  matrix[Tlen, nCovs] covs;                    // design matrix without intercept
  int<lower=1> Nsharks;                        // number of unique individuals
  array[Tlen] int<lower=1, upper=Nsharks> shark_index;
}

parameters {
  positive_ordered[Nstates] mu;
  vector<lower=0>[Nstates] sigma;
  vector<lower=0, upper=1>[Nstates] mixp;
  vector[Nstates] xangle;
  vector[Nstates] yangle;
  simplex[Nstates] initial_dist;

  // Fixed covariate effects on tpm off-diagonal rows (no intercept)
  matrix[Nstates * (Nstates - 1), nCovs] beta;

  // Population-level intercepts for tpm
  vector[Nstates * (Nstates - 1)] mu_tpm;

  // Non-centred random effects
  vector<lower=0>[Nstates * (Nstates - 1)] sigma_re;
  matrix[Nsharks, Nstates * (Nstates - 1)] z_re;   // raw (standard-normal) random effects
}

transformed parameters {
  vector<lower=0>[Nstates] shape;
  vector<lower=0>[Nstates] rate;
  vector<lower=-pi(), upper=pi()>[Nstates] loc;
  vector<lower=0>[Nstates] kappa;

  // Individual random effects in tpm (non-centred)
  matrix[Nsharks, Nstates * (Nstates - 1)] randeff_tpm;
  for (s in 1:Nsharks)
    for (k in 1:(Nstates * (Nstates - 1)))
      randeff_tpm[s, k] = mu_tpm[k] + sigma_re[k] * z_re[s, k];

  for (n in 1:Nstates) {
    shape[n] = mu[n] * mu[n] / (sigma[n] * sigma[n]);
    rate[n]  = mu[n]       / (sigma[n] * sigma[n]);
    loc[n]   = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]^2 + yangle[n]^2);
  }
}

model {
  // --- Priors ---
  mu      ~ normal(0.2, 0.3);
  sigma   ~ student_t(3, 0, 0.5);
  mixp    ~ beta(1, 5);
  xangle  ~ normal(0, 2);
  yangle  ~ normal(0, 1);
  mu_tpm  ~ normal(0, 2);
  sigma_re ~ student_t(3, 0, 1);
  to_vector(z_re) ~ std_normal();
  for (r in 1:(Nstates * (Nstates - 1)))
    beta[r] ~ normal(0, 2);
  initial_dist ~ dirichlet(rep_vector(1.0, Nstates));

  // --- Forward algorithm with individual-specific time-varying tpm ---
  vector[Nstates] lp;
  vector[Nstates] lp_p1;

  for (t in 1:Tlen) {
    int si = shark_index[t];
    // Build tpm: intercept from randeff_tpm[si,] + fixed covariate effects
    matrix[Nstates, Nstates] tpm_t;
    matrix[Nstates, Nstates] log_tpm_tr_t;
    {
      int beta_row = 1;
      for (i in 1:Nstates) {
        for (j in 1:Nstates) {
          if (i == j)
            tpm_t[i, j] = 1.0;
          else {
            tpm_t[i, j] = exp(randeff_tpm[si, beta_row]
                              + covs[t] * beta[beta_row]');
            beta_row += 1;
          }
        }
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
