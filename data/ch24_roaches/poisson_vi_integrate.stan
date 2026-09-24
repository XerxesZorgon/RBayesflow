// Poisson regression with hierarchical intercept (varying effect)
// Integrated LOO log-likelihood via 1D adaptive quadrature.
// v4: fmin() for scalar cap (Stan min() is array-only).
functions {
  real integrand(real z, real notused, array[] real theta,
                 array[] real X_i, array[] int y_i) {
    real sigmaz = theta[1];
    real mu_i   = theta[2];
    real lp = normal_lpdf(z | 0, sigmaz) +
              poisson_log_lpmf(y_i | z + mu_i);
    if (is_inf(lp) || is_nan(lp)) return 0.0;
    real p = exp(lp);
    return (is_inf(p) || is_nan(p)) ? 0.0 : p;
  }
}
data {
  int<lower=0> N;
  int<lower=0> P;
  matrix[N, P] X;
  array[N] int<lower=0> y;
  vector[N] offsett;
  real integrate_1d_reltol;
}
parameters {
  real alpha;
  vector[P] beta;
  vector[N] z;
  real<lower=0> sigmaz;
}
model {
  alpha  ~ normal(0, 1);
  beta   ~ normal(0, 1);
  z      ~ normal(0, sigmaz);
  sigmaz ~ normal(0, 1);
  y ~ poisson_log_glm(X, z + offsett + alpha, beta);
}
generated quantities {
  vector[N] log_lik;
  vector[N] y_rep;
  vector[N] y_loorep;
  for (i in 1:N) {
    real mu_i = offsett[i] + alpha + X[i,] * beta;

    // Cap log-rate at 20 (e^20 ~ 5e8) to keep poisson_log_rng within int range
    y_rep[i]    = poisson_log_rng(fmin(z[i] + mu_i, 20.0));

    // Integrated LOO log-likelihood
    real integral = integrate_1d(
      integrand,
      negative_infinity(), positive_infinity(),
      append_array({sigmaz}, {mu_i}),
      {0},
      {y[i]},
      integrate_1d_reltol
    );
    if (integral > 0 && !is_nan(integral) && !is_inf(integral)) {
      log_lik[i] = log(integral);
    } else {
      log_lik[i] = poisson_log_lpmf(y[i] | z[i] + mu_i);
    }

    // LOO predictive replicate (z drawn from prior)
    y_loorep[i] = poisson_log_rng(fmin(normal_rng(0, sigmaz) + mu_i, 20.0));
  }
}
