
data {
  int<lower=0> y_sample;
  int<lower=0> n_sample;
  int<lower=0> y_spec;
  int<lower=0> n_spec;
  int<lower=0> y_sens;
  int<lower=0> n_sens;
}
parameters {
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> spec;
  real<lower=0, upper=1> sens;
}
transformed parameters {
  real p_sample = p * sens + (1 - p) * (1 - spec);
}
model {
  y_sample ~ binomial(n_sample, p_sample);
  y_spec   ~ binomial(n_spec, spec);
  y_sens   ~ binomial(n_sens, sens);
  // uniform(0,1) priors are implicit via the parameter bounds
}
generated quantities {
  real log_lik   = binomial_lpmf(y_sample | n_sample, p_sample);
  real log_prior = binomial_lpmf(y_spec | n_spec, spec)
                 + binomial_lpmf(y_sens | n_sens, sens);
}

