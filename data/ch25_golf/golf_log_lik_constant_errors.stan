// Standalone generated quantities for PSIS-LOO (constant error models 9/10/11)
// epsilon is a scalar, so the binomial likelihood is available directly.
// Call via: gq_ll_c$generate_quantities(
//   fit_11$draws(variables = "p"), data = golf_new_data)
data {
  int J;
  array[J] int n;
  array[J] int y;
}
parameters {
  vector<lower=0, upper=1>[J] p;
}
generated quantities {
  vector[J] log_lik;
  for (j in 1:J) {
    log_lik[j] = binomial_lpmf(y[j] | n[j], p[j]);
  }
}
