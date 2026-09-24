// Standalone generated quantities for integrated PSIS-LOO (proportional error models 6/7/8)
// Integrates out epsilon ~ exponential(rate = 1/sigma_epsilon) analytically.
// Call via: model$generate_quantities(
//   fit$draws(variables = c("sigma_epsilon", "p_angle", "p_distance")),
//   data = golf_new_data)
//
// Analytical form: for each j, the marginal log-likelihood integrates
//   y_j ~ Binomial(n_j, p_angle_j * p_distance_j * (1 - epsilon_j))
// over epsilon_j ~ Exponential(1/sigma_epsilon).
// This is computed via numerical integration over epsilon in [0, 1].
data {
  int J;
  array[J] int n;
  array[J] int y;
}
parameters {
  real<lower=0> sigma_epsilon;
  vector<lower=0, upper=1>[J] p_angle;
  vector<lower=0, upper=1>[J] p_distance;
}
generated quantities {
  vector[J] log_lik;
  for (j in 1:J) {
    // Integrate out epsilon ~ Exponential(1/sigma_epsilon) over [0, 1]
    // using a 20-point Gauss-Legendre approximation
    real lambda = 1.0 / sigma_epsilon;
    real pj = p_angle[j] * p_distance[j];
    int K = 20;
    array[20] real nodes = {
      0.0765265211334973, 0.2277858511416451, 0.3737060887154195,
      0.5108670019508271, 0.6360536807265150, 0.7463064833401390,
      0.8391169718222188, 0.9122344282513259, 0.9639719272779138,
      0.9931285991850949,
      // negative (symmetric) nodes
      -0.0765265211334973, -0.2277858511416451, -0.3737060887154195,
      -0.5108670019508271, -0.6360536807265150, -0.7463064833401390,
      -0.8391169718222188, -0.9122344282513259, -0.9639719272779138,
      -0.9931285991850949
    };
    array[20] real weights = {
      0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
      0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
      0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
      0.0176140071391521,
      0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
      0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
      0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
      0.0176140071391521
    };
    // Transform from [-1,1] to [0,1]
    real half = 0.5;
    real integral = 0.0;
    for (k in 1:K) {
      real eps = half + half * nodes[k];
      real w   = half * weights[k];
      if (eps > 0 && eps < 1) {
        real log_binom = binomial_lpmf(y[j] | n[j], pj * (1 - eps));
        real log_exp   = exponential_lpdf(eps | lambda);
        integral += w * exp(log_binom + log_exp);
      }
    }
    log_lik[j] = log(integral + 1e-300);
  }
}
