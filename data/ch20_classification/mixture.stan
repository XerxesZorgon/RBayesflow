
data {
  int<lower=1> K;
  int<lower=1> N;
  array[N] real y;
  array[K] real mu;
}
parameters {
  simplex[K] theta;
  real<lower=0> sigma;
}
model {
  array[K] real ps;
  sigma ~ exponential(0.1);           // rate = 0.1, mean = 10 (weak)
  for (n in 1:N) {
    for (k in 1:K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma);
    }
    target += log_sum_exp(ps);        // marginalise the latent z_n
  }
}

