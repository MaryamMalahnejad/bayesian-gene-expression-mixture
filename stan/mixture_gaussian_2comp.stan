data {
  int<lower=1> N;
  vector[N] x;
}
parameters {
  real mu1;
  real mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0, upper=1> theta;
}
model {
  mu1 ~ normal(0, 5);
  mu2 ~ normal(0, 5);
  sigma1 ~ exponential(1);
  sigma2 ~ exponential(1);
  theta ~ beta(2, 2);

  for (n in 1:N) {
    target += log_mix(theta,
      normal_lpdf(x[n] | mu1, sigma1),
      normal_lpdf(x[n] | mu2, sigma2));
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = log_mix(theta,
      normal_lpdf(x[n] | mu1, sigma1),
      normal_lpdf(x[n] | mu2, sigma2));
  }
}