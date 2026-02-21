data {
  int<lower=1> N;
  vector[N] x;
}

parameters {
  ordered[2] mu;
  vector<lower=0>[2] sigma;
  simplex[2] theta;
}

model {
  mu ~ normal(0, 5);
  sigma ~ exponential(1);
  theta ~ dirichlet(rep_vector(2.0, 2));

  for (n in 1:N) {
    target += log_mix(theta[1],
      normal_lpdf(x[n] | mu[1], sigma[1]),
      normal_lpdf(x[n] | mu[2], sigma[2]));
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] x_rep;
  array[N] int z_rep;

  for (n in 1:N) {
    log_lik[n] = log_mix(theta[1],
      normal_lpdf(x[n] | mu[1], sigma[1]),
      normal_lpdf(x[n] | mu[2], sigma[2]));

    z_rep[n] = categorical_rng(theta);
    x_rep[n] = normal_rng(mu[z_rep[n]], sigma[z_rep[n]]);
  }
}