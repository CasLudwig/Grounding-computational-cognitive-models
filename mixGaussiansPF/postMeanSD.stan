// C.J.H. Ludwig, June 2024
// Posterior density of mean and standard deviation of a normal distribution.

data {
  int<lower = 1> N;  // Total number of trials
  vector[N] y;  // Observations
}

parameters {
  real<lower = -3, upper = 3> mu;
  real<lower = 1e-6, upper = 25> sigma2;
}

transformed parameters {
  real sigma;
  sigma = sqrt(sigma2);
}

model {
  // Priors:
  target += normal_lpdf(mu | 0, 1);
  target += scaled_inv_chi_square_lpdf(sigma2 | 5, 5);
  // Likelihood:
  target += normal_lpdf(y | mu, sigma);
}
