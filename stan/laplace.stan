data {
  // number of observations
  int n;
  // response vector
  vector[n] y;
  // number of columns in the design matrix X
  int p;
  // design matrix X
  matrix [n, p] x;
  // regularization parameter
  real<lower = 0.> lambda;
}
transformed data {
  real<lower = 0.> tau = 1 / lambda;
}
parameters {
  // regression coefficient vector
  real alpha;
  vector[p] beta;
  // scale of the regression errors
  real<lower = 0.> sigma;
}
transformed parameters {
  // mu is the observation fitted/predicted value
  // also called yhat
  vector[n] mu;
  mu = alpha + x * beta;
}
model {
  // priors
  sigma ~ cauchy(0., 5);
  // alpha ~ double_exponential(0, tau * sigma);
  alpha ~ double_exponential(0, tau);
  // beta ~ double_exponential(0, tau * sigma);
  beta ~ double_exponential(0, tau);
  // likelihood
  y ~ normal(mu, sigma);
}
// generated quantities {
//   // simulate data from the posterior
//   vector[n] y_rep;
//   // log-likelihood posterior
//   vector[n] log_lik;
//   for (i in 1:n) {
//     y_rep[i] = normal_rng(mu[i], sigma);
//     log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
//   }
// }
