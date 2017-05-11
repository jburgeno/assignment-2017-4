data {
  // number of observations
  int n;
  // response vector
  vector[n] y;
  // number of columns in the design matrix X
  int p;
  // design matrix X
  matrix [n, p] x;
  // degress of freedom for local and global scales
  real<lower = 0.> df_local;
  real<lower = 0.> df_global;
}
transformed data {
  real<lower = 0.> y_sd;
  real a_pr_scale;
  real sigma_pr_scale;
  y_sd = sd(y);
  sigma_pr_scale = y_sd * 5.;
  a_pr_scale = 10.;
}
parameters {
  // regression coefficient vector
  real alpha;
  vector[p] b_raw;
  // scale of the regression errors
  real<lower = 0.> sigma;
  // guess at non-zero coefficients
  real<lower = 0.> p0;
  // local scales of coefficients
  vector<lower = 0.>[p] lambda;
  // glboal scale of coefficients
  real<lower = 0.> tau;
}
transformed parameters {
  // mu is the observation fitted/predicted value
  // also called yhat
  vector[n] mu;
  vector[p] beta;
  real<lower = 0.> tau0;
  // tau0 from Piironen and Vehtari (2017)
  tau0 = p0 * (p - p0) * sigma * pow(n, 0.5);
  beta = b_raw * tau .* lambda;
  mu = alpha + x * beta;
}
model {
  // priors
  lambda ~ student_t(df_local, 0., 1.);
  alpha ~ normal(0., a_pr_scale);
  b_raw ~ normal(0., 1.);
  tau ~ student_t(df_global, 0., tau0);
  sigma ~ cauchy(0., sigma_pr_scale);
  // likelihood
  y ~ normal(mu, sigma);
}
generated quantities {
  // simulate data from the posterior
  vector[n] y_rep;
  // log-likelihood posterior
  vector[n] log_lik;
  for (i in 1:n) {
    y_rep[i] = normal_rng(mu[i], sigma);
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
  }
}
