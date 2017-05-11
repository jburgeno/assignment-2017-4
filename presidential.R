# Sheridan Grant
# 5/9/2017
# Regularization Example

## Packages
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(glmnet)

## Data
# Read
dat <- read.csv('data/presdata.csv')
head(dat)
# Clean
dat <- dat[complete.cases(dat),-1]
head(dat)

# For quick model fits/easy interpretation, we'll select a subset of variables to predict with
x <- as.matrix(cbind(dat$septpoll, dat$gdpqtr2half, dat$q2gdp, dat$julypop, dat$fatalities, dat$term))
y <- dat$dv
n <- dim(dat)[1]
p <- dim(x)[2]

# Standardize all variables so we can put the same prior on all regression coefficients
# We also need the coefficients to be centered around zero 
# if the regularization priors are to work straightforwardly
y <- (y - mean(y))/sd(y)
for (i in 1:p) {
  x[,i] <- (x[,i] - mean(x[,i]))/sd(x[,i])
}

## Stan Models
# Compilation
MLE_mod <- stan_model('stan/improper.stan')
L2_mod <- stan_model('stan/normal.stan')
L1_mod <- stan_model('stan/laplace.stan')

# Fits
MLE_fit <- sampling(MLE_mod, data = list(n = n, p = p, y = y, x = x))
lambda <- 1
L2_fit <- sampling(L2_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda))
L1_fit <- sampling(L1_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda))

# Examine
summary(MLE_fit, pars = c("alpha", "beta"))$summary
summary(L2_fit, pars = c("alpha", "beta"))$summary
summary(L1_fit, pars = c("alpha", "beta"))$summary

# What does a change in lambda do to the coefficients?
lambda_seq <- seq(0.1, 30.1, by = 6)
L2_estimates <- matrix(0, nrow = length(lambda_seq), ncol = p + 1)
L1_estimates <- matrix(0, nrow = length(lambda_seq), ncol = p + 1)
for (i in 1:length(lambda_seq)) {
  L2_estimates[i,] <- summary(sampling(L2_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda_seq[i]),
                                       refresh = -1), pars = c("alpha", "beta"))$summary[,1]
  L1_estimates[i,] <- summary(sampling(L1_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda_seq[i]),
                                       refresh = -1), pars = c("alpha", "beta"))$summary[,1]
}
matplot(lambda_seq, L2_estimates, type = 'l', main = "L2")
matplot(lambda_seq, L1_estimates, type = 'l', main = "L1")

# MAPs
# These are not the preferred summary of a Bayesian posterior
# 
L2_map <- optimizing(L2_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda))
L1_map <- optimizing(L1_mod, data = list(n = n, p = p, y = y, x = x, lambda = lambda))
MLE_map <- optimizing(MLE_mod, data = list(n = n, p = p, y = y, x = x))

# Frequentist
MLE_net <- lm(y ~ x)
MLE_net$coefficients
MLE_map$par[1:(p + 1)]
# summary(MLE_fit, pars = c("alpha", "beta"))$summary[,1]

# Haven't gotten these to work--the parameterizations are different
# and it's a hassle to work out how to make them equivalent
###################################################################
L2_net <- glmnet(x, y, alpha = 0, lambda = 1)
L2_net$beta
L2_map$par[1:(p + 1)]
# summary(L2_fit, pars = c("alpha", "beta"))$summary[,1]
plot(glmnet(x, y, alpha = 0))
abline(0,0)
print(glmnet(x, y, alpha = 0))

L1_net <- glmnet(x,y)
L1_net$beta
L1_map$par[1:(p + 1)]
plot(glmnet(x, y, alpha = 1))
print(glmnet(x, y))
###################################################################

## Prediction
# More advanced models
L1_hyper_mod <- stan_model('stan/laplace_cauchy.stan')
horse_mod <- stan_model('stan/horseshoe.stan')

# Does regularization improve predictive performance?
results <- data.frame(matrix(0, nrow = 14, ncol = 5))
colnames(results) <- c('improper', 'normal', 'laplace', 'laplace_cauchy', 'horseshoe')
# Sequential predictions
for (i in 9:n) {
  Y <- y[1:(i-1)]
  X <- x[1:(i-1),]
  improper_fit <- sampling(MLE_mod, data = list(n = i - 1, p = p, y = Y, x = X), refresh = -1)
  # improper_fit <- sampling(MLE_mod, data = list(n = i - 1, p = p, y = Y, x = X))
  improper_pred <- mean(extract(improper_fit, 'alpha')$alpha) + x[i,]%*%apply(extract(improper_fit, 'beta')$beta, 2, mean)
  results$improper[i-1] <- (y[i] - improper_pred)^2
  normal_fit <- sampling(L2_mod, data = list(n = i - 1, p = p, y = Y, x = X, lambda = 0.3), refresh = -1)
  normal_pred <- mean(extract(normal_fit, 'alpha')$alpha) + x[i,]%*%apply(extract(normal_fit, 'beta')$beta, 2, mean)
  results$normal[i-1] <- (y[i] - normal_pred)^2
  laplace_fit <- sampling(L1_mod, data = list(n = i - 1, p = p, y = Y, x = X, lambda = 0.3), refresh = -1)
  laplace_pred <- mean(extract(laplace_fit, 'alpha')$alpha) + x[i,]%*%apply(extract(laplace_fit, 'beta')$beta, 2, mean)
  results$laplace[i-1] <- (y[i] - laplace_pred)^2
  laplace_cauchy_fit <- sampling(L1_hyper_mod, data = list(n = i - 1, p = p, y = Y, x = X), refresh = -1)
  laplace_cauchy_pred <- mean(extract(laplace_cauchy_fit, 'alpha')$alpha) + x[i,]%*%apply(extract(laplace_cauchy_fit, 'beta')$beta, 2, mean)
  results$laplace_cauchy[i-1] <- (y[i] - laplace_cauchy_pred)^2
  horse_fit <- sampling(horse_mod, data = list(n = i - 1, p = p, y = Y, x = X, df_local = 3, df_global = 3), refresh = -1)
  horse_pred <- mean(extract(horse_fit, 'alpha')$alpha) + x[i,]%*%apply(extract(horse_fit, 'beta')$beta, 2, mean)
  results$horseshoe[i-1] <- (y[i] - horse_pred)^2
}