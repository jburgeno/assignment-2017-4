# Sheridan Grant
# 5/9/2017
# Regularization Example

# Packages
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# library(rstanarm)
library(rubbish)

# Data
dat <- read.csv('data/presdata.csv')
head(dat)
dat <- dat[complete.cases(dat),-1]
head(dat)
x <- as.matrix(cbind(dat$septpoll, dat$gdpqtr2half, dat$q2gdp, dat$julypop, dat$fatalities, dat$term))
n <- dim(dat)[1]
# p <- dim(dat)[2] - 2
p <- dim(x)[2]

improper_mod <- stan_model('stan/improper.stan')
improper_predict <- numeric(n-1)

improper_fit <- sampling(improper_mod, data = list(n = n, p = p, y = dat$dv, x = x))
prediction <- apply(extract(improper_fit, "mu")$mu, 2, mean)



for (i in 2:n) {
  improper_fit <- sampling(improper_mod, data = list(n = i - 1, p = p, y = dat$dv[1:(i-1)], x = dat[1:(i-1),3:(p+2)]))
}