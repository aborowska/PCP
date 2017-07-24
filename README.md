# Partially Censored Posterior project

_*\_MC\_fun.m_ are the functions of (T, sigma2, S, II) to be called for S simulations of time series of length T with the split normal distribution for the DGP where the left tail has standard deviation sigma2, II (a "thinning" parameter) has default value 10 (every 10th draw from the regular posterior is taken to conditionally generate 10 draws from the partially censored posterior).

_*\_run.m_ and _*\_run\_var.m_ are called by _*\_MC\_fun.m_ to run a single simulation with time-constant and time-varying threshold, respectively.