require(rstan)


rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

expose_stan_functions("models/Miscellaneous/skewtest.stan")

distr <- data.frame("obded" = normalskew_rng(100000, 0.42, 0.1, -5))
#1/quantile(distr$obded, probs = seq(0.025, 1, 0.05))
plot(density(distr$obded), xlim = c(0, 1))


