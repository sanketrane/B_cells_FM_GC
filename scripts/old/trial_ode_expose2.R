rm(list = ls())  
gc()    
setwd("/opt/mesh/eigg/sanket/FM_GC_ageCorrected")

require(rstan)
require(tidyverse)

expose_stan_functions("models/Miscellaneous/TDD_trial.stan")

set.seed(4009)

Chi_T1_vec <- Vectorize(Chi_T1)

# time sequence for prediction
ts_pred1 <- seq(48, to = 600, by = 1)
num_pred1 <- length(ts_pred1)
tb_pred1 <- c(40, rep(48, (num_pred1 -1)))

ts_pred2 <- seq(72, to = 600, by = 1)
num_pred2 <- length(ts_pred2)
tb_pred2 <- c(40, rep(72, (num_pred2 -1)))

ts_pred3 <- seq(130, to = 600, by = 1)
num_pred3 <- length(ts_pred3)
tb_pred3 <- c(40, rep(130, (num_pred3 -1)))


tb_time_pred1 <- c(40, rep(48, (num_pred1 -1)))
tb_time_pred2 <- c(40, rep(72, (num_pred1 -1)))
tb_time_pred3 <- c(40, rep(130, (num_pred1 -1)))

tb_index_pred1 <- seq(1, num_pred1)
tb_index_pred2 <- seq(1, num_pred2)
tb_index_pred3 <- seq(1, num_pred3)

parms_list <- list(
  psi = exp(rnorm(1, log(0.01), 0.5)),
  f_fast = exp(rnorm(1,log(0.3), 0.01)),
  rhoFast = exp(rnorm(1,log(0.1), 0.01)),
  rhoSlow = exp(rnorm(1,log(0.01), 0.01)),
  lambdaSlow = exp(rnorm(1, log(0.001), 0.01)),
  lambdaFast = exp(rnorm(1, log(0.01), 0.01)),
  deltaYFP = exp(rnorm(1, log(0.4), 0.01)),
  Beta = exp(rnorm(1, log(3.5), 0.1)),
  
  y0_Log = rnorm(1, 11 , 0.1),
  kappaF_0 = exp(rnorm(1, log(0.9), 0.01)),
  kappaS_0 = exp(rnorm(1, log(0.9), 0.01)),
  
  sigma1 = exp(rnorm(1,log(1.5), 1)),
  sigma2 = exp(rnorm(1,log(1.5), 1)),
  sigma3 = exp(rnorm(1,log(1.5), 1)),
  sigma4 = exp(rnorm(1,log(1.5), 1))
)

alpha = (parms_list$deltaYFP - exp(-parms_list$lambdaSlow * 58)) /(exp(-parms_list$lambdaFast * 58) - exp(-parms_list$lambdaSlow * 58))


init_cond <- c(0, 0, 0, 0,  alpha * exp(parms_list$y0_Log) * parms_list$kappaF_0, alpha * exp(parms_list$y0_Log) * (1 - parms_list$kappaF_0),
               (1 - alpha) * exp(parms_list$y0_Log) * parms_list$kappaS_0, (1 - alpha) * exp(parms_list$y0_Log) * (1 - parms_list$kappaS_0))

parms <- c(psi = parms_list$psi, f_fast = parms_list$f_fast, rhoFast = parms_list$rhoFast, rhoSlow = parms_list$rhoSlow,
           Beta = parms_list$Beta, lambdaFast = parms_list$lambdaFast, lambdaSlow = parms_list$lambdaSlow)

rdata <- c(theta0 = exp(14.13), nu = -5.53e-05, chiEst = 0.73, qEst = 0.06, eps_donor = 0.82, eps_host = 0.73)

#predictions
ode_df1 <- solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, num_pred1, tb_time_pred1)
ode_df2 <- solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, num_pred2, tb_time_pred2)
ode_df3 <- solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, num_pred3, tb_time_pred3)

stan_pred_df1 <- data.frame("time" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_df1), nrow = length(ode_df1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/(total_counts * Chi_T1_vec(time)),
         ki_d = y_pred.1/(y_pred.1 + y_pred.2),
         ki_h = y_pred.3/(y_pred.3 + y_pred.4))

stan_pred_df2 <- data.frame("time" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_df2), nrow = length(ode_df2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/(total_counts * Chi_T1_vec(time)),
         ki_d = y_pred.1/(y_pred.1 + y_pred.2),
         ki_h = y_pred.3/(y_pred.3 + y_pred.4))

stan_pred_df3 <- data.frame("time" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_df3), nrow = length(ode_df3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 + y_pred.7 + y_pred.8,
         fd = (y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)/(total_counts * Chi_T1_vec(time)),
         ki_d = y_pred.1/(y_pred.1 + y_pred.2),
         ki_h = y_pred.3/(y_pred.3 + y_pred.4))


ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= total_counts), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= total_counts), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= total_counts), col =3,  size = 2) + scale_y_log10()

ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= fd), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= fd), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= fd), col =3,  size = 2) + ylim(0,1)

ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= ki_d), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= ki_d), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= ki_d), col =3,  size = 2) + ylim(0,1)

ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= ki_h), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= ki_h), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= ki_h), col =3,  size = 2) + ylim(0,1)


###################################################################################################################################################
###################################################################################################################################################
