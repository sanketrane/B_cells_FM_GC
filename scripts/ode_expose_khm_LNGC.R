require(rstan)

expose_stan_functions("models/YFP_informed_models/LNGC/ki67_KHM_LNGC_FM.stan")

parms <- c(psi = quantile(matrix_of_draws$psi, probs = 0.975),
           f_fast = quantile(matrix_of_draws$f_fast, probs = 0.975),
           rhofast = quantile(matrix_of_draws$rhoFast, probs = 0.975),
           rhoSlow = quantile(matrix_of_draws$rhoSlow, probs = 0.975),
           Beta = quantile(matrix_of_draws$Beta, probs = 0.975),
           lambdaFast = quantile(matrix_of_draws$lambdaFast, probs = 0.975),
           lambdaSlow = quantile(matrix_of_draws$lambdaSlow, probs = 0.975))

alpha_c <- quantile(matrix_of_draws$alpha_calc, probs = 0.975)

# setting data input according to the precurosr pop used
if (grepl("T2", data_derived1) == TRUE){
  source_list = data.frame("nu" = -5.53e-05, "theta0" = 14.13, "chiEst" = 0.73, "qEst" = 0.06, "eps_donor" = 0.82, "eps_host" = 0.73)
} else if (grepl("T1", data_derived1) == TRUE) {
  source_list = data.frame("nu" = 1e-04, "theta0" = 13.81, "chiEst" = 0.78, "qEst" = 0.07, "eps_donor" = 0.97, "eps_host" = 0.94)
} else if (grepl("FM", data_derived1) == TRUE){
  source_list = data.frame("nu" = -1.1e-03 , "theta0" = 17.20, "chiEst" = 0.74, "qEst" = 0.023, "eps_donor" = 0.31, "eps_host" = 0.14)
} else {
  source_list = data.frame("nu" = 3.30e-04, "theta0" = 14.83, "chiEst" = 0.75, "qEst" = 0.066, "eps_donor" = 0.88, "eps_host" = 0.81)
}

rdata <- c(theta0 = exp(source_list$theta0), nu = source_list$nu, chiEst = source_list$chiEst, qEst = source_list$qEst, 
           eps_donor = source_list$eps_donor, eps_host = source_list$eps_host)

init_cond <- c(0, 0, 0, 0,
               exp(quantile(matrix_of_draws$y0_Log, probs = 0.975)) * alpha_c * quantile(matrix_of_draws$kappaF_0, probs = 0.975),
               exp(quantile(matrix_of_draws$y0_Log, probs = 0.975)) * alpha_c * (1 - quantile(matrix_of_draws$kappaF_0, probs = 0.975)),
               exp(quantile(matrix_of_draws$y0_Log, probs = 0.975)) * (1 - alpha_c) * quantile(matrix_of_draws$kappaS_0, probs = 0.975),
               exp(quantile(matrix_of_draws$y0_Log, probs = 0.975)) * (1 - alpha_c) * (1 - quantile(matrix_of_draws$kappaS_0, probs = 0.975)))


ode_df1 <- solve_ode_pred(data$ts_pred1, init_cond, parms, rdata, data$tb_pred1, data$numPred1, data$tb_time_pred1)
ode_df2 <- solve_ode_pred(data$ts_pred2, init_cond, parms, rdata, data$tb_pred2, data$numPred2, data$tb_time_pred2)
ode_df3 <- solve_ode_pred(data$ts_pred3, init_cond, parms, rdata, data$tb_pred3, data$numPred3, data$tb_time_pred3)


stan_pred_df1 <- data.frame("time" = ts_pred1,
                            "y_pred" = matrix(unlist(ode_df1), nrow = length(ode_df1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4 + y_pred.5 + y_pred.6 +  y_pred.7 + y_pred.8,
         f_fast = (y_pred.1 + y_pred.2 + y_pred.5 + y_pred.6)/ total_counts,
         f_slow = (y_pred.3 + y_pred.4 + y_pred.7 + y_pred.8)/ total_counts)

stan_pred_df1$f_fast[stan_pred_df1$time == 75]
stan_pred_df1$f_slow[stan_pred_df1$time == 75]
stan_pred_df1$f_fast[stan_pred_df1$time == 300]
stan_pred_df1$f_slow[stan_pred_df1$time == 300]


stan_pred_df2 <- data.frame("time" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_df2), nrow = length(ode_df2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)

stan_pred_df3 <- data.frame("time" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_df3), nrow = length(ode_df3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)
