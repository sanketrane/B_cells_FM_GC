require(rstan)

expose_stan_functions("models/Miscellaneous/TDD_trial.stan")

# time sequence for prediction
ts_pred1 <- seq(48, to = 600, by = 1)
tb_pred1 <- rep(48, length(ts_pred1))
ts_pred2 <- seq(72, to = 600, by = 1)
tb_pred2 <- rep(72, length(ts_pred2))
ts_pred3 <- seq(130, to = 600, by = 1)
tb_pred3 <- rep(130, length(ts_pred3))

num_pred1 <- length(ts_pred1)
num_pred2 <- length(ts_pred2)
num_pred3 <- length(ts_pred3)

tb_time_pred1 <- c(40, 48)
tb_time_pred2 <- c(40, 72)
tb_time_pred3 <- c(40, 130)

init_cond <- c(0, 0, 1e6, 1e6)
parms <- c(psi = 0.5, rho = 0.001, Beta = 3.5,  delta = 0.031, r_psi = 0.001)
rdata <- c(theta0 = exp(14.13), nu = -5.53e-05, chiEst = 0.73, qEst = 0.06, eps_donor = 0.82, eps_host = 0.73)

#shm_df <- solve_init(tb_time_pred1, init_cond, parms, rdata, 40, 2)

#stan_pred_df <- data.frame("time" = tb_time_pred1,
#                            "y_pred" = matrix(unlist(shm_df), nrow = length(shm_df), byrow = TRUE))%>%
#  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)


ode_df1 <- solve_ode(ts_pred1, init_cond, parms, rdata, tb_pred1, num_pred1, tb_time_pred1)
ode_df2 <- solve_ode(ts_pred2, init_cond, parms, rdata, tb_pred2, num_pred2, tb_time_pred2)
ode_df3 <- solve_ode(ts_pred3, init_cond, parms, rdata, tb_pred3, num_pred3, tb_time_pred3)


stan_pred_df1 <- data.frame("time" = ts_pred1,
                           "y_pred" = matrix(unlist(ode_df1), nrow = length(ode_df1), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)

stan_pred_df2 <- data.frame("time" = ts_pred2,
                            "y_pred" = matrix(unlist(ode_df2), nrow = length(ode_df2), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)

stan_pred_df3 <- data.frame("time" = ts_pred3,
                            "y_pred" = matrix(unlist(ode_df3), nrow = length(ode_df3), byrow = TRUE))%>%
  mutate(total_counts = y_pred.1 + y_pred.2 + y_pred.3 + y_pred.4)


ggplot()+
  geom_point(data = stan_pred_df1, aes(x = time, y= total_counts), size = 2) +
  geom_point(data = stan_pred_df2, aes(x = time, y= total_counts), col =2, size = 2) +
  geom_point(data = stan_pred_df3, aes(x = time, y= total_counts), col =3,  size = 2)
 
