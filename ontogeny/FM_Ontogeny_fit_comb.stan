functions{
  real theta_spline(real t, real theta0, real r1, real n1) {

    real theta = theta0 * (1 + t^n1 * exp(-r1 * t));
    return theta;
   }

   real lambda(real t, real gamma){

     real value;

     real r_d = 0.00116935;          // estimated from TDT model FM fit
     real lambda0 = 0.02803697;       // estimated from TDT model FM fit
     int tb  = 40;                    // t0 for BMT (age of BMT in youngest BM recip[ient])

     real r_l;

     if (t <= tb) {
       r_l = r_d * exp(- gamma * (t-tb));
     } else {
       r_l = r_d;
     }

     value = lambda0 * exp(-r_l * (t-tb));

     return(value);
   }

  real[] tdt(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {
     real gamma   = parms[1];
     real r_psi   = 0;

     real theta0  = rdata[1];
     real r1      = rdata[2];
     real n1      = rdata[3];
     real eps     = rdata[4];
     real psi     = rdata[5];
     real rho     = rdata[6];
     real Beta    = rdata[7];

     real dkdt[2];

     dkdt[1] = (psi * (1 - exp(-r_psi * t)) * theta_spline(t, theta0, r1, n1) * eps) + rho * (2 * k[2] + k[1]) - ((1/Beta) + (lambda(t, gamma) + rho)) * k[1];
     dkdt[2] = (psi * (1 - exp(-r_psi * t)) * theta_spline(t, theta0, r1, n1) * (1-eps)) + (1/Beta) * k[1] - (rho + (lambda(t, gamma) + rho)) * k[2];

     return dkdt;
   }
}

data{
  int<lower  = 1> numObs;
  int<lower  = 1> num_index;
  real<lower = 0> solve_time[num_index];
  int<lower  = 0> time_index[numObs];
  int<lower  = 1> numPred;
  real ts_pred[numPred];
  real counts[numObs];
  real theta0;
  real r1;
  real n1;
  real eps;
  real psi;
  real rho;
  real Beta;
  real y0_Log;
  real kappa0;
}

transformed data{
  real y[numObs];
  real rdata[9];
  int idata[0];

  for (i in 1:numObs){
    y[i] = log(counts[i]);
  }

  rdata[1] = theta0;
  rdata[2] = r1;
  rdata[3] = n1;
  rdata[4] = eps;
  rdata[5] = psi;
  rdata[6] = rho;
  rdata[7] = Beta;
  rdata[8] = y0_Log;
  rdata[9] = kappa0;
}

parameters{
  real gamma_Log;
  //real y0_Log;

  real<lower = 0> sigma;
}

transformed parameters{
  real k_hat[num_index, 2];
  real ymean[numObs];
  real init_cond[2];
  real parms[1];

  real y0 = exp(y0_Log);
  real gamma = exp(gamma_Log);

  init_cond[1] = y0 * kappa0;
  init_cond[2] = y0 * (1 - kappa0);

  parms[1] = gamma;

  k_hat[1, ] = init_cond;
  k_hat[2:num_index, ] = integrate_ode_rk45(tdt, init_cond, solve_time[1], solve_time[2:num_index], parms, rdata, idata);

  for (i in 1:numObs){
    ymean[i] = k_hat[time_index[i], 1] + k_hat[time_index[i], 2];
  }
}

model{
  gamma_Log ~ normal(-4, 1);
  //y0_Log ~ normal(14, 1);

  sigma ~ normal(0.5, 1);

  y ~ normal(log(ymean), sigma);
}

generated quantities{
  real k_hat_pred[numPred, 2];
  real ymean_pred[numPred];
  real countspred[numPred];
  vector[numObs] log_lik;

  // Predictions
  k_hat_pred[1, ] = init_cond;
  k_hat_pred[2:numPred, ] = integrate_ode_rk45(tdt, init_cond, ts_pred[1], ts_pred[2:numPred], parms, rdata, idata);


  ymean_pred[1] = y0;
  countspred[1] = exp(normal_rng(y0, sigma));
  // post-predicive check
  for (i in 2:numPred){
    // mean
    ymean_pred[i]  = k_hat_pred[i, 1] + k_hat_pred[i, 2];

    // intervals
    countspred[i] = exp(normal_rng(log(ymean_pred[i]), sigma));
  }

  //  log likelihood
  for (n in 1:numObs) {
    log_lik[n] = normal_lpdf(y[n] | log(ymean[n]), sigma);
  }
}
