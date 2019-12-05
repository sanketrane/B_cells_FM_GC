functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real spline_int) {     // spline_int gives the initial counts of the source compartment

      real theta;
      int ta = 40;         // the earliest age at BMT

      theta = spline_int * exp(-nu * (t-ta));
      return theta;
   }

   real Chi_spline( real t,
     real chiEst,         // chiEst is the level if stabilised chimerism in the source compartment
     real qEst) {         // qEst is the rate with which cimerism chnages in the source compartment

       real chi;

       if (t < 0){
         chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
       } else {
         chi = chiEst * (1 - exp(-qEst * t));
       }

         return chi;
    }

    real Chi_T1(real t) {
      real chi;
      real chiEst = 0.78;
      real qEst = 0.07;


        if (t < 0){
          chi = 0;                       // conditioning the function to adapt to the timepoints before BMT
        } else {
          chi = chiEst * (1 - exp(-qEst * t));
        }
        return chi;
     }


    real[] shm(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

     real psi    = parms[1];
     real rho    = parms[2];
     real Beta   = parms[3];
     real lambda = parms[4];

     real theta0    = rdata[1];
     real nu        = rdata[2];
     real chiEst    = rdata[3];
     real qEst      = rdata[4];
     real eps_donor = rdata[5];
     real eps_host  = rdata[6];

     real tb = idata[1];

     real dkdt[4];

     dkdt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rho * (2 * k[2] + k[1]) - ((1/Beta) + (lambda + rho)) * k[1];
     dkdt[2] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[1] - (rho + (lambda + rho)) * k[2];

     dkdt[3] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rho * (2 * k[4] + k[3]) - ((1/Beta) + (lambda + rho)) * k[3];
     dkdt[4] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[3] - (rho + (lambda + rho)) * k[4];

     return dkdt;
   }

  real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){

    return to_array_1d(integrate_ode_rk45(shm, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
  }

  real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
   real y_init[num_tb, 4];

   y_init[1] = init_cond;
   for (i in 2:num_tb){
     y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
   }

   return y_init;
  }

  real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {

    return to_array_1d(integrate_ode_rk45(shm, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
   }

  real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){

     real y_hat[num_index, 4];
     real y0[num_tb, 4];
     real init_tb[4];

     y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

     init_tb[1] = init_cond[1];                                           //at tbmt - donor ki67Hi subset size is zero
     init_tb[2] = init_cond[2];                                           //at tbmt - donor ki67Lo subset size is zero

     for (i in 1:num_index){
       init_tb[3] = y0[tb_index[i], 1] + y0[tb_index[i], 3];              //at tbmt - all ki67Hi cells would be host
       init_tb[4] = y0[tb_index[i], 2] + y0[tb_index[i], 4];              //at tbmt - all ki67Lo cells would be host
       y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
     }

     return y_hat;
  }

  real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){

     real y_hat[num_index, 4];
     real y0[2, 4];
     real init_tb[4];

     y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

     init_tb[1] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero
     init_tb[2] = init_cond[1];                                         //at tbmt - donor ki67Lo subset size is zero

     init_tb[3] = y0[2, 1] + y0[2, 3];                                  //at tbmt - all ki67Hi cells would be host
     init_tb[4] = y0[2, 2] + y0[2, 4];                                  //at tbmt - all ki67Hi cells would be host

     y_hat[1] = init_tb;
     for (i in 2:num_index){
       y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
     }

     return y_hat;
  }

  real expit(real x){
    real ans;

    ans = exp(x)/(1+exp(x));

    return ans;
  }

  real[] logit_boundary_array(real[] x){
    int ndims = size(x);
    real answer[ndims];
    real b = 1.5;

    for (i in 1: ndims){
      answer[i] = log(x[i]/(b-x[i]));
    }

    return answer;
  }

  real logit_boundary(real x){
    real answer;
    real b = 1.5;

    answer = log(x/(b-x));

    return answer;
  }

  real[] expit_boundary_array(real[] x){
    int ndims = size(x);
    real answer[ndims];
    real b = 1.5;

    for (i in 1: ndims){
      answer[i] =  (b * exp(x[i])) /( 1+ exp(x[i]));
    }

    return answer;
  }

  real expit_boundary(real x){
    real answer;
    real b = 1.5;

    answer =  (b * exp(x)) /( 1+ exp(x));

    return answer;
  }
}

data{
  int<lower  = 1> numObs;                           // number of observations for donor fractions may be different that cell counts
  int<lower  = 1> num_index;
  int<lower  = 1> num_tb;
  real<lower = 0> solve_time[num_index];
  int<lower  = 1> time_index[numObs];
  int<lower  = 1> tb_index[num_index];
  int<lower = 0> ageAtBMT[num_index];
  real<lower =0> tb_time[num_tb];
  int dpBMT[numObs];
  real Nd_0;
  real<lower = 0> counts[numObs];
  real<lower = 0> Nfd[numObs];
  real<lower = 0> ki_donor[numObs];
  real<lower = 0> ki_host[numObs];
  int<lower  = 1> numPred1;
  int<lower  = 1> numPred2;
  int<lower  = 1> numPred3;
  real ts_pred1[numPred1];
  real ts_pred2[numPred2];
  real ts_pred3[numPred3];
  int tb_pred1[numPred1];
  int tb_pred2[numPred2];
  int tb_pred3[numPred3];
  real tb_time_pred1[2];
  real tb_time_pred2[2];
  real tb_time_pred3[2];
  real theta0;
  real nu;
  real chiEst;
  real qEst;
  real eps_donor;
  real eps_host;
  }

transformed data{
  real y1[numObs];
  real y2[numObs];
  real y3[numObs];
  real y4[numObs];
  real rdata[6];
  int tb[num_index];

  rdata[1] = theta0;
  rdata[2] = nu;
  rdata[3] = chiEst;
  rdata[4] = qEst;
  rdata[5] = eps_host;
  rdata[6] = eps_donor;

  for (i in 1:num_index){
    tb[i] = ageAtBMT[i];
  }

  y1 = log(counts);                 // transforming cell counts of donor compartments to feed in to ODEs
  y2 = (Nfd);                       // untransfored cell counts of donor fractions normalised to source chimerism to feed in to ODEs
  y3 = logit_boundary_array(ki_donor);             // transforming counts of ki67 positive cells in the  donor compartments to feed in to ODEs
  y4 = logit_boundary_array(ki_host);              // transforming counts of ki67 positive cells in the  host compartments to feed in to ODEs
}

parameters{
  real y0_Log;
  real<lower = 0, upper=1> kappa_0;
  real<lower = 0, upper=1> psi;
  real<lower = 0> rho;
  //real lambda;
  real<lower = 0, upper=1> deltaYFP;
  real<lower = 0> Beta;

  real<lower = 0> sigma1;
  real<lower = 0> sigma2;
  real<lower = 0> sigma3;
  real<lower = 0> sigma4;
  }

transformed parameters{
  real k_hat[num_index, 4];
  real y1_mean[numObs];
  real y2_mean[numObs];
  real y3_mean[numObs];
  real y4_mean[numObs];
  real parms[4];
  real init_cond[4];
  real y0 = exp(y0_Log);
  real lambda_calc = -log(deltaYFP)/ 41;

  init_cond[1] = Nd_0;
  init_cond[2] = Nd_0;
  init_cond[3] = y0 * kappa_0;
  init_cond[4] = y0 * (1 - kappa_0);
  parms[1] = psi;
  parms[2] = rho;
  parms[3] = Beta;
  parms[4] = lambda_calc;

  k_hat = solve_ode(solve_time, init_cond, parms, rdata, tb, num_index, tb_time, tb_index, num_tb);

  for (i in 1:numObs){
    // total counts
    y1_mean[i] = k_hat[time_index[i], 1] + k_hat[time_index[i], 2] + k_hat[time_index[i], 3] + k_hat[time_index[i], 4];

    // donor fractions normalised with chimerism in the source
    y2_mean[i] = (k_hat[time_index[i], 1] + k_hat[time_index[i], 2])/(y1_mean[i] * Chi_T1(dpBMT[i]));

    // fractions of ki67 positive cells in the donor compartment
    y3_mean[i] = (k_hat[time_index[i], 1])/(k_hat[time_index[i], 1] + k_hat[time_index[i], 2]);

    // fractions of ki67 positive cells in the host compartment
    y4_mean[i] = (k_hat[time_index[i], 3])/(k_hat[time_index[i], 3] + k_hat[time_index[i], 4]);
  }
}

model{
  psi ~ normal(0.5, 0.25);
  y0_Log ~ normal(17, 2);
  kappa_0 ~  normal(0.2, 0.15);
  rho ~ normal(0.01, 0.5);
  deltaYFP ~ normal(0.3, 0.1);
  Beta ~ normal(3.5, 1);

  sigma1 ~ normal(1, 1);
  sigma2 ~ normal(1.5, 1);
  sigma3 ~ normal(1, 1);
  sigma4 ~ normal(0.5, 1);

  y1 ~ normal(log(y1_mean), sigma1);
  y2 ~ normal((y2_mean), sigma2);
  y3 ~ normal(logit(y3_mean), sigma3);
  y4 ~ normal(logit(y4_mean), sigma4);
}

generated quantities{
  real k_hat_pred_age1[numPred1, 4];
  real k_hat_pred_age2[numPred2, 4];
  real k_hat_pred_age3[numPred3, 4];
  real y1_mean_pred_age1[numPred1];
  real y2_mean_pred_age1[numPred1];
  real y1_mean_pred_age2[numPred2];
  real y2_mean_pred_age2[numPred2];
  real y1_mean_pred_age3[numPred3];
  real y2_mean_pred_age3[numPred3];
  real y3_mean_pred1[numPred1];
  real y4_mean_pred1[numPred1];
  real y3_mean_pred2[numPred2];
  real y4_mean_pred2[numPred2];
  real y3_mean_pred3[numPred3];
  real y4_mean_pred3[numPred3];
  real countspred_age1[numPred1];
  real fdpred_age1[numPred1];
  real countspred_age2[numPred2];
  real fdpred_age2[numPred2];
  real countspred_age3[numPred3];
  real fdpred_age3[numPred3];
  real donor_kiprop_pred1[numPred1];
  real host_kiprop_pred1[numPred1];
  real donor_kiprop_pred2[numPred2];
  real host_kiprop_pred2[numPred2];
  real donor_kiprop_pred3[numPred3];
  real host_kiprop_pred3[numPred3];
  vector[numObs] log_lik;
  vector[numObs] log_lik1;
  vector[numObs] log_lik2;
  vector[numObs] log_lik3;
  vector[numObs] log_lik4;

  //parameters of interest
  real lambda_inv = 1/lambda_calc;
  real rho_inv = 1/rho;
  real delta = lambda_calc + rho;
  real delta_inv = 1/delta;

  //ODE solution for different age bins
  k_hat_pred_age1 = solve_ode_pred(ts_pred1, init_cond, parms, rdata, tb_pred1, numPred1, tb_time_pred1);
  k_hat_pred_age2 = solve_ode_pred(ts_pred2, init_cond, parms, rdata, tb_pred2, numPred2, tb_time_pred2);
  k_hat_pred_age3 = solve_ode_pred(ts_pred3, init_cond, parms, rdata, tb_pred3, numPred3, tb_time_pred3);

  // Total cell counts (donor + host)
  for (i in 1:numPred1){
    y1_mean_pred_age1[i] = k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2] + k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4];
    countspred_age1[i] = exp(normal_rng(log(y1_mean_pred_age1[i]), sigma1));
  }

  for (i in 1: numPred2){
    y1_mean_pred_age2[i] = k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2] + k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4];
    countspred_age2[i] = exp(normal_rng(log(y1_mean_pred_age2[i]), sigma1));
  }

  for (i in 1: numPred3){
    y1_mean_pred_age3[i] = k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2] + k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4];
    countspred_age3[i] = exp(normal_rng(log(y1_mean_pred_age3[i]), sigma1));
  }

  // Initial conditions for fd and frcations of ki67host and ki67 donor
  y2_mean_pred_age1[1] = Nd_0;
  fdpred_age1[1] = Nd_0;
  y2_mean_pred_age2[1] = Nd_0;
  fdpred_age2[1] = Nd_0;
  y2_mean_pred_age3[1] = Nd_0;
  fdpred_age3[1] = Nd_0;

  y3_mean_pred1[1] = eps_donor;
  donor_kiprop_pred1[1] = expit(normal_rng(logit(y3_mean_pred1[1]), sigma3));
  y3_mean_pred2[1] = eps_donor;
  donor_kiprop_pred2[1] = expit(normal_rng(logit(y3_mean_pred2[1]), sigma3));
  y3_mean_pred3[1] = eps_donor;
  donor_kiprop_pred3[1] = expit(normal_rng(logit(y3_mean_pred3[1]), sigma3));

  y4_mean_pred1[1] = kappa_0;
  host_kiprop_pred1[1] = expit(normal_rng(logit(y4_mean_pred1[1]), sigma4));
  y4_mean_pred2[1] = kappa_0;
  host_kiprop_pred2[1] = expit(normal_rng(logit(y4_mean_pred2[1]), sigma4));
  y4_mean_pred3[1] = kappa_0;
  host_kiprop_pred3[1] = expit(normal_rng(logit(y4_mean_pred3[1]), sigma4));

  for (i in 2:numPred1){
    y2_mean_pred_age1[i] = (k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2])/(y1_mean_pred_age1[i] * Chi_T1(ts_pred1[i] - tb_pred1[i]));
    fdpred_age1[i] = (normal_rng((y2_mean_pred_age1[i]), sigma2));

    y3_mean_pred1[i] = k_hat_pred_age1[i, 1]/(k_hat_pred_age1[i, 1] + k_hat_pred_age1[i, 2]);
    donor_kiprop_pred1[i] =  expit(normal_rng(logit(y3_mean_pred1[i]), sigma3));

    y4_mean_pred1[i] = k_hat_pred_age1[i, 3]/(k_hat_pred_age1[i, 3] + k_hat_pred_age1[i, 4]);
    host_kiprop_pred1[i] = expit(normal_rng(logit(y4_mean_pred1[i]), sigma4));
  }

  for (i in 2:numPred2){
    y2_mean_pred_age2[i] = (k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2])/(y1_mean_pred_age2[i] * Chi_T1(ts_pred2[i] - tb_pred2[i]));
    fdpred_age2[i] = (normal_rng((y2_mean_pred_age2[i]), sigma2));

    y3_mean_pred2[i] = k_hat_pred_age2[i, 1]/(k_hat_pred_age2[i, 1] + k_hat_pred_age2[i, 2]);
    donor_kiprop_pred2[i] =  expit(normal_rng(logit(y3_mean_pred2[i]), sigma3));

    y4_mean_pred2[i] = k_hat_pred_age2[i, 3]/(k_hat_pred_age2[i, 3] + k_hat_pred_age2[i, 4]);
    host_kiprop_pred2[i] = expit(normal_rng(logit(y4_mean_pred2[i]), sigma4));
  }

  for (i in 2:numPred3){
    y2_mean_pred_age3[i] = (k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2])/(y1_mean_pred_age3[i] * Chi_T1(ts_pred3[i] - tb_pred3[i]));
    fdpred_age3[i] = (normal_rng((y2_mean_pred_age3[i]), sigma2));

    y3_mean_pred3[i] = k_hat_pred_age3[i, 1]/(k_hat_pred_age3[i, 1] + k_hat_pred_age3[i, 2]);
    donor_kiprop_pred3[i] =  expit(normal_rng(logit(y3_mean_pred3[i]), sigma3));

    y4_mean_pred3[i] = k_hat_pred_age3[i, 3]/(k_hat_pred_age3[i, 3] + k_hat_pred_age3[i, 4]);
    host_kiprop_pred3[i] = expit(normal_rng(logit(y4_mean_pred3[i]), sigma4));
  }


  // calculating the log predictive accuracy for each point
  for (n in 1:numObs) {
    log_lik1[n] = normal_lpdf(y1[n] | log(y1_mean[n]), sigma1);
    log_lik2[n] = normal_lpdf(y2[n] | (y2_mean[n]), sigma2);
    log_lik3[n] = normal_lpdf(y3[n] | logit(y3_mean[n]), sigma3);
    log_lik4[n] = normal_lpdf(y4[n] | logit(y4_mean[n]), sigma4);
    log_lik[n] = log_lik1[n] + log_lik2[n] + log_lik3[n] + log_lik4[n];
  }
}
