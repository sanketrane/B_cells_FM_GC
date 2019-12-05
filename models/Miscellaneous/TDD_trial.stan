functions{
  real theta_spline(real t,
    real nu,               // nu is the rate of decline of source influx assuming that its proportional to the counts of source compartment
    real theta0) {     // theta0 gives the initial counts of the source compartment

      real theta;
      int ta = 40;         // the earliest age at BMT

      theta = theta0 * exp(-nu * (t-ta));
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


  real[] KHM(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {
    real psi        = parms[1];
    real f_fast     = parms[2];
    real rhoFast    = parms[3];
    real rhoSlow    = parms[4];
    real Beta       = parms[5];
    real lambdaFast = parms[6];
    real lambdaSlow = parms[7];

    real theta0    = rdata[1];
    real nu        = rdata[2];
    real chiEst    = rdata[3];
    real qEst      = rdata[4];
    real eps_donor = rdata[5];
    real eps_host  = rdata[6];

    real tb = idata[1];

    real dkdt[8];

    dkdt[1] = (psi * f_fast * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rhoFast * (2 * k[2] + k[1]) - ((1/Beta) + (lambdaFast+rhoFast)) * k[1];
    dkdt[2] = (psi * f_fast * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[1] - (rhoFast + (lambdaFast+rhoFast)) * k[2];
    dkdt[3] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rhoSlow * (2 * k[4] + k[3]) - ((1/Beta) + (lambdaSlow+rhoSlow)) * k[3];
    dkdt[4] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + (1/Beta) * k[3] - (rhoSlow + (lambdaSlow+rhoSlow)) * k[4];

    dkdt[5] = (psi * f_fast * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rhoFast * (2 * k[6] + k[5]) - ((1/Beta) + (lambdaFast+rhoFast)) * k[5];
    dkdt[6] = (psi * f_fast * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[5] - (rhoFast + (lambdaFast+rhoFast)) * k[6];
    dkdt[7] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rhoSlow * (2 * k[8] + k[7]) - ((1/Beta) + (lambdaSlow+rhoSlow)) * k[7];
    dkdt[8] = (psi * (1-f_fast) * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + (1/Beta) * k[7] - (rhoSlow + (lambdaSlow+rhoSlow)) * k[8];


    return dkdt;
  }

  real[] foreach_init(real tb, real ta, real[] init_cond, real[] parms, real[] rdata, int x_i){

    return to_array_1d(integrate_ode_bdf(KHM, init_cond, ta, rep_array(tb, 1), parms, rdata, rep_array(x_i, 1)));
  }

  real[,] solve_init(real[] tb_time, real[] init_cond, real[] parms, real[] rdata, int x_i, int num_tb){
   real y_init[num_tb, 8];

   y_init[1] = init_cond;
   for (i in 2:num_tb){
     y_init[i] = foreach_init(tb_time[i], tb_time[1], init_cond, parms, rdata, x_i);
   }

   return y_init;
  }

  real[] foreach_ode(real ts, real t0, real[] init_cond, real[] parms, real[] rdata, int x_i) {

    return to_array_1d(integrate_ode_bdf(KHM, init_cond, t0, rep_array(ts, 1), parms, rdata, rep_array(x_i, 1)));
   }

   real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time, int[] tb_index, int num_tb){
    real y_hat[num_index, 8];
    real y0[num_tb, 8];
    real init_tb[8];

    y0 = solve_init(tb_time, init_cond, parms, rdata, 40, num_tb);

    init_tb[1] = init_cond[1];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[2] = init_cond[2];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[3] = init_cond[3];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[4] = init_cond[4];                                          //at tbmt - donor ki67Hi subset size is zero

    for (i in 1:num_index){
      init_tb[5] = y0[tb_index[i], 1] + y0[tb_index[i], 5];              //at tbmt - all ki67Hi cells would be host
      init_tb[6] = y0[tb_index[i], 2] + y0[tb_index[i], 6];              //at tbmt - all ki67Lo cells would be host
      init_tb[7] = y0[tb_index[i], 3] + y0[tb_index[i], 7];              //at tbmt - all ki67Hi cells would be host
      init_tb[8] = y0[tb_index[i], 4] + y0[tb_index[i], 8];              //at tbmt - all ki67Lo cells would be host
      y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
    }

    return y_hat;
   }

   real[,] solve_ode_pred(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
    real y_hat[num_index, 8];
    real y0[2, 8];
    real init_tb[8];

    y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

    init_tb[1] = init_cond[1];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[2] = init_cond[2];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[3] = init_cond[3];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[4] = init_cond[4];                                          //at tbmt - donor ki67Hi subset size is zero
    init_tb[5] = y0[2, 1] + y0[2, 5];                                   //at tbmt - all ki67Hi cells would be host
    init_tb[6] = y0[2, 2] + y0[2, 6];                                   //at tbmt - all ki67Lo cells would be host
    init_tb[7] = y0[2, 3] + y0[2, 7];                                   //at tbmt - all ki67Hi cells would be host
    init_tb[8] = y0[2, 4] + y0[2, 8];                                   //at tbmt - all ki67Lo cells would be host

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
