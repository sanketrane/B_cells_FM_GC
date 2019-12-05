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
     real psi = parms[1];
     real rho = parms[2];
     real delta = parms[3];

     real theta0 = rdata[1];
     real nu = rdata[2];
     real chiEst = rdata[3];
     real qEst = rdata[4];
     real beta = rdata[5];
     real eps_donor = rdata[6];
     real eps_host = rdata[7];

     real tb = idata[1];

     real dkdt[4];

     // since theta_spline doesnt change with time using the mean value of source counts as spline_int

     dkdt[1] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * eps_donor) + rho * (2 * k[2] + k[1]) - (beta + delta) * k[1];
     dkdt[2] = (psi * theta_spline(t, nu, theta0) * Chi_spline(t-tb, chiEst, qEst) * (1-eps_donor)) + beta * k[1] - (rho + delta) * k[2];

     dkdt[3] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * eps_host) + rho * (2 * k[4] + k[3]) - (beta + delta) * k[3];
     dkdt[4] = (psi * theta_spline(t, nu, theta0) * (1-Chi_spline(t-tb, chiEst, qEst)) * (1-eps_host)) + beta * k[3] - (rho + delta) * k[4];

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

  real[,] solve_ode(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int[] tb, int num_index, real[] tb_time){
   real y_hat[num_index, 4];
   real y0[2, 4];
   real init_tb[4];

   y0 = solve_init(tb_time, init_cond, parms, rdata, 40, 2);

   init_tb[1] = init_cond[1];
   init_tb[2] = init_cond[1];
   init_tb[3] = y0[2, 1] + y0[2, 3];
   init_tb[4] = y0[2, 2] + y0[2, 4];

   y_hat[1] = init_tb;
   for (i in 2:num_index){
     y_hat[i] = foreach_ode(solve_time[i], tb[i], init_tb, parms, rdata, tb[i]);
   }

   return y_hat;
  }
}
