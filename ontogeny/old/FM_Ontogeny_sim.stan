functions{
  real theta_spline(real t, real theta0, real r1, real n1) {

    real theta = exp(theta0) * (1 + t^n1 * exp(-r1 * t));
    return theta;
   }

  real theta_Piecewise(real t, real theta0, real r1, real r3, real alpha) {
     int ta = 8;
     real theta;
     real theta_ta = exp(theta0) * exp(r1 * ta);

     if (t <= ta){
       theta = exp(theta0) * exp(r1 * t);
     } else {
       theta = (theta_ta * alpha * exp(-r1 * (t-ta))) + (theta_ta * (1 - alpha) * exp(-r3 * (t-ta)));
     }
     return theta;
   }


  real[] tdt_pw(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

     real psi     = parms[1];
     real rho     = parms[2];
     real delta0  = parms[3];
     real Beta    = parms[4];
     real r_d     = parms[5];

     real theta0  = rdata[1];
     real r1      = rdata[2];
     real r3      = rdata[3];
     real alpha   = rdata[4];
     real eps     = rdata[5];

     int tb  = 40;
     real dkdt[2];

     dkdt[1] = (psi * theta_Piecewise(t, theta0, r1, r3, alpha) * eps) + rho * (2 * k[2] + k[1]) - ((1/Beta) + delta0 * exp(-r_d * (t-tb))) * k[1];
     dkdt[2] = (psi * theta_Piecewise(t, theta0, r1, r3, alpha) * (1-eps)) + (1/Beta) * k[1] - (rho + delta0 * exp(-r_d * (t-tb))) * k[2];

     return dkdt;
   }

   real[] tdt_exp(real t,  real[] k, real[] parms, real[] rdata, int[] idata) {

      real psi     = parms[1];
      real rho     = parms[2];
      real delta0 = parms[3];
      real Beta    = parms[4];
      real r_d     = parms[5];

      real theta0  = rdata[1];
      real r1      = rdata[2];
      real n1      = rdata[3];
      real eps     = rdata[4];

      int tb  = 40;
      real dkdt[2];

      dkdt[1] = (psi * theta_spline(t, theta0, r1, n1) * eps) + rho * (2 * k[2] + k[1]) - ((1/Beta) + delta0 * exp(-r_d * (t-tb))) * k[1];
      dkdt[2] = (psi * theta_spline(t, theta0, r1, n1) * (1-eps)) + (1/Beta) * k[1] - (rho + delta0 * exp(-r_d * (t-tb))) * k[2];

      return dkdt;
    }

 real[,] solve_ode_pw(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int num_index){
  real y_hat[num_index, 2];
  int idata[0];

  y_hat[1, ] = init_cond;
  y_hat[2:num_index, ] = integrate_ode_rk45(tdt_pw, init_cond, solve_time[1], solve_time[2:num_index], parms, rdata, idata);

  return y_hat;
 }

 real[,] solve_ode_exp(real[] solve_time, real[] init_cond, real[] parms, real[] rdata, int num_index){
  real y_hat[num_index, 2];
  int idata[0];

  y_hat[1, ] = init_cond;
  y_hat[2:num_index, ] = integrate_ode_rk45(tdt_exp, init_cond, solve_time[1], solve_time[2:num_index], parms, rdata, idata);

  return y_hat;
 }
}
