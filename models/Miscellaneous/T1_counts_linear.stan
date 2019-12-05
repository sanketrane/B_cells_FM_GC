functions{

  real theta_spline(real t,
    real nu,
    real theta0) {

      real theta;

      theta = theta0 * exp(-nu * t);
      return theta;
    }
  }


  data{
  int<lower  = 1> numObs;
  int<lower  = 0> Time[numObs];
  int<lower  = 1> numPred;
  int<lower  = 0> ts_pred[numPred];
  real<lower = 0> counts[numObs];
}

transformed data{
  real y[numObs];

  for (i in 1:numObs){
    y[i] = log(counts[i]);
  }
}

parameters{
  real theta0Log;
  real nu;
  real<lower = 0> sigma;
}

transformed parameters{
  real ymean[numObs];
  real theta0 = exp(theta0Log);

 for(i in 1:numObs){
   ymean[i] = theta_spline(Time[i], nu, theta0);
 }
}

model{
  nu ~ normal(0.001, 0.05);
  theta0Log ~ normal(16, 2);
  sigma ~ cauchy(0.2, 2);

  y ~ normal(log(ymean), sigma);
}

generated quantities{
  real ymean_pred[numPred];
  real countspred[numPred];

  for(i in 1:numPred){
    ymean_pred[i] = theta_spline(ts_pred[i], nu, theta0);
  }

  for(i in 1:numPred){
    countspred[i] = normal_rng((ymean_pred[i]), sigma);
  }
}
