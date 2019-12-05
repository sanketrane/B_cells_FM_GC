functions{
  real Chi_spline( real t,
      real chiEst,
      real qEst) {

        real chi;

        chi = chiEst * (1 - exp(-qEst * t));
        return chi;
    }
}

data{
int<lower = 1> numObs;
int<lower = 0>  Time[numObs];
int<lower  = 1> numPred;
int<lower  = 0> ts_pred[numPred];
real<lower = 0> chiT2[numObs];
}

transformed data{
real y[numObs];

for (i in 1:numObs){
y[i] = chiT2[i];
}
}

parameters{
  real<lower = 0, upper = 1> chiEst;
  real<lower =0> qEst;
  real<lower = 0> sigma;
}

transformed parameters{
 real ymean[numObs];

 for(i in 1:numObs){
 ymean[i] = Chi_spline(Time[i], chiEst, qEst);
 }
}

model{
  chiEst ~ uniform(0, 1);
  qEst ~ normal(0.01, 2);
  sigma ~ normal(0, 1);

  y ~ normal(ymean, sigma);
}

generated quantities{
  real ymean_pred[numPred];
  real chipred[numPred];

  for(i in 1:numPred){
  ymean_pred[i] = Chi_spline(ts_pred[i], chiEst, qEst);
  }

for(i in 1:numPred){
    chipred[i] = normal_rng(ymean_pred[i], sigma);
    }
}
