functions{
  real[] normalskew_rng(int n, real mu, real sigma,  real alpha) {

      real dist[n];
      for (i in 1:n){
        dist[i] = skew_normal_rng(mu, sigma, alpha);
      }
    return dist;
   }
}
