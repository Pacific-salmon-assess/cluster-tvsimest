data{
  int<lower=1> N;//number of annual samples (time-series length)
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T

 }
parameters{
  real log_a;// initial productivity (on log scale)
  real<upper=0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = -1, upper = 1> rho;

}
transformed parameters{
real b;
vector[N] mu;
vector[N] epsilon; //residuals
real sigma_AR;

b = exp(log_b);
mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]);
  }

sigma_AR = sigma*sqrt(1-rho^2);
}
model{
  //priors
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
  //variance terms
   target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero   
  
  //autocorrelation term
  rho ~ uniform(-1,1);

 R_S[1] ~ normal(mu[1], sigma);
 
 for(t in 2:N) R_S[t] ~ normal(mu[t], sigma_AR);
}
generated quantities{
  real log_lik_oos;

  log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b+rho*epsilon[N], sigma_AR);
 } 