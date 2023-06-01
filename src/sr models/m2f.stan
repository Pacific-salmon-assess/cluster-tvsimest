data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
}
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

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
    mu[t] = mu[t] + (rho^(ii[t]-ii[t-1])*epsilon[t-1]); //rho raised the power of the number of time-steps between successive productivity estimates
  }
  sigma_AR = sigma*sqrt(1-rho^2);
}
model{
  //priors
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
   
  
  //autocorrelation term
  rho ~ uniform(-1,1);
  
  R_S[1] ~ normal(mu[1], sigma);
  for(t in 2:N) R_S[t] ~ normal(mu[t], sigma_AR);
  
}
generated quantities{
  vector[N] log_lik;
  real S_max;
  real U_msy;
  real S_msy;
  
  log_lik[1] = normal_lpdf(R_S[1]|mu[1], sigma);
  for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|mu[n], sigma_AR);
   
  S_max = 1/b;
  U_msy = 1-lambert_w0(exp(1-log_a));
  S_msy = (1-lambert_w0(exp(1-log_a)))/b;
}
    
