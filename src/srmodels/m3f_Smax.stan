data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
parameters{
  real log_a0;// initial productivity (on log scale)
  real<lower = 0> Smax; // rate capacity - fixed in this

  //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L] log_a; //year-to-year estimates of a
}
transformed parameters{
  
  
}  
model{
  //priors
  log_a0 ~ normal(1.5,2.5); //initial productivity - wide prior
  Smax ~ normal(pSmax_mean,pSmax_sig); //per capita capacity parameter - wide prior
     
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_a ~ normal(0,1); //half normal on variance (lower limit of zero)
   
  log_a[1] ~ normal(log_a0, sigma_a); 
  for(n in 2:N) log_a[n] ~ normal(log_a[n-1],  sigma_a); 
  
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]/Smax, sigma); 


  
  
}
 generated quantities{
     vector[N] log_lik;
     real b;
     vector[L] U_msy;
     vector[L] S_msy;
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[ii[n]] - S[n]/Smax, sigma);
   
    b = 1/Smax;
    U_msy = 1-lambert_w0(exp(1-log_a));
    S_msy = (1-lambert_w0(exp(1-log_a)))/b;
    }
