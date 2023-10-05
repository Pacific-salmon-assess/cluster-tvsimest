data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
  real psig_b;
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<lower = 0> Smax0; // rate capacity -

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] smax_dev; //year-to-year deviations in smax

}

transformed parameters{
  vector<lower=0>[L] Smax; //year-to-year deviations in a
  
  Smax[1]=Smax0;
  
  for(t in 2:L){
    Smax[t]=Smax[t-1] + smax_dev[t-1]*sigma_b;
  }

  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(1.5,2.5); //productivity
  Smax0 ~ normal(pSmax_mean,pSmax_sig); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,psig_b); //half normal on variance (lower limit of zero)
  
  smax_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a-S[n]/Smax[ii[n]], sigma);
}
generated quantities{
  real Smax_3b;
  real b_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  Smax_3b = (Smax[ii[N]]+Smax[ii[N-1]]+Smax[ii[N-2]])/3;
  Smax_5b = (Smax[ii[N]]+Smax[ii[N-1]]+Smax[ii[N-2]]+Smax[ii[N-3]]+Smax[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos/Smax[ii[N]], sigma);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos/Smax_3b, sigma);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos/Smax_5b, sigma);
 }