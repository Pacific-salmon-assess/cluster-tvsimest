data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
}
parameters {
  real log_a;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this
    
//variance components  
  real<lower = 0> sigma;
    
}
transformed parameters{
  real b;
    
  b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  //variance terms
  //target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)

   R_S ~ normal(log_a - S*b, sigma);
}
generated quantities{
 vector[N] log_lik;
 real S_max;
 real U_msy;
 real S_msy;
 
 for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b, sigma);
 
S_max = 1/b;
U_msy = 1-lambert_w0(exp(1-log_a));
S_msy = (1-lambert_w0(exp(1-log_a)))/b;
}
    
