data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
parameters {
  real<upper = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<lower = 0> Smax0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_sm;
  
  //time-varying parameters
  vector[L-1] smax_dev; //year-to-year deviations in a

}

transformed parameters{
  vector<lower=0>[L] b; //b in each year
  vector<lower=0>[L] Smax; //b in each year
  
  Smax[1] = Smax0;
  for(t in 2:L){
    Smax[t] = Smax[t-1] + smax_dev[t-1]*sigma_sm;
  } 
  b=1 ./ Smax;
}  

model{
  //priors
  log_a ~ normal(1.5,2.5); //productivity
  Smax0 ~ normal(pSmax_mean,pSmax_sig); //capacity
  
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_sm ~ normal(0,pSmax_mean/10); //half normal on variance (lower limit of zero)
  
   
  smax_dev ~ std_normal();
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma);
}
generated quantities{
     vector[N] log_lik;
     real U_msy;
     vector[L] S_msy;
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a - b[ii[n]]*S[n], sigma);
     
    for(l in 1:L){S_msy[l] = (1-lambert_w0(exp(1-log_a)))/b[l];
    }
    U_msy = 1-lambert_w0(exp(1-log_a));
}
 
