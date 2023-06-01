data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(1.5,2.5); //initial productivity
  log_b0 ~ normal(-12,3); //initial capacity
  
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_a ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,1); //half normal on variance (lower limit of zero)
  
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma);
}
 generated quantities{
     vector[N] log_lik;
     vector[L] S_max;
     vector[L] U_msy;
     vector[L] S_msy;
     
   for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[ii[n]] - S[n]*b[ii[n]], sigma);
   
   for(l in 1:L){ S_max[l] = 1/b[l];
    U_msy[l] = 1-lambert_w0(exp(1-log_a[l]));
    S_msy[l] = (1-lambert_w0(exp(1-log_a[l])))/b[l];
   }
}
