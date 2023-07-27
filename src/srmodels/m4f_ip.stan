data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real pSmax_mean;
  real pSmax_sig;
}
transformed data{
real logbeta_pr;
real logbeta_pr_sig;

logbeta_pr_sig=sqrt(log(1+((1/pSmax_sig)*(1/pSmax_sig))/((1/pSmax_mean)*(1/pSmax_mean)))); //this converts sigma on the untransformed scale to a log scale
logbeta_pr=log(1/pSmax_mean)-0.5*logbeta_pr_sig*logbeta_pr_sig; //convert smax prior to per capita slope - transform to log scale with bias correction

}
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(1.5,2.5); //productivity
  b0 ~ normal(logbeta_pr,logbeta_pr_sig); //capacity
  
  //variance terms
  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  sigma_b ~ normal(0,1); //half normal on variance (lower limit of zero)
  
   
  b_dev ~ std_normal();
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma);
}
generated quantities{
     vector[N] log_lik;
     vector[L] S_max;
     real U_msy;
     vector[L] S_msy;
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a - b[ii[n]]*S[n], sigma);
     
    for(l in 1:L){ S_max[l] = 1/b[l];
                   S_msy[l] = (1-lambert_w0(exp(1-log_a)))/b[l];
    }
    U_msy = 1-lambert_w0(exp(1-log_a));
}
 
