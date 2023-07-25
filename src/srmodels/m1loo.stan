data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      real y_oos; //log(recruits per spawner)
      real x_oos; //spawners in time T
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
      log_a ~ normal(1.5,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
      
       sigma ~ normal(0,1); //half normal on variance (lower limit of zero)  
      
      R_S ~ normal(log_a - S*b, sigma);
    }
    generated quantities{
     real log_lik_oos;
    	log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma);
    }