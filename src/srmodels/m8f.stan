functions {
      vector normalize(vector x) {
        return x / sum(x);
      }
    }
    data {
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      int<lower=1> K; //number of hidden regime states
      matrix[K,K] alpha_dirichlet; //prior inputs for dirichlet 
    }
    parameters {
      // Discrete state model
      simplex[K] A[K]; // transition probabilities
      
      // A[i][j] = p(z_t = j | z_{t-1} = i)
      // Continuous observation model
      ordered[K] log_a; // regime max. productivity
      vector[K] log_b; // regime rate capacity 
      real<lower=0> sigma; // observation standard deviations
    }
    
    transformed parameters {
	 simplex[K] pi1; // initial state probabilities
      vector[K] logalpha[N];
      vector[K] b; //
        
	  pi1=rep_vector(1.0/K,K);

        b=exp(log_b);
        
        { // Forward algorithm log p(z_t = j | y_{1:t})
          real accumulator[K];
          
          logalpha[1] = log(pi1) + normal_lpdf(R_S[1]|log_a - b*S[1], sigma);
          for (t in 2:N) {
            for (j in 1:K) { // j = current (t)
            for (i in 1:K) { // i = previous (t-1)
            // Murphy (2012) p. 609 eq. 17.48
            // belief state + transition prob + local evidence at t
            accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma);
            }
            logalpha[t, j] = log_sum_exp(accumulator);
            }
          }
        } // Forward
    }
    model{
     
      log_a ~ normal(1.5,2.5);
      log_b ~ normal(-12,3);

      sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
      
      for(k in 1:K){
        A[k,] ~ dirichlet(alpha_dirichlet[k,]);
      }
      
      target += log_sum_exp(logalpha[N]);
    }
generated quantities {
vector[N] log_lik;
//HMM estimators
int<lower=1, upper=K> zstar[N]; //most-likely regime state sequence
real logp_zstar;
vector[K] alpha[N]; //forward state probabilities
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N]; //backward state probabilities
vector[K] gamma[N]; //forward-backward state probabilities

//reference points
vector[K] S_max;
vector[K] U_msy;
vector[K] S_msy;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward


{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
}

for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|log_a[zstar[n]] - S[n]*b[zstar[n]], sigma);

for(k in 1:K){
S_max[k] = 1/b[k];
U_msy[k] = 1-lambert_w0(exp(1-log_a[k]));
S_msy[k] = (1-lambert_w0(exp(1-log_a[k])))/b[k];
}

}

