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
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
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
// Discrete state model
simplex[K] A[K]; // transition probabilities
simplex[K] pi1; // initial state probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K] log_b; // rate capacity - fixed in this
real<lower=0> sigma; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
vector[K] b;

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{

log_a ~ normal(1.5,2.5);
for(k in 1:K) log_b[k] ~ normal(logbeta_pr,logbeta_pr_sig); //capacity

  sigma ~ normal(0,1); //half normal on variance (lower limit of zero)
  pi1~ dirichlet(rep_vector(1,K));
for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet[k,]);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted

//slope based on regime in year N
real b_1b;
real b_3b;
real b_5b;

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
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma);
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
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma);
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

b_1b = b[zstar[N]]; //slope based on most probable state in sample N
b_3b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]])/3); //intercept
b_5b = exp((log_b[zstar[N]]+log_b[zstar[N-1]]+log_b[zstar[N-2]]+log_b[zstar[N-3]]+log_b[zstar[N-4]])/5); //intercept

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b_1b, sigma);
log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma);
log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma);
}