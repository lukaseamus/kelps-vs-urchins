data{
  int n;
  vector<lower=0>[n] Consumption;
  vector<lower=0>[n] Defecation;
  array[n] int Season;
  int n_Season;
}

parameters{
  // Hyperparameters
  real<lower=0> alpha_mu;
  real<lower=0> alpha_theta;
  real<lower=0, upper=1> beta_mu;
  real<lower=0> beta_nu;
  
  // Seasonal parameters
  vector<lower=0>[n_Season] alpha;
  vector<lower=0, upper=1>[n_Season] beta;
  
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  alpha_mu ~ gamma( square(1) / square(0.5) , 1 / square(0.5) );
  alpha_theta ~ exponential( 5 );
  beta_mu ~ beta( 0.352 * 20 , (1 - 0.352) * 20 );
  beta_nu ~ gamma( square(15) / square(10) , 15 / square(10) );
  
  // Seasonal priors
  alpha ~ gamma( alpha_mu / alpha_theta , 1 / alpha_theta );
  beta ~ beta( beta_mu * beta_nu , (1 - beta_mu) * beta_nu );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu = alpha[Season] + beta[Season] .* Consumption;

  // Truncated Gaussian likelihood
  Defecation ~ normal( mu , sigma ) T[0,];
}