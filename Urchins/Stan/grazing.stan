data{
  int n;
  vector<lower=0>[n] Consumption;
  vector<lower=0>[n] Defecation;
  array[n] int Season;
  int n_Season;
}

parameters{
  // Hyperparameters
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
  beta_mu ~ beta( 0.44 * 30 , (1 - 0.44) * 30 );
  beta_nu ~ gamma( square(30) / square(20) , 30 / square(20) );
  
  // Seasonal priors
  beta ~ beta( beta_mu * beta_nu , (1 - beta_mu) * beta_nu );
  
  // Likelihood uncertainty prior
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu = beta[Season] .* Consumption;

  // Truncated Gaussian likelihood
  Defecation ~ normal( mu , sigma ) T[0,];
}