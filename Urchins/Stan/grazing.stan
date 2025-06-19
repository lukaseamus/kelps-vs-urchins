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
  vector<lower=0, upper=1>[n_Season] beta;

  // Likelihood uncertainty parameter
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  beta_mu ~ beta( 0.44 * 10 , (1 - 0.44) * 10 );
  beta_nu ~ gamma( square(10) / square(7) , 10 / square(7) );
  
  // Seasonal priors
  beta ~ beta( beta_mu * beta_nu , (1 - beta_mu) * beta_nu );
  
  // Likelihood uncertainty
  sigma ~ exponential( 1 );
  
  // Model
  vector[n] mu = beta[Season] .* Consumption;
  
  // Truncated Gaussian likelihood
  Defecation ~ normal( mu , sigma ) T[0,];
}