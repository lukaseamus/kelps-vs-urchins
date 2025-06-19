data{ 
  int n;
  vector<lower=0>[n] Consumption;
  array[n] int Season;
  int n_Season;
}

parameters{
  // Hyperparameters
  real alpha_mu; 
  real<lower=0> alpha_sigma;
      
  // Seasonal parameters
  vector[n_Season] alpha; 

  // Likelihood uncertainty (scale)
  real<lower=0> theta; 
}

model{
  // Hyperpriors
  alpha_mu ~ normal( log(15) , 1 );
  alpha_sigma ~ exponential( 2 );
      
  // Seasonal priors
  alpha ~ normal( alpha_mu , alpha_sigma );
  
  // Likelihood uncertainty prior
  theta ~ exponential( 1 );
      
  // Model with link function
  vector[n] mu = exp( alpha[Season] );

  // Gamma likelihood
  Consumption ~ gamma( mu / theta , 1 / theta );
}