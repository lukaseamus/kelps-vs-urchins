data{
  int n;
  vector[n] Pheopigments; // NNLS estimates
  vector[n] CV; // NNLS coefficient of variation
  real Pheopigments_mean; // Calculated mean
}

parameters{
  // True estimates for measurement error
  vector<lower=0>[n] Pheopigments_true;
  
  // Likelihood mean and uncertainty
  real<lower=0> mu;
  real<lower=0> sigma;
}

model{
  // Priors
  mu ~ normal( Pheopigments_mean , 0.08 * Pheopigments_mean ) T[0,];
  sigma ~ exponential( 15 / ( Pheopigments_mean + 0.4 ) );

  // Truncated normal likelihood
  Pheopigments_true ~ normal( mu , sigma ) T[0,];
  
  // Normal measurement error
  Pheopigments ~ normal( Pheopigments_true , CV .* Pheopigments_true );
}
