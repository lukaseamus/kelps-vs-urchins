data{
  int n;
  vector[n] Chlorophyll; // NNLS estimates
  vector[n] CV; // NNLS coefficient of variation
  real Chlorophyll_mean; // Calculated mean
}

parameters{
  // True estimates for measurement error
  vector<lower=0>[n] Chlorophyll_true;
  
  // Likelihood mean and uncertainty
  real<lower=0> mu;
  real<lower=0> sigma;
}

model{
  // Priors
  mu ~ normal( Chlorophyll_mean , 0.08 * Chlorophyll_mean ) T[0,];
  sigma ~ exponential( 15 / ( Chlorophyll_mean + 0.4 ) );

  // Truncated normal likelihood
  Chlorophyll_true ~ normal( mu , sigma ) T[0,];
  
  // Normal measurement error
  Chlorophyll ~ normal( Chlorophyll_true , CV .* Chlorophyll_true );
}