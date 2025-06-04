data{
  int n;
  vector[n] Total; // NNLS estimates
  vector[n] CV; // NNLS coefficient of variation
  real Total_mean; // Calculated mean
}

parameters{
  // True estimates for measurement error
  vector<lower=0>[n] Total_true;
  
  // Likelihood mean and uncertainty
  real<lower=0> mu;
  real<lower=0> sigma;
}

model{
  // Priors
  mu ~ normal( Total_mean , 0.08 * Total_mean ) T[0,];
  sigma ~ exponential( 15 / ( Total_mean + 0.4 ) );

  // Truncated normal likelihood
  Total_true ~ normal( mu , sigma ) T[0,];
  
  // Normal measurement error
  Total ~ normal( Total_true , CV .* Total_true );
}