data{
  int n;
  vector<lower=0,upper=2.8>[n] Absorbance;
  real Absorbance_mean;
}

parameters{
  real<lower=0,upper=2.8> mu;
  real<lower=0> sigma;
}

model{
  // Priors
  mu ~ normal( Absorbance_mean , 0.08 * Absorbance_mean ) T[ 0 , 2.8 ];
  sigma ~ exponential( 15 / ( Absorbance_mean + 0.4 ) );

  // Likelihood
  Absorbance ~ normal( mu , sigma ) T[ 0 , 2.8 ];
}