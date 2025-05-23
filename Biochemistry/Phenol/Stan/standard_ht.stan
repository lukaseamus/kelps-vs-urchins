data{ 
  int n;
  vector<lower=0>[n] Absorbance;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> A0;
  real<lower=0> Amax;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  A0 ~ exponential( 10 );
  Amax ~ gamma( square(3) / square(2) , 3 / square(2) );
  beta ~ gamma( square(2.5) / square(2) , 2.5 / square(2) );
  sigma ~ exponential( 5 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = A0 + Amax * tanh( beta * Concentration[i] / Amax );
  }

  // Truncated normal likelihood
  Absorbance ~ normal( mu , sigma ) T[0,];
}
      
generated quantities{
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = A0 + Amax * tanh( beta * Concentration[i] / Amax );
  }
      
  vector[n] log_lik;
  for ( i in 1:n ) {
    log_lik[i] = normal_lpdf( Absorbance[i] | mu[i] , sigma ) - 
                 normal_lccdf( 0 | mu[i] , sigma );
  }
}