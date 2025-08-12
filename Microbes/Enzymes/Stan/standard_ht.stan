data{ 
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
}

parameters{
  real<lower=0> F0;
  real<lower=0> Fmax;
  real<lower=0> beta;
  real<lower=0> sigma;
}

model{
  // Priors
  F0 ~ exponential( 0.01 );
  Fmax ~ gamma( square(8e4) / square(4e4) , 8e4 / square(4e4) );
  beta ~ gamma( square(600) / square(300) , 600 / square(300) );
  sigma ~ exponential( 1e-3 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = F0 + Fmax * tanh( beta * Concentration[i] / Fmax );
  }

  // Truncated normal likelihood
  Fluorescence ~ normal( mu , sigma ) T[0,];
}
      
generated quantities{
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = F0 + Fmax * tanh( beta * Concentration[i] / Fmax );
  }
      
  vector[n] log_lik;
  for ( i in 1:n ) {
    log_lik[i] = normal_lpdf( Fluorescence[i] | mu[i] , sigma ) - 
                 normal_lccdf( 0 | mu[i] , sigma );
  }
}