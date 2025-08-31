data{ 
  int n;
  vector<lower=0>[n] Fluorescence;
  vector<lower=0>[n] Concentration;
  array[n] int Group;
  int n_Group;
}

parameters{
  vector<lower=0>[n_Group] F0;
  vector<lower=0>[n_Group] Fmax;
  vector<lower=0>[n_Group] beta;
  real<lower=0> sigma;
}

model{
  // Priors
  F0 ~ exponential( 0.01 );
  Fmax ~ gamma( square(15e4) / square(8e4) , 15e4 / square(8e4) );
  beta ~ gamma( square(1e3) / square(500) , 1e3 / square(500) );
  sigma ~ exponential( 1e-4 );

  // Model
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = F0[Group[i]] + Fmax[Group[i]] * 
            ( 1 - exp( -beta[Group[i]] * Concentration[i] / Fmax[Group[i]] ) );
  }

  // Truncated normal likelihood
  Fluorescence ~ normal( mu , sigma ) T[0,];
}
      
generated quantities{
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = F0[Group[i]] + Fmax[Group[i]] *
            ( 1 - exp( -beta[Group[i]] * Concentration[i] / Fmax[Group[i]] ) );
  }

  vector[n] log_lik;
  for ( i in 1:n ) {
    log_lik[i] = normal_lpdf( Fluorescence[i] | mu[i] , sigma ) -
                 normal_lccdf( 0 | mu[i] , sigma );
  }
}