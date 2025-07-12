data{
  int n;
  vector[n] Speed_c;
  vector<lower=0>[n] Depth;
}

parameters{
  // Parameters describing mu
  real alpha;
  real beta;

  // Likelihood uncertainty parameter
  real<lower=0> theta;
}

model{
  // Priors
  alpha ~ normal( log(50) , 1 );
  beta ~ normal( -1 , 0.5 );
  theta ~ normal( 2 , 5 ) T[0,];
  
  // Model
  vector[n] mu = exp( alpha + beta * Speed_c );
  
  // Gamma likelihood
  Depth ~ gamma( mu / theta , 1 / theta );
}