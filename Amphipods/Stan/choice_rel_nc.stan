data{ 
  int n;
  vector[n] Proportion;
  array[n] int Season;
  int n_Season;
}

parameters{
  // Likelihood mean in logit space
  real alpha_mu;
  real<lower=0> alpha_sigma;
  vector[n_Season] alpha_z; // z-scores

  // Likelihood precision
  real<lower=0> nu;
}

transformed parameters{
  // Convert z-scores
  vector[n_Season] alpha = alpha_z * alpha_sigma + alpha_mu;
}

model{
  // Priors
  alpha_mu ~ normal( logit(0.5) , 1 );
  alpha_sigma ~ normal( 0 , 1 ) T[0,]; // half-normal priors
  alpha_z ~ normal( 0 , 1 );
  nu ~ gamma( square(20) / square(10) , 20 / square(10) );
  
  // Model with link function
  vector[n] mu = inv_logit( alpha[Season] );
  
  // Beta likelihood
  Proportion ~ beta( mu * nu , (1 - mu) * nu );
}