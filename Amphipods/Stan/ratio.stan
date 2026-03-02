data{
  int n;
  vector[n] Buoyant;
  vector[n] Dry;
  array[n] int Experiment;
  int n_Experiment;
}

parameters{
  real logit_beta_mu; // Hyperparameters
  real<lower=0> logit_beta_sigma;
  real log_sigma_mu;
  real<lower=0> log_sigma_sigma;
  vector[n_Experiment] logit_beta_z; // z-scores
  vector[n_Experiment] log_sigma_z;
}

transformed parameters{
  // Convert z-scores
  vector[n_Experiment] logit_beta = logit_beta_z * logit_beta_sigma + logit_beta_mu;
  vector[n_Experiment] log_sigma = log_sigma_z * log_sigma_sigma + log_sigma_mu;
}

model{
  // Priors
  logit_beta_mu ~ normal( logit(0.1) , 0.5 );
  logit_beta_sigma ~ normal( 0 , 0.5 ) T[0,]; // half-normal prior
  log_sigma_mu ~ normal( log(1) , 0.5 );
  log_sigma_sigma ~ normal( 0 , 0.5 ) T[0,];
  logit_beta_z ~ normal( 0 , 1 ); // standard normal priors for z-scores
  log_sigma_z ~ normal( 0 , 1 );
  
  // Model
  vector[n] beta = inv_logit( logit_beta[Experiment] );
  vector[n] mu = beta .* Buoyant;
  vector[n] sigma = exp( log_sigma[Experiment] );
  
  // Truncated normal likelihood
  Dry ~ normal( mu , sigma ) T[0,];
}