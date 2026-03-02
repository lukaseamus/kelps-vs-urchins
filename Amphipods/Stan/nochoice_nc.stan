data{ 
  int n;
  vector[n] Consumption;
  array[n] int Food;
  int n_Food;
  array[n] int Season;
  int n_Season;
  array[n] int Individual;
  int n_Individual;
}

parameters{
  // Mean
  real alpha_mu;
  real<lower=0> alpha_sigma_f;
  real<lower=0> alpha_sigma_s;
  real<lower=0> alpha_sigma_fs;
  real<lower=0> alpha_sigma_i;
  vector[n_Food] alpha_z_f; // z-scores
  vector[n_Season] alpha_z_s;
  matrix[n_Food, n_Season] alpha_z_fs;
  vector[n_Individual] alpha_z_i;
  
  // Scale
  real log_theta_mu;
  real<lower=0> log_theta_sigma;
  vector[n_Food] log_theta_z;
}

transformed parameters{
  // Convert z-scores
  vector[n_Food] alpha_f = alpha_z_f * alpha_sigma_f + alpha_mu;
  vector[n_Season] alpha_s = alpha_z_s * alpha_sigma_s + 0;
  matrix[n_Food, n_Season] alpha_fs = alpha_z_fs * alpha_sigma_fs + 0;
  vector[n_Individual] alpha_i = alpha_z_i * alpha_sigma_i + 0;
  vector[n_Food] log_theta = log_theta_z * log_theta_sigma + log_theta_mu;
}

model{
  // Mean
  alpha_mu ~ normal( log(0.5) , 2 );
  alpha_sigma_f ~ normal( 0 , 1 ) T[0,];
  alpha_sigma_s ~ normal( 0 , 1 ) T[0,];
  alpha_sigma_fs ~ normal( 0 , 1 ) T[0,];
  alpha_sigma_i ~ normal( 0 , 1 ) T[0,];
  
  alpha_z_f ~ normal( 0 , 1 );
  alpha_z_s ~ normal( 0 , 1 );
  to_vector(alpha_z_fs) ~ normal( 0 , 1 );
  alpha_z_i ~ normal( 0 , 1 );
  
  // Scale
  log_theta_mu ~ normal( log(0.5) , 1 );
  log_theta_sigma ~ normal( 0 , 1 ) T[0,];
  log_theta_z ~ normal( 0 , 1 );
  
  // Model
  vector[n] mu;
  for ( i in 1:n ) mu[i] = exp(
    alpha_f[ Food[i] ] + 
    alpha_s[ Season[i] ] + 
    alpha_fs[ Food[i] , Season[i] ] + 
    alpha_i[ Individual[i] ]
  );
  
  vector[n] theta = exp( log_theta[Food] );

  // Gamma likelihood
  Consumption ~ gamma( mu ./ theta , 1 / theta );
}