data{ 
  int n;
  vector[n] Concentration_mean;
  vector[n] Concentration_sd;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Season;
  int n_Season;
  array[n] int Individual;
  int n_Individual;
}

parameters{
  // True estimates for measurement error
  vector<lower=0>[n] Concentration; 
      
  // Intercepts in log space
  vector[n_Treatment] alpha_t; 
  array[n_Treatment] sum_to_zero_vector[n_Season] z_s;
  array[n_Treatment] sum_to_zero_vector[n_Individual] z_i;
      
  // Seasonal and inter-individual variability
  vector<lower=0>[n_Treatment] sigma_s;
  vector<lower=0>[n_Treatment] sigma_i;
      
  // Likelihood uncertainty (scale)
  vector<lower=0>[n_Treatment] theta; 
}

model{
  // Hyperpriors
  sigma_s ~ exponential( 5 );
  sigma_i ~ exponential( 5 );
      
  // Priors
  alpha_t ~ normal( log(1.07) , 0.4 );
  for (i in 1:n_Treatment) {
    z_s[i][] ~ normal( 0 , 1 );
    z_i[i][] ~ normal( 0 , 1 );
  }
  theta ~ exponential( 5 );
      
  // Convert z-scores
  matrix[n_Treatment, n_Season] alpha_s;
  matrix[n_Treatment, n_Individual] alpha_i;
      
  for (i in 1:n_Treatment) {
    for (j in 1:n_Season) {
      alpha_s[i, j] = z_s[i][j] * sqrt( n_Season * inv( n_Season - 1 ) ) * 
                      sigma_s[i] + 0;
    }
    for (j in 1:n_Individual) {
      alpha_i[i, j] = z_i[i][j] * sqrt( n_Individual * inv( n_Individual - 1 ) ) * 
                      sigma_i[i] + 0;
    }
  }
      
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( alpha_t[ Treatment[i] ] + 
                 alpha_s[ Treatment[i] , Season[i] ] + 
                 alpha_i[ Treatment[i] , Individual[i] ] );
  }

  // Gamma likelihood
  Concentration ~ gamma( mu ./ theta[Treatment] , 
                         1 ./ theta[Treatment] );
      
  // Normal measurement error
  Concentration_mean ~ normal( Concentration , Concentration_sd );
}
    
generated quantities{
  // Save converted z-scores
  matrix[n_Treatment, n_Season] alpha_s;
  matrix[n_Treatment, n_Individual] alpha_i;
      
  for (i in 1:n_Treatment) {
    for (j in 1:n_Season) {
      alpha_s[i, j] = z_s[i][j] * sqrt( n_Season * inv( n_Season - 1 ) ) * 
                      sigma_s[i] + 0;
    }
    for (j in 1:n_Individual) {
      alpha_i[i, j] = z_i[i][j] * sqrt( n_Individual * inv( n_Individual - 1 ) ) * 
                      sigma_i[i] + 0;
    }
  }
}