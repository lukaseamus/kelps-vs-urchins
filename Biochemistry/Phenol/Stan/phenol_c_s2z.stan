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
  array[n_Treatment] sum_to_zero_vector[n_Season] alpha_s;
  array[n_Treatment] sum_to_zero_vector[n_Individual] alpha_i; 
      
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
    alpha_s[i][] ~ normal( 0 , sigma_s[i] * 
                               sqrt( n_Season * inv( n_Season - 1 ) ) );
    target += log(sigma_s[i]); // Note this addition
    alpha_i[i][] ~ normal( 0 , sigma_i[i] * 
                               sqrt( n_Individual * inv( n_Individual - 1 ) ) );
    target += log(sigma_i[i]);
  }
  theta ~ exponential( 5 );
      
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( alpha_t[ Treatment[i] ] + 
                 alpha_s[ Treatment[i] ][ Season[i] ] + 
                 alpha_i[ Treatment[i] ][ Individual[i] ] );
  }

  // Gamma likelihood
  Concentration ~ gamma( mu ./ theta[Treatment] , 
                         1 ./ theta[Treatment] );
      
  // Normal measurement error
  Concentration_mean ~ normal( Concentration , Concentration_sd );
}