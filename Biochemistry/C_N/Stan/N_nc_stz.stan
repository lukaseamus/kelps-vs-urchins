data{ 
  int n;
  vector[n] N;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Season;
  int n_Season;
  array[n] int Individual;
  int n_Individual;
}

parameters{
  // Intercepts in log space
  vector[n_Treatment] alpha_t; 
  array[n_Treatment] sum_to_zero_vector[n_Season] z_s;
  array[n_Treatment] sum_to_zero_vector[n_Individual] z_i;
      
  // Seasonal and inter-individual variability
  vector<lower=0>[n_Treatment] sigma_s;
  vector<lower=0>[n_Treatment] sigma_i;
      
  // Likelihood uncertainty
  real<lower=0> sigma; 
}

model{
  // Hyperpriors
  sigma_s ~ exponential( 5 );
  sigma_i ~ exponential( 5 );
      
  // Priors
  alpha_t ~ normal( log(2) , 0.5 );
  for (i in 1:n_Treatment) {
    z_s[i][] ~ normal( 0 , 1 );
    z_i[i][] ~ normal( 0 , 1 );
  }
  sigma ~ exponential( 5 );
      
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
  N ~ gamma( square( mu ) / square( sigma ) ,
             mu / square( sigma ) );
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