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
  matrix[n_Treatment, n_Season] alpha_s;
  matrix[n_Treatment, n_Individual] alpha_i; 
      
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
  alpha_t ~ normal( log(1.74) , 0.4 );
  for (i in 1:n_Treatment) {
    alpha_s[i,] ~ normal( 0 , sigma_s[i] );
    alpha_i[i,] ~ normal( 0 , sigma_i[i] );
  }
  theta ~ exponential( 5 );
      
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( alpha_t[ Treatment[i] ] + 
                 alpha_s[ Treatment[i] , Season[i] ] + 
                 alpha_i[ Treatment[i] , Individual[i] ] );
  }

  // Gamma likelihood
  N ~ gamma( mu ./ theta[Treatment] ,
             1 ./ theta[Treatment] );
}