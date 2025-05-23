data{ 
  int n;
  vector[n] C;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Season;
  int n_Season;
  array[n] int Individual;
  int n_Individual;
}
  
transformed data{
  vector[n] C_prop;
  C_prop = C / 100;
}

parameters{
  // Intercepts in logit space
  vector[n_Treatment] alpha_t; 
  matrix[n_Treatment, n_Season] alpha_s;
  matrix[n_Treatment, n_Individual] alpha_i;
      
  // Seasonal and inter-individual variability
  vector<lower=0>[n_Treatment] sigma_s;
  vector<lower=0>[n_Treatment] sigma_i;
      
  // Likelihood uncertainty
  real<lower=0> sigma;
}

model{
  // Hyperpriors
  sigma_s ~ exponential( 20 );
  sigma_i ~ exponential( 20 );
      
  // Priors
  alpha_t ~ normal( logit(0.3) , 0.3 );
  for (i in 1:n_Treatment) {
    alpha_s[i,] ~ normal( 0 , sigma_s[i] );
    alpha_i[i,] ~ normal( 0 , sigma_i[i] );
  }
  sigma ~ exponential( 20 );
      
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = inv_logit( 
              alpha_t[ Treatment[i] ] + 
              alpha_s[ Treatment[i] , Season[i] ] +
              alpha_i[ Treatment[i] , Individual[i] ]
            );
  }

  // Beta likelihood
  vector[n] nu = ( mu .* (1 - mu) / square(sigma) - 1 );
  C_prop ~ beta( mu .* nu , (1 - mu) .* nu );
}