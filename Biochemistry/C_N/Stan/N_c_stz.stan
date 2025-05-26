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
  array[n_Treatment] sum_to_zero_vector[n_Season] alpha_s;
  array[n_Treatment] sum_to_zero_vector[n_Individual] alpha_i; 
      
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
    alpha_s[i][] ~ normal( 0 , sigma_s[i] * 
                               sqrt( n_Season * inv( n_Season - 1 ) ) );
    alpha_i[i][] ~ normal( 0 , sigma_i[i] * 
                               sqrt( n_Individual * inv( n_Individual - 1 ) ) );
  }
  sigma ~ exponential( 5 );
      
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = exp( alpha_t[ Treatment[i] ] + 
                 alpha_s[ Treatment[i] ][ Season[i] ] + 
                 alpha_i[ Treatment[i] ][ Individual[i] ] );
  }

  // Gamma likelihood
  N ~ gamma( square( mu ) / square( sigma ) ,
             mu / square( sigma ) );
}