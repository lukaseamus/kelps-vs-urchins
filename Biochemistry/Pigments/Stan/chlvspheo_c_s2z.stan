data{ 
  int n;
  vector[n] Chlorophyll_mean;
  vector[n] Chlorophyll_sd;
  vector[n] Pheopigments_mean;
  vector[n] Pheopigments_sd;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Season;
  int n_Season;
  array[n] int Individual;
  int n_Individual;
}

parameters{
  // True estimates for measurement error
  vector<lower=0>[n] Chlorophyll;
  vector<lower=0>[n] Pheopigments;
  
  // Intercepts in logit space
  vector[n_Treatment] alpha_t; 
  array[n_Treatment] sum_to_zero_vector[n_Season] alpha_s;
  array[n_Treatment] sum_to_zero_vector[n_Individual] alpha_i;
      
  // Seasonal and inter-individual variability
  vector<lower=0>[n_Treatment] sigma_s;
  vector<lower=0>[n_Treatment] sigma_i;
      
  // Likelihood precision: nu = ( mu * (1 - mu) / square(sigma) - 1 )
  vector<lower=0>[n_Treatment] nu;
}

transformed parameters{
  vector<lower=0,upper=1>[n] Chl_prop;
  Chl_prop = Chlorophyll ./ ( Chlorophyll + Pheopigments );
}

model{
  // Hyperpriors
  sigma_s ~ exponential( 2 );
  sigma_i ~ exponential( 2 );
      
  // Priors
  alpha_t ~ normal( logit(0.5) , 1 );
  for (i in 1:n_Treatment) {
    alpha_s[i][] ~ normal( 0 , sigma_s[i] * 
                               sqrt( n_Season * inv( n_Season - 1 ) ) );
    target += log(sigma_s[i]);
    alpha_i[i][] ~ normal( 0 , sigma_i[i] * 
                               sqrt( n_Individual * inv( n_Individual - 1 ) ) );
    target += log(sigma_i[i]);
  }
  nu ~ gamma( square(30) / square(20) , 30 / square(20) );
  
  // Measurement error (specified like predictor error)
  Chlorophyll ~ gamma( square(0.87) / square(0.5) , 0.87 / square(0.5) );
  Chlorophyll_mean ~ normal( Chlorophyll , Chlorophyll_sd );
  Pheopigments ~ gamma( square(0.87) / square(0.5) , 0.87 / square(0.5) );
  Pheopigments_mean ~ normal( Pheopigments , Pheopigments_sd );
  
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) {
    mu[i] = inv_logit( 
              alpha_t[ Treatment[i] ] + 
              alpha_s[ Treatment[i] ][ Season[i] ] + 
              alpha_i[ Treatment[i] ][ Individual[i] ]
            );
  }

  // Beta likelihood
  Chl_prop ~ beta( mu .* nu[Treatment] , 
                   (1 - mu) .* nu[Treatment] );
}