functions{
  // Beta prime log probability density function
  real betap_lpdf( real y , real alpha , real beta ) {
    return ( alpha - 1 ) * log( y )
    - ( alpha + beta ) * log1p( y ) -
    lbeta( alpha , beta );
  }
}

data{ 
  int n;
  vector[n] Ratio;
  array[n] int Compound;
  int n_Compound;
  array[n] int Function;
  int n_Function;
  array[n] int Genus;
  int n_Genus;
  array[n] int Species;
  int n_Species;
  array[n] int Reference;
  int n_Reference;
}

parameters{
  // Global parameters
  real alpha;
  real<lower=0> sigma_c;

  real log_nu_mu;
  real<lower=0> log_nu_sigma;
  
  // Compound parameters
  vector[n_Compound] z_c;
  vector[n_Compound] log_nu_z;
  
  vector<lower=0>[n_Compound] sigma_f;
  vector<lower=0>[n_Compound] sigma_g;
  vector<lower=0>[n_Compound] sigma_s;
  vector<lower=0>[n_Compound] sigma_r;
  
  // Parameters nested in compound
  matrix[n_Compound, n_Function] z_f;
  matrix[n_Compound, n_Genus] z_g;
  matrix[n_Compound, n_Species] z_s;
  matrix[n_Compound, n_Reference] z_r;
}

transformed parameters{
  // Convert z-scores
  vector[n_Compound] c = z_c * sigma_c + alpha;
  vector[n_Compound] log_nu = log_nu_z * log_nu_sigma + log_nu_mu;
  
  matrix[n_Compound, n_Function] f;
  matrix[n_Compound, n_Genus] g;
  matrix[n_Compound, n_Species] s;
  matrix[n_Compound, n_Reference] r;
  
  for ( i in 1:n_Compound ) {
    for ( j in 1:n_Function ) {
      f[i,j] = z_f[i,j] * sigma_f[i] + 0;
    }
    for ( j in 1:n_Genus ) {
      g[i,j] = z_g[i,j] * sigma_g[i] + 0;
    }
    for ( j in 1:n_Species ) {
      s[i,j] = z_s[i,j] * sigma_s[i] + 0;
    }
    for ( j in 1:n_Reference ) {
      r[i,j] = z_r[i,j] * sigma_r[i] + 0;
    }
  }
}

model{
  // Priors
  alpha ~ normal( log(1) , 0.5 );
  sigma_c ~ normal( 0 , 0.5 )T[0,];
  
  log_nu_mu ~ normal( log(20) , 0.5 );
  log_nu_sigma ~ normal( 0 , 0.5 )T[0,];
  
  z_c ~ normal( 0 , 1 );
  log_nu_z ~ normal( 0 , 1 );
  
  sigma_f ~ normal( 0 , 0.5 )T[0,];
  sigma_g ~ normal( 0 , 0.5 )T[0,];
  sigma_s ~ normal( 0 , 0.5 )T[0,];
  sigma_r ~ normal( 0 , 0.5 )T[0,];
  
  to_vector(z_f) ~ normal( 0 , 1 );
  to_vector(z_g) ~ normal( 0 , 1 );
  to_vector(z_s) ~ normal( 0 , 1 );
  to_vector(z_r) ~ normal( 0 , 1 );
  
  // Model with link function
  vector[n] mu;
  for ( i in 1:n ) mu[i] = exp(
    c[Compound[i]] + 
    f[Compound[i], Function[i]] + 
    g[Compound[i], Genus[i]] + 
    s[Compound[i], Species[i]] + 
    r[Compound[i], Reference[i]]
  );
  
  vector[n] nu = exp( log_nu[Compound] );
  
  // Beta prime likelihood
  for ( i in 1:n ) Ratio[i] ~ betap( mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
}