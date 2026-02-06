data{ 
  int n;
  vector[n] Proportion;
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
  // Intercepts in logit space
  real alpha;
  real<lower=0> sigma_f;
  real<lower=0> sigma_g;
  real<lower=0> sigma_s;
  real<lower=0> sigma_r;
  
  vector[n_Function] f;
  vector[n_Genus] g;
  vector[n_Species] s;
  vector[n_Reference] r;

  // Likelihood precision
  real<lower=0> nu;
}

model{
  // Priors
  alpha ~ normal( logit(0.5) , 0.8 );
  sigma_f ~ normal( 0 , 0.5 )T[0,]; // half-normal priors
  sigma_g ~ normal( 0 , 0.5 )T[0,];
  sigma_s ~ normal( 0 , 0.5 )T[0,];
  sigma_r ~ normal( 0 , 0.5 )T[0,];
  
  f ~ normal( alpha , sigma_f );
  g ~ normal( 0 , sigma_g );
  s ~ normal( 0 , sigma_s );
  r ~ normal( 0 , sigma_r );
  
  nu ~ gamma( square(30) / square(20) , 30 / square(20) );
  
  // Model with link function
  vector[n] mu = inv_logit(
    f[Function] + g[Genus] + s[Species] + r[Reference]
  );
  
  // Beta likelihood
  Proportion ~ beta( mu * nu , (1 - mu) * nu );
}