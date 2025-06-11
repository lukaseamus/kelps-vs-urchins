data{ 
  int n;
  vector[n] Speed;
  array[n] int Tissue;
  int n_Tissue;
}

parameters{
  // Intercepts in log space
  vector[n_Tissue] alpha;

  // Likelihood uncertainty (scale)
  vector<lower=0>[n_Tissue] theta; 
}

model{
  // Priors
  alpha ~ normal( log(0.1) , 1.2 );
  theta ~ exponential( 100 );    
      
  // Model with link function
  vector[n] mu = exp( alpha[Tissue] );

  // Gamma likelihood
  Speed ~ gamma( mu ./ theta[Tissue] ,
                 1 ./ theta[Tissue] );
}