data{ 
  int n;
  vector[n] Distance;
  array[n] int Tissue;
  int n_Tissue;
}

transformed data{
  vector[n] Distance_nz = Distance + 1e-3; // Add 1 m
}

parameters{
  // Intercepts in log space
  vector[n_Tissue] alpha;

  // Likelihood uncertainty (scale)
  vector<lower=0>[n_Tissue] theta; 
}

model{
  // Priors
  alpha ~ normal( log(100) , 1 );
  theta ~ exponential( 1 );
  
  // Model with link function
  vector[n] mu = exp( alpha[Tissue] );

  // Gamma likelihood
  Distance_nz ~ gamma( mu ./ theta[Tissue] ,
                       1 ./ theta[Tissue] );
}