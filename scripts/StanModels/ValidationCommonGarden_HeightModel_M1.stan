data {
  int N;                                    // Number of individuals
  vector[N] Y;                              // Response variable (individual tree height)
  vector[N] X;                              // Genomic offset or climatic transfer distances
  int<lower=0> nb_bloc;                     // Number of blocks
  int<lower=0, upper=nb_bloc> bloc[N];      // Blocks
}

parameters {
  vector[nb_bloc] alpha_bloc;               // Intercepts of the blocks
  real beta_X1;                             // Linear coefficent of GO or CTD
  real beta_X2;                             // Quadratic coefficent of GO or CTD
  real<lower = 0>  sigma;                   // Residual variance of the model
}

transformed parameters {
  vector[N] mu;                             // Linear predictor
  real R_squared;                           // R^2 to evaluate the goodness of fit of the model
  
  mu = alpha_bloc[bloc] + beta_X1 * X + beta_X2 * square(X);
  R_squared = 1 - variance(Y - mu) / variance(Y);
}

model {
  
  // Likelihood
  Y ~ normal(mu, sigma);
  
  // Priors
  sigma ~ exponential(1);
  alpha_bloc ~ std_normal();
  beta_X1 ~ std_normal();
  beta_X2 ~ std_normal();
}

// generated quantities{
//   // log likelihood for loo
//   vector[N] log_lik;
//   vector[N] muhat;
//   for (n in 1:N) {
//     log_lik[n] = normal_lpdf(Y[n] |mu[n],sigma);
//     muhat[n] = normal_rng(mu[n], sigma);
//   }
// }
