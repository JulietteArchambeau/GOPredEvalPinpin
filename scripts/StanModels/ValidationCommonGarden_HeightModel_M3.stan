data {
  int N;
  vector[N] Y;                              // Response variable (individual height)
  vector[N] X;                              // Genomic offset or climatic transfer distances
  vector[N] H;                              // Initial height of the populations
  int<lower=0> nb_bloc;                     // Number of blocks
  int<lower=0, upper=nb_bloc> bloc[N];      // Blocks
  int<lower=0> nb_clon;                     // Number of clones
  int<lower=0, upper=nb_clon> clon[N];      // Clones
}

parameters {
  vector[nb_bloc] alpha_bloc;               // Intercepts of the blocks
  vector[nb_clon] alpha_clon;               // Intercepts of the clones
  real beta_0;                              // Global intercept
  real beta_X1;                             // Linear coefficent of GO or CTD
  real beta_X2;                             // Quadratic coefficent of GO or CTD
  real beta_H;                              // Linear coefficient of the initial height of the populations
  real<lower = 0>  sigma;                   // Residual variance of the model
}

transformed parameters {
  vector[N] mu;                             // Linear predictor
  real R_squared;                           // R^2 to evaluate the goodness of fit of the model
  
  mu = beta_0 + alpha_bloc[bloc] + alpha_clon[clon] + beta_X1 * X + beta_X2 * square(X) + beta_H * H;
  R_squared = 1 - variance(Y - mu) / variance(Y);
}

model {
  
  // Likelihood
  Y ~ normal(mu, sigma);
  
  // Priors
  sigma ~ exponential(1);
  alpha_bloc ~ std_normal();
  alpha_clon ~ std_normal();
  beta_0 ~ std_normal();
  beta_X1 ~ std_normal();
  beta_X2 ~ std_normal();
  beta_H ~ std_normal();
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
