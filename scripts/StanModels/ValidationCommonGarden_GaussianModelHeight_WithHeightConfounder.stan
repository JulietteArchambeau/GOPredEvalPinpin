data {
  int N;
  vector[N] Y;                              // Response variable (individual height)
  vector[N] X;                              // Genomic offset or climatic transfer distances
  vector[N] H;                              // Initial height of the populations
  int<lower=0> nb_bloc;                     // Number of blocks
  int<lower=0, upper=nb_bloc> bloc[N];      // Blocks
}

parameters {
  vector[nb_bloc] alpha_bloc;
  real beta_X1;
  real beta_X2;
  real beta_H;
  real<lower = 0>  sigma_r;
}

transformed parameters {
  vector[N] mu;    // linear predictor
  real R_squared;  // R^2 to evaluate the goodness of fit of the model
  
  mu = alpha_bloc[bloc] + beta_X1 * X + beta_X2 * square(X) + beta_H * H;
  R_squared = 1 - variance(H - mu) / variance(H);
}

model {
  
  Y ~ normal(mu, sigma_r); // Likelihood
  
  sigma_r ~ exponential(1);
  alpha_bloc ~ std_normal();
  beta_H ~ std_normal();
  beta_X1 ~ std_normal();
  beta_X2 ~ std_normal();
}

generated quantities{
  // log likelihood for loo
  vector[N] log_lik;
  vector[N] muhat;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(Y[n] |mu[n],sigma_r);
    muhat[n] = normal_rng(mu[n], sigma_r);
  }
}
