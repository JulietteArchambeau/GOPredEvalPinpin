data {
  int N;
  vector[N] Y;                              // Response variable (individual height)
  vector[N] H;                              // Initial tree height of the populations
  int<lower=0> nb_bloc;                     // Number of blocks
  int<lower=0, upper=nb_bloc> bloc[N];      // Blocks
}

parameters {
  vector[nb_bloc] alpha_bloc;               // Intercepts of the blocks
  real beta_H;                              // Coefficient of the initial height of the populations
  real<lower = 0>  sigma;                   // Residual variance of the model
  
}

transformed parameters {
  vector[N] mu;    // Linear predictor
  real classic_R2;  // classic R2 to evaluate the goodness of fit of the model
  
  mu = alpha_bloc[bloc] + beta_H * H;
  classic_R2 = 1 - variance(Y - mu) / variance(Y); // classical R2
}

model {
  
  // Likelihood
  Y ~ normal(mu, sigma);
  
  // Priors
  sigma ~ exponential(1);
  alpha_bloc ~ std_normal();
  beta_H ~ std_normal();
}

generated quantities{
  // Bayesian R2
  real bayes_R2_res; // Residual-based based R2 (uses draws from the residual distribution)
  real bayes_R2_mod; // Model-based R2 (uses draws from the modeled residual variances)
  
  bayes_R2_res = variance(mu) / (variance(mu) + variance(Y-mu));
  bayes_R2_mod = variance(mu) / (variance(mu) + square(sigma));
}
