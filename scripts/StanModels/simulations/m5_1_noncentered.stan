
functions{
  matrix cov_GPL2(matrix x, real alpha, real rho, real delta) {
    real sq_alpha = square(alpha);
    real sq_rho = square(rho);
    int N = dims(x)[1];
    matrix[N, N] K;
    K = sq_alpha * exp(-0.5 *(square(x)/sq_rho)) + diag_matrix(rep_vector(delta, N));;
    return cholesky_decompose(K);
  }
}

data {
  int N;
  int nb_dead[N];                                                                 // Number of dead trees in the plot    
  int nb_tot[N];                                                                  // Total number of trees in the plot
  matrix[N, N] Dmat;                                                              // Distance matrix (in km)
}

transformed data {
  real delta = 1e-9; // Small offset to ensure the covariance matrix is positive definite
}


parameters {
  vector[N] eta;
  real beta_0;
  real<lower=0> alpha;
  real<lower=0> rho;
}

transformed parameters {
  matrix[N, N] L_K;
  vector[N] mu;
  L_K = cov_GPL2(Dmat, alpha, rho, delta);
  mu = L_K * eta;
  
}

model {
  rho ~ normal(500, 200);
  alpha ~ normal(0,1);
  beta_0 ~ normal(-4,2);
  eta ~ std_normal(); // Multiplier for non-centred GP parameterisation
  
  nb_dead ~ binomial(nb_tot,inv_cloglog(beta_0 + mu));

}

