// This model does not work
functions{
  matrix cov_GPL2(matrix x, real alpha, real rho, real delta) {
    real sq_alpha = square(alpha);
    real sq_rho = square(rho);
    int N = dims(x)[1];
    matrix[N, N] K;
    K = sq_alpha * exp(-0.5 *(square(x)/sq_rho)) + diag_matrix(rep_vector(delta, N));;
    return K;
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
  vector[N] mu;
  real beta_0;
  real<lower=0> alpha;
  real<lower=0> rho;
}

model {
  matrix[N, N] K;
  rho ~ normal(500, 200);
  alpha ~ normal(0,0.8);
  K = cov_GPL2(Dmat, alpha, rho, delta);
  mu ~ multi_normal( rep_vector(0,N) , K );
  
  
  nb_dead ~ binomial(nb_tot,inv_cloglog(beta_0 + mu));
  

}

