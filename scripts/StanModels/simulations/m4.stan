// cov_GPL2 macro extracted from ulam object with get_stancode
functions{
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
}
data {
  int N;
  vector[N] log_nb_years;                                                         // Offset to account for different census intervals
  vector[N] GO;                                                                   // Genomic offset
  vector[N] C;                                                                    // Proxy of the competition among trees (i.e. basal area)
  int nb_dead[N];                                                                 // Number of dead trees in the plot    
  int nb_tot[N];                                                                  // Total number of trees in the plot
  int<lower=0> nb_country;                                                        // Number of countries
  int<lower=0, upper=nb_country> country[N];                                      // Countries
  matrix[N, N] Dmat;                                                              // Matrix distance
}
parameters {
  vector[nb_country] alpha_country;
  real beta_GO;
  real beta_C;
  vector[N] mu;
  real<lower=0> etasq;
  real<lower=0> rhosq;
}

model {
  matrix[N, N] K;
  rhosq ~ exponential(0.5);
  etasq ~ exponential(2);
  K = cov_GPL2(Dmat, etasq, rhosq, 0.01);
  mu ~ multi_normal( rep_vector(0,N) , K );
  
  
  nb_dead ~ binomial(nb_tot,inv_cloglog(alpha_country[country] + beta_GO * GO +  beta_C * C  + log_nb_years + mu));
  
  alpha_country ~ normal(0, 1);
  beta_GO ~ normal(0, 1);
  beta_C ~ normal(0, 1);
}

