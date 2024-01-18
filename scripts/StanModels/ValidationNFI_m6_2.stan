data {
  int N;
  vector[N] log_nb_years;                                                           // Offset to account for different census intervals
  vector[N] C;                                                                      // Proxy of the competition among trees (i.e. basal area)
  vector[N] DBH;                                                                    // Proxy of the average tree age (i.e. mean DBH)
  vector[N] GO;                                                                     // Genomic offset
  int nb_dead[N];                                                                   // Number of dead trees in the plot    
  int nb_tot[N];                                                                    // Total number of trees in the plot
  int<lower=0> nb_country;                                                          // Number of countries
  int<lower=0, upper=nb_country> country[N];                                        // Countries
}
parameters {
  vector[nb_country] alpha_country;
  real beta_GO;
  real beta_C;
  real beta_DBH;
  real beta_C_DBH;
}

model {
  nb_dead ~ binomial(nb_tot,inv_cloglog(alpha_country[country] + beta_GO * GO + beta_C * C  + beta_DBH * DBH + beta_C_DBH * C .* DBH + log_nb_years));
  
  alpha_country ~ normal(0, 1);
  beta_GO ~ normal(0, 1);
  beta_C ~ normal(0, 1);
  beta_DBH ~ normal(0, 1);
  beta_C_DBH ~ normal(0, 1);
}

