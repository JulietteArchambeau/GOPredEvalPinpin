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

transformed parameters {
  vector[N] p;                               // probability of mortality during the inventory period of each plot
  vector[N] p_annual;                        // annual probability of mortality
  vector[N] p_annual_go;                     // annual probability of mortality with all predictors except GO fixed
  
  p = inv_cloglog(alpha_country[country] + beta_GO * GO + beta_C * C  + beta_DBH * DBH + beta_C_DBH * C .* DBH + log_nb_years);
  p_annual = inv_cloglog(alpha_country[country] + beta_GO * GO + beta_C * C  + beta_DBH * DBH + beta_C_DBH * C .* DBH + 1);
  p_annual_go  = inv_cloglog(alpha_country[1] + beta_GO * GO + beta_C * mean(C)  + beta_DBH * mean(DBH) + beta_C_DBH * mean(C) .* mean(DBH) + 1);
}

model {
  nb_dead ~ binomial(nb_tot,p);
  
  alpha_country ~ normal(0, 1);
  beta_GO ~ normal(0, 1);
  beta_C ~ normal(0, 1);
  beta_DBH ~ normal(0, 1);
  beta_C_DBH ~ normal(0, 1);
}

generated quantities{
  vector[N] nb_dead_pred;
  vector[N] nb_dead_res;
  for (n in 1:N) {
    nb_dead_pred[n] = binomial_rng(nb_tot[n], p[n]);
    nb_dead_res[n] = nb_dead[n] - nb_dead_pred[n];
  }
}
