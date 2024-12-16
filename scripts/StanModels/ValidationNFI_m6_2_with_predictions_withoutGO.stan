data {
  int N;
  vector[N] log_nb_years;                                                           // Offset to account for different census intervals
  vector[N] C;                                                                      // Proxy of the competition among trees (i.e. basal area)
  vector[N] DBH;                                                                    // Proxy of the average tree age (i.e. mean DBH)
  int nb_dead[N];                                                                   // Number of dead trees in the plot    
  int nb_tot[N];                                                                    // Total number of trees in the plot
  int<lower=0> nb_country;                                                          // Number of countries
  int<lower=0, upper=nb_country> country[N];                                        // Countries
}
parameters {
  vector[nb_country] alpha_country;
  real beta_C;
  real beta_DBH;
  real beta_C_DBH;
}

transformed parameters {
  vector[N] p;                               // probability of mortality during the inventory period of each plot
  vector[N] p_annual;                        // annual probability of mortality that can be explained
  //vector[N] res_mortality;                   // residual mortality
  
  p = inv_cloglog(alpha_country[country] + beta_C * C  + beta_DBH * DBH + beta_C_DBH * C .* DBH + log_nb_years);
  p_annual = inv_cloglog(alpha_country[country] + beta_C * C  + beta_DBH * DBH + beta_C_DBH * C .* DBH + log(1));
  //res_mortality = (nb_dead/nb_tot)/exp(log_nb_years) - p_annual;
  }

model {
  nb_dead ~ binomial(nb_tot,p);
  
  alpha_country ~ normal(0, 1);
  beta_C ~ normal(0, 1);
  beta_DBH ~ normal(0, 1);
  beta_C_DBH ~ normal(0, 1);
}
