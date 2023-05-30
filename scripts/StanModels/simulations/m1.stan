data {
  int N;
  vector[N] log_nb_years;                                                           // Offset to account for different census intervals
  vector[N] GO;                                                                     // Genomic offset
  int nb_dead[N];                                                                   // Number of dead trees in the plot    
  int nb_tot[N];                                                                    // Total number of trees in the plot
}
parameters {
  real beta_0;
  real beta_GO;
}

model {
  nb_dead ~ binomial(nb_tot,inv_cloglog(beta_0 + beta_GO * GO + log_nb_years));
  
  beta_0 ~ normal(0, 1);
  beta_GO ~ normal(0, 1);
}

