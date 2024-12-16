data {
  int N;
  int nb_tot[N];    // Total number of trees in the population
  int nb_dead[N];   // Number of dead trees in the population
  vector[N] H;      // Mean tree height of the population
}

parameters {
  real beta_0;
  real beta_H;
}

model {
  nb_dead ~ binomial_logit(nb_tot,beta_0 + beta_H * H);
  
  beta_0 ~ normal(0,5);//std_normal();
  beta_H ~ normal(0,5);//std_normal();
}

