data {
  int N;
  int nb_tot[N];    // Total number of trees in the population
  int nb_dead[N];   // Number of dead trees in the population
  vector[N] H;      // Mean tree height of the population
  vector[N] X;      // Genomic offset or climatic transfer distance of the population
}

parameters {
  real beta_0;
  real beta_H;
  real beta_X;
}

model {
  nb_dead ~ binomial_logit(nb_tot,beta_0 + beta_H * H + beta_X * X);
  
  beta_0 ~ normal(0,5);//std_normal();
  beta_H ~ normal(0,5);//std_normal();
  beta_X ~ normal(0,5);//std_normal();
}

// generated quantities{
//   vector[N] log_lik;
//   // log likelihood for loo
//   for (n in 1:N) {
//     log_lik[n] = binomial_logit_lpmf( nb_dead[n] | nb_tot[n] , beta_0 + beta_H * H[n] + beta_X * X[n]);
//   }
// }
