functions {
  // 1D Matern spectral density, angular-frequency convention (Riutort-Mayol et al. 2020).
  real spd_matern(real alpha, real rho, real nu, real w) {
    real numerator = 2 * sqrt(pi()) * tgamma(nu + 0.5) * pow(2 * nu, nu);
    real denominator = tgamma(nu) * pow(rho, 2 * nu);
    return square(alpha) * (numerator / denominator)
      * pow(2 * nu / square(rho) + square(w), -(nu + 0.5));
  }
}

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> P;
  int<lower=1> D;
  int<lower=1> M;
  int<lower=0, upper=1> y[N];
  int<lower=1, upper=J> athlete[N];
  int<lower=1, upper=2> period[N];
  vector<lower=0>[N] survey_weight;
  vector<lower=0, upper=1>[N] any_carbon;
  matrix[N, M] gp_basis;
  vector[M] sqrt_lambda;
  real<lower=0> matern_nu;
  matrix[N, P] X;
  matrix[N, D] discipline;
  real<lower=0> prior_exposure_sd;
  real<lower=0> prior_covariate_sd;
  real<lower=0> prior_gp_alpha_sd;
  real<lower=0> gp_rho_shape;
  real<lower=0> gp_rho_rate;
}

parameters {
  vector[2] alpha;
  vector[2] beta_any;
  matrix[2, M] beta_share_raw;
  vector<lower=0>[2] gp_alpha;
  vector<lower=0>[2] gp_rho;
  vector[P] beta_common;
  vector[D] beta_discipline_raw;
  real<lower=0> tau_discipline;
  vector[J] athlete_raw;
  real<lower=0> sigma_athlete;
}

transformed parameters {
  vector[D] beta_discipline = tau_discipline * beta_discipline_raw;
  vector[J] athlete_effect = sigma_athlete * athlete_raw;
  matrix[2, M] beta_share;
  for (p in 1:2) {
    for (m in 1:M) {
      beta_share[p, m] = beta_share_raw[p, m]
        * sqrt(spd_matern(gp_alpha[p], gp_rho[p], matern_nu, sqrt_lambda[m]));
    }
  }
}

model {
  alpha ~ normal(logit(0.4), 1.5);
  beta_any ~ normal(0, prior_exposure_sd);
  to_vector(beta_share_raw) ~ std_normal();
  gp_alpha ~ normal(0, prior_gp_alpha_sd);
  gp_rho ~ inv_gamma(gp_rho_shape, gp_rho_rate);
  beta_common ~ normal(0, prior_covariate_sd);
  beta_discipline_raw ~ std_normal();
  tau_discipline ~ normal(0, 0.5);
  athlete_raw ~ std_normal();
  sigma_athlete ~ normal(0, 1);

  for (n in 1:N) {
    real exposure = any_carbon[n] * beta_any[period[n]]
      + gp_basis[n] * beta_share[period[n]]';
    real eta = alpha[period[n]] + X[n] * beta_common
      + discipline[n] * beta_discipline + athlete_effect[athlete[n]] + exposure;
    target += survey_weight[n] * bernoulli_logit_lpmf(y[n] | eta);
  }
}
