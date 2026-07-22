data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> P;
  int<lower=1> D;
  int<lower=1> K;
  int<lower=0, upper=1> y[N];
  int<lower=1, upper=J> athlete[N];
  int<lower=1, upper=2> period[N];
  vector<lower=0>[N] survey_weight;
  vector<lower=0, upper=1>[N] any_carbon;
  matrix[N, K] share_basis;
  matrix[N, P] X;
  matrix[N, D] discipline;
}

parameters {
  vector[2] alpha;
  vector[2] beta_any;
  matrix[2, K] beta_share;
  vector[P] beta_common;
  vector[D] beta_discipline_raw;
  real<lower=0> tau_discipline;
  vector[J] athlete_raw;
  real<lower=0> sigma_athlete;
}

transformed parameters {
  vector[D] beta_discipline = tau_discipline * beta_discipline_raw;
  vector[J] athlete_effect = sigma_athlete * athlete_raw;
}

model {
  alpha ~ normal(logit(0.4), 1.5);
  beta_any ~ normal(0, 0.7);
  to_vector(beta_share) ~ normal(0, 0.7);
  beta_common ~ normal(0, 0.5);
  beta_discipline_raw ~ std_normal();
  tau_discipline ~ normal(0, 0.5);
  athlete_raw ~ std_normal();
  sigma_athlete ~ normal(0, 1);

  for (n in 1:N) {
    real exposure = any_carbon[n] * beta_any[period[n]]
      + share_basis[n] * beta_share[period[n]]';
    real eta = alpha[period[n]] + X[n] * beta_common
      + discipline[n] * beta_discipline + athlete_effect[athlete[n]] + exposure;
    target += survey_weight[n] * bernoulli_logit_lpmf(y[n] | eta);
  }
}
