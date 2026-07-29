options(stringsAsFactors = FALSE)
source("experiments/multilevel/toolchain.R")

# Primary analysis using a Hilbert-space approximate Gaussian process (HSGP) for the
# carbon-share exposure-response curve, replacing the fixed 3-df natural spline used in
# analysis.R. This script is self-contained (data loading is duplicated rather than
# sourced, matching the convention already used by fit_exposure_comparison.R) and does
# not modify analysis.R or its outputs.
#
# Design choices (see analyses-flow.md for the spline version):
#  - Kernel: Matern 5/2 (finite smoothness, closer to the natural cubic spline it
#    replaces than an infinitely smooth squared-exponential kernel would be).
#  - Lengthscale and marginal SD (amplitude) are period-specific: prep and competition
#    get separate gp_alpha/gp_rho, but share the same basis functions (same boundary L,
#    same M), mirroring how the spline basis was already pooled across periods while
#    beta_share coefficients were period-specific.
#  - The GP hyperparameter priors let the data decide curve flexibility instead of
#    fixing df=3 a priori; a post-fit diagnostic checks that M is large enough for the
#    HSGP approximation to be accurate given the posterior for the lengthscale.

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
seed <- 20260719L
set.seed(seed)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

root <- "experiments/multilevel"
run_root <- file.path(root, "runs", "07_hsgp")
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}

z_score <- function(x) {
  s <- sd(x)
  stopifnot(is.finite(s), s > 0)
  as.numeric((x - mean(x)) / s)
}

quantile_row <- function(x) {
  c(mean = mean(x), median = median(x), sd = sd(x),
    lower_95 = unname(quantile(x, 0.025)), upper_95 = unname(quantile(x, 0.975)),
    probability_positive = mean(x > 0), probability_gt_5pp = mean(x > 0.05),
    probability_rope_2pp = mean(abs(x) < 0.02))
}

estimand_row <- function(x, measure) {
  base <- c(mean = mean(x), median = median(x), sd = sd(x),
            lower_95 = unname(quantile(x, 0.025)), upper_95 = unname(quantile(x, 0.975)))
  probabilities <- switch(
    measure,
    prevalence_difference = c(
      probability_above_null = mean(x > 0), probability_above_threshold = mean(x > 0.05),
      probability_in_rope = mean(abs(x) < 0.02)
    ),
    prevalence_ratio = c(
      probability_above_null = mean(x > 1), probability_above_threshold = mean(x > 1.10),
      probability_in_rope = mean(x > 0.95 & x < 1.05)
    ),
    c(probability_above_null = NA, probability_above_threshold = NA, probability_in_rope = NA)
  )
  c(base, probabilities)
}

softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
log_sum_exp <- function(x) {
  top <- max(x)
  top + log(sum(exp(x - top)))
}

gauss_hermite_normal <- function(n = 15L) {
  jacobi <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  list(nodes = sqrt(2) * eig$values, weights = eig$vectors[1, ]^2)
}

# ---- HSGP machinery (Riutort-Mayol, Vehtari, Gelman, Andersen 2020) ------------------

matern_nu <- 2.5

spd_matern_r <- function(alpha, rho, nu, w) {
  numerator <- 2 * sqrt(pi) * gamma(nu + 0.5) * (2 * nu)^nu
  denominator <- gamma(nu) * rho^(2 * nu)
  alpha^2 * (numerator / denominator) * (2 * nu / rho^2 + w^2)^(-(nu + 0.5))
}

hsgp_basis_matrix <- function(x, L, sqrt_lambda) sin(outer(x + L, sqrt_lambda)) / sqrt(L)

# Solve inv_gamma(shape, rate) so that its 1%/99% quantiles match lower/upper targets.
# Standard recipe for weakly-informative GP lengthscale priors (Betancourt-style):
# avoids piling prior mass on implausibly short lengthscales the data can't resolve.
solve_inv_gamma <- function(lower, upper) {
  target_hi <- 1 / lower
  target_lo <- 1 / upper
  objective <- function(log_par) {
    shape <- exp(log_par[1]); rate <- exp(log_par[2])
    hi <- qgamma(0.99, shape = shape, rate = rate)
    lo <- qgamma(0.01, shape = shape, rate = rate)
    (log(hi) - log(target_hi))^2 + (log(lo) - log(target_lo))^2
  }
  fit <- optim(c(log(2), log(2)), objective, method = "BFGS")
  stopifnot(fit$convergence == 0, fit$value < 1e-4)
  list(shape = exp(fit$par[1]), rate = exp(fit$par[2]))
}

# ---- Data (duplicated from analysis.R; see that file for the annotated version) ------

raw_all <- read.csv(source_file, check.names = FALSE)
required <- c(
  "consenso", "eleggibile_calc", "survey_atleta_it_complete", "infortunio_stagione",
  "inj_sing_periodo", "inj_multi1_periodo", "inj_multi2_periodo",
  "eta", "anni_atletica", "genere_gara", "statura_cm", "peso_kg",
  "infortunio_ultime_due_stagioni", paste0("disciplina___", 1:11)
)
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
training_columns <- unlist(lapply(c("prep_carbon_", "prep_no_carbon_", "gara_carbon_", "gara_no_carbon_"),
                                  paste0, domains))
stopifnot(all(c(required, training_columns, "prep_no_carbon_palestra", "gara_no_carbon_palestra") %in%
                names(raw_all)))

eligible <- raw_all$consenso == 1 & raw_all$eleggibile_calc == 1 &
  raw_all$survey_atleta_it_complete == 2
raw <- raw_all[eligible, , drop = FALSE]
sum_named <- function(prefix) rowSums(raw[paste0(prefix, domains)])
prep_carbon <- sum_named("prep_carbon_")
prep_noncarbon <- sum_named("prep_no_carbon_")
competition_carbon <- sum_named("gara_carbon_")
competition_noncarbon <- sum_named("gara_no_carbon_")

derive_period_injury <- function(period_codes) {
  count <- raw$infortunio_stagione
  multiple <- cbind(raw$inj_multi1_periodo, raw$inj_multi2_periodo)
  outcome <- rep(NA_integer_, nrow(raw))
  outcome[count == 0] <- 0L
  outcome[count == 1] <- as.integer(raw$inj_sing_periodo[count == 1] %in% period_codes)
  multiple_hit <- rowSums((multiple == period_codes[1]) | (multiple == period_codes[2]), na.rm = TRUE) > 0
  outcome[count == 2] <- as.integer(multiple_hit[count == 2])
  outcome
}

prep_injury <- derive_period_injury(c(1, 3))
competition_injury <- derive_period_injury(c(2, 3))
complete_timing <- raw$infortunio_stagione == 0 |
  (raw$infortunio_stagione == 1 & !is.na(raw$inj_sing_periodo)) |
  (raw$infortunio_stagione == 2 & !is.na(raw$inj_multi1_periodo) & !is.na(raw$inj_multi2_periodo))
valid <- complete_timing & prep_carbon + prep_noncarbon > 0 &
  competition_carbon + competition_noncarbon > 0
raw <- raw[valid, , drop = FALSE]
prep_carbon <- prep_carbon[valid]
prep_noncarbon <- prep_noncarbon[valid]
competition_carbon <- competition_carbon[valid]
competition_noncarbon <- competition_noncarbon[valid]
prep_injury <- prep_injury[valid]
competition_injury <- competition_injury[valid]

J <- nrow(raw)
stopifnot(J == 352L)
period <- c(rep(1L, J), rep(2L, J))
athlete <- rep(seq_len(J), 2L)
y <- c(prep_injury, competition_injury)
carbon <- c(prep_carbon, competition_carbon)
noncarbon <- c(prep_noncarbon, competition_noncarbon)
share <- carbon / (carbon + noncarbon)
any_carbon <- as.numeric(carbon > 0)

height <- as.numeric(raw$statura_cm)
height[height > 0 & height < 3] <- height[height > 0 & height < 3] * 100
athlete_covariates <- data.frame(
  age_z = z_score(as.numeric(raw$eta)),
  sex_female = as.numeric(raw$genere_gara == 2),
  height_z = z_score(height),
  weight_z = z_score(as.numeric(raw$peso_kg)),
  experience_z = z_score(as.numeric(raw$anni_atletica)),
  prior_injury = as.numeric(raw$infortunio_ultime_due_stagioni)
)

domain_totals <- do.call(cbind, lapply(domains, function(domain) c(
  raw[[paste0("prep_carbon_", domain)]] + raw[[paste0("prep_no_carbon_", domain)]],
  raw[[paste0("gara_carbon_", domain)]] + raw[[paste0("gara_no_carbon_", domain)]]
)))
colnames(domain_totals) <- paste0(domains, "_total_z")
domain_totals <- apply(domain_totals, 2L, z_score)
palestra <- z_score(c(raw$prep_no_carbon_palestra, raw$gara_no_carbon_palestra))
X <- cbind(athlete_covariates[rep(seq_len(J), 2L), ], domain_totals, palestra_z = palestra)
X <- as.matrix(X)
storage.mode(X) <- "double"

discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
active_discipline_codes <- which(discipline_counts > 0)
discipline <- as.matrix(raw[paste0("disciplina___", active_discipline_codes)])
colnames(discipline) <- paste0("discipline_", active_discipline_codes)
discipline <- discipline[rep(seq_len(J), 2L), , drop = FALSE]
storage.mode(discipline) <- "double"

stopifnot(J > 0L, length(y) == 2L * J, all(y %in% 0:1), all(is.finite(X)),
          all(is.finite(discipline)), all(share >= 0 & share <= 1))

# ---- Exposure bases --------------------------------------------------------------

share_center <- median(share[any_carbon == 1])

# HSGP basis, pooled across periods (same convention as the ns() spline it replaces):
# boundary factor c and number of basis functions M are fixed generous defaults; the
# post-fit spd_ratio diagnostic below verifies M is large enough for the posterior
# lengthscale, rather than trusting the defaults blindly.
hsgp_c <- 1.5
hsgp_M <- 40L
x_users <- share[any_carbon == 1] - share_center
hsgp_L <- hsgp_c * max(abs(x_users))
sqrt_lambda <- seq_len(hsgp_M) * pi / (2 * hsgp_L)

range_x <- diff(range(x_users))
gp_rho_prior <- solve_inv_gamma(lower = range_x / 10, upper = range_x * 2)

basis_gp <- function(values, used = rep(1, length(values))) {
  out <- matrix(0, length(values), hsgp_M)
  index <- used == 1
  if (any(index)) out[index, ] <- hsgp_basis_matrix(values[index] - share_center, hsgp_L, sqrt_lambda)
  out
}
basis_linear <- function(values, used = rep(1, length(values))) {
  out <- matrix(0, length(values), 1L)
  out[used == 1, 1] <- values[used == 1] - share_center
  out
}

write_table(data.frame(
  boundary_factor_c = hsgp_c, boundary_L = hsgp_L, num_basis_functions = hsgp_M,
  matern_nu = matern_nu, share_center = share_center,
  gp_rho_prior_shape = gp_rho_prior$shape, gp_rho_prior_rate = gp_rho_prior$rate,
  gp_rho_prior_1pct = 1 / qgamma(0.99, gp_rho_prior$shape, gp_rho_prior$rate),
  gp_rho_prior_99pct = 1 / qgamma(0.01, gp_rho_prior$shape, gp_rho_prior$rate)
), file.path(run_root, "hsgp_configuration.csv"))

# ---- Model compilation ------------------------------------------------------------

model_hsgp <- rstan::stan_model(file = file.path(root, "model_hsgp.stan"))
model_plain <- rstan::stan_model(file = file.path(root, "model.stan"))

sampling_config <- if (quick) {
  list(chains = 2L, iter = 600L, warmup = 300L)
} else {
  list(chains = 4L, iter = 2000L, warmup = 1000L)
}

configs <- list(
  primary_gp_association = list(
    folder = "01_primary/association", stan_model = model_hsgp, model_type = "hsgp",
    basis = "gp", include = 1L, exposure_sd = 0.7, covariate_sd = 0.5, gp_alpha_sd = 0.7
  ),
  primary_gp_baseline = list(
    folder = "01_primary/baseline", stan_model = model_plain, model_type = "plain",
    basis = "linear", include = 0L, exposure_sd = 0.7, covariate_sd = 0.5
  ),
  sensitivity_linear_share = list(
    folder = "02_sensitivity/linear_share", stan_model = model_plain, model_type = "plain",
    basis = "linear", include = 1L, exposure_sd = 0.7, covariate_sd = 0.5
  ),
  sensitivity_gp_wide_prior = list(
    folder = "02_sensitivity/wide_prior", stan_model = model_hsgp, model_type = "hsgp",
    basis = "gp", include = 1L, exposure_sd = 0.7, covariate_sd = 1.0, gp_alpha_sd = 1.2
  )
)

gh <- gauss_hermite_normal(15L)

stan_data_for <- function(config) {
  basis <- if (config$basis == "gp") basis_gp(share, any_carbon) else basis_linear(share, any_carbon)
  data_list <- list(
    N = length(y), J = J, P = ncol(X), D = ncol(discipline),
    y = y, athlete = athlete, period = period, any_carbon = any_carbon,
    X = X, discipline = discipline,
    include_exposure = config$include,
    prior_exposure_sd = config$exposure_sd,
    prior_covariate_sd = config$covariate_sd
  )
  if (config$model_type == "hsgp") {
    data_list$M <- hsgp_M
    data_list$gp_basis <- basis
    data_list$sqrt_lambda <- sqrt_lambda
    data_list$matern_nu <- matern_nu
    data_list$prior_gp_alpha_sd <- config$gp_alpha_sd
    data_list$gp_rho_shape <- gp_rho_prior$shape
    data_list$gp_rho_rate <- gp_rho_prior$rate
  } else {
    data_list$K <- ncol(basis)
    data_list$share_basis <- basis
  }
  list(data = data_list, basis = basis)
}

fixed_eta <- function(draws, config, stan_data, s) {
  eta <- draws$alpha[s, period] + as.numeric(X %*% draws$beta_common[s, ]) +
    as.numeric(discipline %*% draws$beta_discipline[s, ])
  if (config$include == 1L) {
    share_effect <- rowSums(stan_data$basis * draws$beta_share[s, period, ])
    eta <- eta + any_carbon * draws$beta_any[s, period] + share_effect
  }
  eta
}

grouped_log_likelihood <- function(draws, config, stan_data) {
  samples <- nrow(draws$alpha)
  output <- matrix(NA_real_, samples, J)
  log_weights <- log(gh$weights)
  for (s in seq_len(samples)) {
    eta <- fixed_eta(draws, config, stan_data, s)
    shift <- draws$sigma_athlete[s] * gh$nodes
    prep_eta <- outer(eta[seq_len(J)], shift, "+")
    competition_eta <- outer(eta[J + seq_len(J)], shift, "+")
    prep_log <- -softplus(prep_eta)
    prep_positive <- y[seq_len(J)] == 1L
    prep_log[prep_positive, ] <- -softplus(-prep_eta[prep_positive, , drop = FALSE])
    competition_log <- -softplus(competition_eta)
    competition_positive <- y[J + seq_len(J)] == 1L
    competition_log[competition_positive, ] <-
      -softplus(-competition_eta[competition_positive, , drop = FALSE])
    joint <- sweep(prep_log + competition_log, 2L, log_weights, "+")
    output[s, ] <- apply(joint, 1L, log_sum_exp)
  }
  output
}

marginal_prevalence <- function(eta, sigma) {
  probabilities <- plogis(outer(eta, sigma * gh$nodes, "+"))
  as.numeric(probabilities %*% gh$weights)
}

summarize_estimands <- function(draws, config, basis_function) {
  samples <- nrow(draws$alpha)
  output <- list()
  counter <- 1L
  for (p in 1:2) {
    index <- which(period == p)
    user_share <- share[index][any_carbon[index] == 1]
    values <- c(low = unname(quantile(user_share, 0.25)), median = median(user_share),
                high = unname(quantile(user_share, 0.75)))
    basis_values <- basis_function(values)
    rownames(basis_values) <- names(values)
    scenario <- matrix(NA_real_, samples, 4L,
                       dimnames = list(NULL, c("nonuser", "median_user", "low_user", "high_user")))
    for (s in seq_len(samples)) {
      base <- draws$alpha[s, p] + as.numeric(X[index, , drop = FALSE] %*% draws$beta_common[s, ]) +
        as.numeric(discipline[index, , drop = FALSE] %*% draws$beta_discipline[s, ])
      scenario[s, "nonuser"] <- mean(marginal_prevalence(base, draws$sigma_athlete[s]))
      for (name in c("median", "low", "high")) {
        exposure <- draws$beta_any[s, p] +
          sum(basis_values[name, ] * draws$beta_share[s, p, ])
        scenario[s, paste0(name, "_user")] <-
          mean(marginal_prevalence(base + exposure, draws$sigma_athlete[s]))
      }
    }
    period_name <- c("preparation", "competition")[p]
    for (name in colnames(scenario)) {
      output[[counter]] <- data.frame(
        period = period_name, estimand = paste0("prevalence_", name), measure = "prevalence",
        share_from = NA, share_to = if (name == "nonuser") 0 else values[sub("_user", "", name)],
        t(estimand_row(scenario[, name], "prevalence")), check.names = FALSE
      )
      counter <- counter + 1L
    }
    contrasts <- list(
      adoption_difference = scenario[, "median_user"] - scenario[, "nonuser"],
      adoption_ratio = scenario[, "median_user"] / scenario[, "nonuser"],
      dose_difference = scenario[, "high_user"] - scenario[, "low_user"],
      dose_ratio = scenario[, "high_user"] / scenario[, "low_user"]
    )
    for (name in names(contrasts)) {
      dose <- grepl("^dose", name)
      measure <- if (grepl("difference$", name)) "prevalence_difference" else "prevalence_ratio"
      output[[counter]] <- data.frame(
        period = period_name, estimand = name,
        measure = measure,
        share_from = if (dose) values["low"] else 0,
        share_to = if (dose) values["high"] else values["median"],
        t(estimand_row(contrasts[[name]], measure)), check.names = FALSE
      )
      counter <- counter + 1L
    }
  }
  do.call(rbind, output)
}

parameter_summary <- function(draws, config) {
  rows <- list()
  add <- function(group, term, period_name, values) {
    rows[[length(rows) + 1L]] <<- data.frame(
      group = group, term = term, period = period_name,
      t(quantile_row(values)), check.names = FALSE
    )
  }
  for (p in 1:2) {
    name <- c("preparation", "competition")[p]
    add("intercept", "intercept", name, draws$alpha[, p])
    if (config$include == 1L) add("exposure", "any_carbon_at_center_share", name, draws$beta_any[, p])
    if (config$include == 1L && config$model_type == "hsgp") {
      add("exposure_gp", "gp_marginal_sd", name, draws$gp_alpha[, p])
      add("exposure_gp", "gp_lengthscale", name, draws$gp_rho[, p])
    }
    if (config$include == 1L && config$model_type == "plain")
      for (k in seq_len(dim(draws$beta_share)[3]))
        add("exposure", paste0("share_basis_", k), name, draws$beta_share[, p, k])
  }
  for (i in seq_len(ncol(X))) add("covariate", colnames(X)[i], "common", draws$beta_common[, i])
  for (i in seq_len(ncol(discipline)))
    add("discipline", colnames(discipline)[i], "common", draws$beta_discipline[, i])
  add("scale", "sigma_athlete", "common", draws$sigma_athlete)
  add("scale", "tau_discipline", "common", draws$tau_discipline)
  do.call(rbind, rows)
}

hsgp_approximation_diagnostic <- function(draws, name) {
  do.call(rbind, lapply(1:2, function(p) {
    ratio <- spd_matern_r(draws$gp_alpha[, p], draws$gp_rho[, p], matern_nu, sqrt_lambda[hsgp_M]) /
      spd_matern_r(draws$gp_alpha[, p], draws$gp_rho[, p], matern_nu, sqrt_lambda[1])
    data.frame(
      model = name, period = c("preparation", "competition")[p],
      max_spd_ratio_M_over_1 = max(ratio), q99_spd_ratio_M_over_1 = unname(quantile(ratio, 0.99)),
      adequate_M = max(ratio) < 1e-2
    )
  }))
}

fit_one <- function(name, config, offset) {
  message("Fitting ", name)
  stan_data <- stan_data_for(config)
  run_dir <- file.path(run_root, config$folder)
  dir.create(file.path(run_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  # Written incrementally by Stan as sampling proceeds (one row per saved
  # iteration, warmup included), so it can be tailed live from another R
  # session with watch_progress.R while this chain is still running -- rstan's
  # own `refresh` progress text does not stream back from parallel PSOCK
  # workers on Windows, so this is the reliable way to see progress firsthand.
  progress_stub <- file.path(run_dir, "models", "chain.csv")
  fit <- rstan::sampling(
    config$stan_model, data = stan_data$data, chains = sampling_config$chains,
    iter = sampling_config$iter, warmup = sampling_config$warmup,
    seed = seed + offset, refresh = if (quick) 0 else 200,
    sample_file = progress_stub,
    control = list(adapt_delta = 0.98, max_treedepth = 12)
  )
  monitored <- c("alpha", "beta_common", "beta_discipline", "tau_discipline", "sigma_athlete")
  if (config$include == 1L) monitored <- c(monitored, "beta_any", "beta_share")
  if (config$model_type == "hsgp") monitored <- c(monitored, "gp_alpha", "gp_rho")
  summary_matrix <- rstan::summary(fit, pars = monitored)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  ebfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2) / var(x[, "energy__"]), numeric(1))
  diagnostics <- data.frame(
    model = name, max_rhat = max(summary_matrix[, "Rhat"], na.rm = TRUE),
    min_bulk_proxy_neff = min(summary_matrix[, "n_eff"], na.rm = TRUE),
    divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
    max_treedepth_hits = sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))),
    min_ebfmi = min(ebfmi), chains = sampling_config$chains,
    iterations = sampling_config$iter, warmup = sampling_config$warmup
  )
  write_table(diagnostics, file.path(run_dir, "tables", "sampling_diagnostics.csv"))

  extract_pars <- c(monitored, "injury_rep")
  if (config$include == 0L) extract_pars <- setdiff(extract_pars, c("beta_any", "beta_share"))
  draws <- rstan::extract(fit, pars = extract_pars, permuted = TRUE)
  compact_pars <- setdiff(extract_pars, "injury_rep")
  compact_draws <- draws[compact_pars]
  saveRDS(compact_draws, file.path(run_dir, "models", "posterior_draws.rds"))
  write_table(parameter_summary(draws, config), file.path(run_dir, "tables", "parameter_summary.csv"))

  if (config$include == 1L && config$model_type == "hsgp") {
    write_table(hsgp_approximation_diagnostic(draws, name),
                file.path(run_dir, "tables", "hsgp_approximation_diagnostics.csv"))
  }

  ppc <- do.call(rbind, lapply(1:2, function(p) data.frame(
    period = c("preparation", "competition")[p], observed = sum(y[period == p]),
    posterior_mean = mean(draws$injury_rep[, p]),
    lower_95 = quantile(draws$injury_rep[, p], 0.025),
    upper_95 = quantile(draws$injury_rep[, p], 0.975),
    bayesian_p_upper = mean(draws$injury_rep[, p] >= sum(y[period == p]))
  )))
  write_table(ppc, file.path(run_dir, "tables", "posterior_predictive_check.csv"))

  grouped_ll <- grouped_log_likelihood(draws, config, stan_data)
  loo_result <- loo::loo(grouped_ll, cores = min(4L, parallel::detectCores()))
  pareto <- loo::pareto_k_values(loo_result)
  write_table(data.frame(
    model = name, elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    looic = loo_result$estimates["looic", "Estimate"],
    pareto_k_le_0p5 = sum(pareto <= 0.5), pareto_k_0p5_0p7 = sum(pareto > 0.5 & pareto <= 0.7),
    pareto_k_0p7_1 = sum(pareto > 0.7 & pareto <= 1), pareto_k_gt_1 = sum(pareto > 1)
  ), file.path(run_dir, "tables", "grouped_psis_loo.csv"))

  estimands <- NULL
  if (config$include == 1L) {
    basis_function <- if (config$basis == "gp") basis_gp else basis_linear
    estimands <- summarize_estimands(draws, config, basis_function)
    write_table(estimands, file.path(run_dir, "tables", "association_estimands.csv"))
  }
  rm(fit)
  gc()
  unlink(Sys.glob(paste0(sub("\\.csv$", "", progress_stub), "*.csv")))
  list(draws = draws, loo = loo_result, grouped_ll = grouped_ll,
       diagnostics = diagnostics, estimands = estimands, ppc = ppc)
}

fits <- vector("list", length(configs))
names(fits) <- names(configs)
for (i in seq_along(configs)) fits[[i]] <- fit_one(names(configs)[i], configs[[i]], i * 100L)

primary_comparison <- loo::loo_compare(
  lapply(fits[c("primary_gp_association", "primary_gp_baseline")], `[[`, "loo")
)
write_table(data.frame(model = rownames(primary_comparison), primary_comparison,
                       row.names = NULL, check.names = FALSE),
            file.path(run_root, "01_primary", "tables", "model_comparison.csv"))

all_comparison <- loo::loo_compare(lapply(fits, `[[`, "loo"))
write_table(data.frame(model = rownames(all_comparison), all_comparison,
                       row.names = NULL, check.names = FALSE),
            file.path(run_root, "02_sensitivity", "tables", "all_model_comparison.csv"))

pointwise_difference <- fits$primary_gp_association$loo$pointwise[, "elpd_loo"] -
  fits$primary_gp_baseline$loo$pointwise[, "elpd_loo"]
write_table(data.frame(
  contrast = "association_minus_baseline", elpd_difference = sum(pointwise_difference),
  se_difference = sqrt(J * var(pointwise_difference))
), file.path(run_root, "01_primary", "tables", "primary_predictive_difference.csv"))

sensitivity_estimands <- do.call(rbind, lapply(
  c("primary_gp_association", "sensitivity_linear_share", "sensitivity_gp_wide_prior"),
  function(name) cbind(model = name, fits[[name]]$estimands)
))
write_table(sensitivity_estimands,
            file.path(run_root, "02_sensitivity", "tables", "all_association_estimands.csv"))

primary_estimands <- read.csv(file.path(
  run_root, "01_primary", "association", "tables", "association_estimands.csv"
))
forest <- primary_estimands[primary_estimands$measure == "prevalence_difference", ]
png(file.path(run_root, "01_primary", "association", "plots", "association_contrasts.png"),
    width = 1000, height = 620, res = 140)
labels <- paste(forest$period, sub("_difference", "", forest$estimand), sep = ": ")
positions <- rev(seq_len(nrow(forest)))
plot(forest$mean, positions, xlim = range(c(forest$lower_95, forest$upper_95, 0)),
     ylim = c(0.5, nrow(forest) + 0.5), yaxt = "n", ylab = "", xlab = "Adjusted prevalence difference",
     pch = 19, col = "#3478A8")
segments(forest$lower_95, positions, forest$upper_95, positions, col = "#3478A8", lwd = 2)
abline(v = 0, lty = 2, col = "grey45")
axis(2, at = positions, labels = labels, las = 1)
dev.off()

ppc <- read.csv(file.path(
  run_root, "01_primary", "association", "tables", "posterior_predictive_check.csv"
))
png(file.path(run_root, "01_primary", "association", "plots", "posterior_predictive_injuries.png"),
    width = 850, height = 600, res = 140)
centers <- barplot(rbind(ppc$observed, ppc$posterior_mean), beside = TRUE,
                   names.arg = ppc$period, col = c("grey65", "#3478A8"),
                   ylab = "Athletes with injury", ylim = c(0, max(ppc$upper_95) * 1.12))
arrows(centers[2, ], ppc$lower_95, centers[2, ], ppc$upper_95,
       angle = 90, code = 3, length = 0.04)
legend("topleft", c("Observed", "Posterior predictive"), fill = c("grey65", "#3478A8"), bty = "n")
dev.off()

write_table(data.frame(
  mode = if (quick) "quick" else "full", seed = seed,
  chains = sampling_config$chains, iterations = sampling_config$iter,
  warmup = sampling_config$warmup, athletes = J, period_records = length(y),
  primary_model = "Bayesian repeated-period logistic model with athlete random intercept",
  primary_exposure = "any carbon use plus Hilbert-space approximate GP (Matern 5/2) of carbon share among users",
  primary_estimands = "marginal adjusted prevalence differences and ratios",
  validation = "athlete-grouped marginal PSIS-LOO using 15-point Gauss-Hermite quadrature",
  hsgp_boundary_L = hsgp_L, hsgp_num_basis = hsgp_M, hsgp_boundary_factor_c = hsgp_c
), file.path(run_root, "run_config.csv"))

stopifnot(all(is.finite(primary_estimands$mean)), nrow(primary_estimands) == 16L)
if (!quick) stopifnot(
  all(vapply(fits, function(x) x$diagnostics$divergences == 0, logical(1))),
  all(vapply(fits, function(x) x$diagnostics$max_rhat < 1.01, logical(1))),
  all(hsgp_approximation_diagnostic(fits$primary_gp_association$draws, "primary_gp_association")$adequate_M),
  all(hsgp_approximation_diagnostic(fits$sensitivity_gp_wide_prior$draws, "sensitivity_gp_wide_prior")$adequate_M)
)
cat("HSGP multilevel experiments complete.\n")
