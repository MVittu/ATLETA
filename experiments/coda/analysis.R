options(stringsAsFactors = FALSE)

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
n_draws <- if (quick) 200L else 1000L
set.seed(20260719)

source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_2026-07-17_1419.csv"
out_root <- "experiments/coda/runs"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}

inv_spd <- function(x) chol2inv(chol(x))

close_composition <- function(x) x / rowSums(x)

replace_zeros <- function(x, pseudocount) {
  replaced <- x
  replaced[replaced == 0] <- pseudocount
  close_composition(replaced)
}

alr_transform <- function(composition, reference) {
  numerators <- setdiff(seq_len(ncol(composition)), reference)
  z <- sweep(log(composition[, numerators, drop = FALSE]), 1L,
             log(composition[, reference]), "-")
  colnames(z) <- paste0(colnames(composition)[numerators], "_vs_",
                        colnames(composition)[reference])
  z
}

alr_inverse <- function(z, reference, component_names) {
  numerators <- setdiff(seq_along(component_names), reference)
  eta <- matrix(0, nrow(z), length(component_names),
                dimnames = list(NULL, component_names))
  eta[, numerators] <- z
  eta <- sweep(eta, 1L, apply(eta, 1L, max), "-")
  close_composition(exp(eta))
}

posterior_parameters <- function(x, y) {
  p <- ncol(x)
  q <- ncol(y)
  b0 <- matrix(0, p, q)
  v0_inv <- diag(1 / 4, p)
  vn_inv <- v0_inv + crossprod(x)
  vn <- inv_spd(vn_inv)
  bn <- vn %*% crossprod(x, y)
  residual <- y - x %*% bn
  sn <- diag(q) + crossprod(residual) + crossprod(bn - b0, v0_inv %*% (bn - b0))
  list(b = bn, v = vn, s = sn, nu = q + 2L + nrow(x))
}

draw_posterior <- function(fit, draws) {
  p <- nrow(fit$b)
  q <- ncol(fit$b)
  b <- array(NA_real_, c(p, q, draws))
  sigma <- array(NA_real_, c(q, q, draws))
  lv <- t(chol(fit$v))
  for (s in seq_len(draws)) {
    sigma[, , s] <- solve(rWishart(1L, fit$nu, solve(fit$s))[, , 1L])
    b[, , s] <- fit$b + lv %*% matrix(rnorm(p * q), p, q) %*% chol(sigma[, , s])
  }
  list(b = b, sigma = sigma)
}

log_likelihood_rows <- function(y, mean, sigma) {
  q <- ncol(y)
  residual <- y - mean
  precision <- solve(sigma)
  log_det <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
  -0.5 * (q * log(2 * pi) + log_det + rowSums((residual %*% precision) * residual))
}

log_mean_exp <- function(x) {
  top <- max(x)
  top + log(mean(exp(x - top)))
}

log_multivariate_t <- function(y, location, scale, df) {
  q <- length(y)
  delta <- as.numeric(crossprod(y - location, solve(scale, y - location)))
  log_det <- as.numeric(determinant(scale, logarithm = TRUE)$modulus)
  lgamma((df + q) / 2) - lgamma(df / 2) -
    0.5 * (q * log(df * pi) + log_det) -
    0.5 * (df + q) * log1p(delta / df)
}

exact_loo <- function(x, y) {
  log_cpo <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    fit <- posterior_parameters(x[-i, , drop = FALSE], y[-i, , drop = FALSE])
    leverage <- 1 + as.numeric(x[i, , drop = FALSE] %*% fit$v %*% x[i, ])
    df <- fit$nu - ncol(y) + 1
    scale <- leverage * fit$s / df
    log_cpo[i] <- log_multivariate_t(y[i, ], as.numeric(x[i, ] %*% fit$b), scale, df)
  }
  list(elpd_loo = sum(log_cpo), lcpo = -mean(log_cpo), pointwise = log_cpo)
}

model_metrics <- function(x, y, draws) {
  fit <- posterior_parameters(x, y)
  posterior <- draw_posterior(fit, draws)
  log_lik <- matrix(NA_real_, draws, nrow(y))
  for (s in seq_len(draws)) {
    log_lik[s, ] <- log_likelihood_rows(y, x %*% posterior$b[, , s],
                                        posterior$sigma[, , s])
  }
  lppd <- sum(apply(log_lik, 2L, log_mean_exp))
  p_waic <- sum(apply(log_lik, 2L, var))
  d_bar <- -2 * mean(rowSums(log_lik))
  sigma_mean <- apply(posterior$sigma, c(1L, 2L), mean)
  d_hat <- -2 * sum(log_likelihood_rows(y, x %*% apply(posterior$b, c(1L, 2L), mean),
                                         sigma_mean))
  loo <- exact_loo(x, y)
  list(
    fit = fit,
    posterior = posterior,
    metrics = c(waic = -2 * (lppd - p_waic), p_waic = p_waic,
                dic = d_bar + (d_bar - d_hat), p_dic = d_bar - d_hat,
                elpd_loo = loo$elpd_loo, lcpo = loo$lcpo),
    loo_pointwise = loo$pointwise
  )
}

z_score <- function(x) {
  scale_value <- sd(x)
  if (!is.finite(scale_value) || scale_value == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x)) / scale_value)
}

raw <- read.csv(source_file, check.names = FALSE)
required <- c(
  "consenso", "eleggibile_calc", "survey_atleta_it_complete",
  "eta", "anni_atletica", "genere_gara", "statura_cm", "peso_kg",
  "infortunio_ultime_due_stagioni", paste0("disciplina___", 1:11)
)
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
part_columns <- c(
  paste0("prep_carbon_", domains), paste0("prep_no_carbon_", domains),
  paste0("gara_carbon_", domains), paste0("gara_no_carbon_", domains)
)
stopifnot(all(c(required, part_columns) %in% names(raw)))

eligible <- raw$consenso == 1 & raw$eleggibile_calc == 1 & raw$survey_atleta_it_complete == 2
raw <- raw[eligible, , drop = FALSE]
sum_columns <- function(prefix) rowSums(raw[paste0(prefix, domains)])
parts <- cbind(
  prep_carbon = sum_columns("prep_carbon_"),
  prep_noncarbon = sum_columns("prep_no_carbon_"),
  competition_carbon = sum_columns("gara_carbon_"),
  competition_noncarbon = sum_columns("gara_no_carbon_")
)
storage.mode(parts) <- "double"

height <- as.numeric(raw$statura_cm)
height[height > 0 & height < 3] <- height[height > 0 & height < 3] * 100
covariates <- data.frame(
  age_z = z_score(as.numeric(raw$eta)),
  sex_female = as.numeric(raw$genere_gara == 2),
  height_z = z_score(height),
  weight_z = z_score(as.numeric(raw$peso_kg)),
  experience_z = z_score(as.numeric(raw$anni_atletica)),
  prior_injury = as.numeric(raw$infortunio_ultime_due_stagioni),
  log_total_z = z_score(log(rowSums(parts)))
)
discipline_names <- paste0("discipline_", 1:11)
covariates[discipline_names] <- raw[paste0("disciplina___", 1:11)]
covariates$age_sq <- covariates$age_z^2 - mean(covariates$age_z^2)
covariates$experience_sq <- covariates$experience_z^2 - mean(covariates$experience_z^2)
covariates$age_sex <- covariates$age_z * covariates$sex_female
covariates$experience_prior <- covariates$experience_z * covariates$prior_injury

complete <- complete.cases(parts, covariates) & rowSums(parts) > 0
parts <- parts[complete, , drop = FALSE]
covariates <- covariates[complete, , drop = FALSE]
stopifnot(nrow(parts) > 0, all(parts >= 0), all(is.finite(as.matrix(covariates))))

full_main <- c("age_z", "sex_female", "height_z", "weight_z", "experience_z",
               "prior_injury", "log_total_z", discipline_names)
model_specs <- list(
  type_I_intercept = character(),
  type_II_demographics = c("age_z", "sex_female"),
  type_III_anthropometrics = c("age_z", "sex_female", "height_z", "weight_z"),
  type_IV_experience_history = c("age_z", "sex_female", "experience_z", "prior_injury"),
  type_V_discipline_profile = discipline_names,
  type_VI_full_main = full_main,
  type_VII_nonlinear = c(full_main, "age_sq", "experience_sq"),
  type_VIII_interactions = c(full_main, "age_sex", "experience_prior")
)

build_x <- function(data, terms) {
  x <- cbind(intercept = 1, as.matrix(data[, terms, drop = FALSE]))
  storage.mode(x) <- "double"
  x
}

fit_model_set <- function(composition, data, reference, run_dir, seed) {
  dir.create(file.path(run_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  y <- alr_transform(composition, reference)
  stopifnot(max(abs(alr_inverse(y, reference, colnames(composition)) - composition)) < 1e-10)
  fits <- vector("list", length(model_specs))
  comparison <- vector("list", length(model_specs))
  set.seed(seed)
  for (i in seq_along(model_specs)) {
    x <- build_x(data, model_specs[[i]])
    fits[[i]] <- model_metrics(x, y, n_draws)
    comparison[[i]] <- data.frame(
      model = names(model_specs)[i], predictors = length(model_specs[[i]]),
      n = nrow(x), reference = colnames(composition)[reference],
      t(fits[[i]]$metrics), check.names = FALSE
    )
  }
  comparison <- do.call(rbind, comparison)
  best_loo <- which.max(comparison$elpd_loo)
  comparison$elpd_diff <- comparison$elpd_loo - comparison$elpd_loo[best_loo]
  comparison$se_elpd_diff <- vapply(seq_along(fits), function(i) {
    difference <- fits[[i]]$loo_pointwise - fits[[best_loo]]$loo_pointwise
    sqrt(length(difference) * var(difference))
  }, numeric(1))
  comparison <- comparison[order(comparison$lcpo), ]
  rownames(comparison) <- NULL
  write_table(comparison, file.path(run_dir, "tables", "model_comparison.csv"))

  best_name <- comparison$model[1]
  best_index <- match(best_name, names(model_specs))
  best <- fits[[best_index]]
  x_best <- build_x(data, model_specs[[best_index]])
  fitted_composition <- alr_inverse(x_best %*% best$fit$b, reference, colnames(composition))
  posterior_average <- matrix(NA_real_, n_draws, ncol(composition))
  for (s in seq_len(n_draws)) {
    posterior_average[s, ] <- colMeans(alr_inverse(
      x_best %*% best$posterior$b[, , s], reference, colnames(composition)
    ))
  }
  component_summary <- data.frame(
    component = colnames(composition),
    observed_mean = colMeans(composition),
    posterior_mean = colMeans(posterior_average),
    lower_95 = apply(posterior_average, 2L, quantile, 0.025),
    upper_95 = apply(posterior_average, 2L, quantile, 0.975),
    rmse = sqrt(colMeans((fitted_composition - composition)^2)),
    mae = colMeans(abs(fitted_composition - composition))
  )
  write_table(component_summary, file.path(run_dir, "tables", "component_summary.csv"))
  coefficient_summary <- do.call(rbind, lapply(seq_len(dim(best$posterior$b)[1]), function(i) {
    do.call(rbind, lapply(seq_len(dim(best$posterior$b)[2]), function(j) {
      values <- best$posterior$b[i, j, ]
      data.frame(
        term = colnames(x_best)[i], coordinate = colnames(y)[j],
        posterior_mean = mean(values), posterior_sd = sd(values),
        lower_95 = quantile(values, 0.025), upper_95 = quantile(values, 0.975),
        probability_positive = mean(values > 0)
      )
    }))
  }))
  write_table(coefficient_summary, file.path(run_dir, "tables", "coefficient_summary.csv"))
  covariance <- as.data.frame(best$fit$s / (best$fit$nu - ncol(y) - 1))
  names(covariance) <- colnames(y)
  covariance$coordinate <- rownames(covariance) <- colnames(y)
  covariance <- covariance[, c("coordinate", colnames(y))]
  write_table(covariance, file.path(run_dir, "tables", "residual_covariance.csv"))
  saveRDS(list(
    model = best_name, terms = model_specs[[best_index]], reference = colnames(composition)[reference],
    component_names = colnames(composition), alr_names = colnames(y), posterior = best$fit,
    source_file = source_file
  ), file.path(run_dir, "models", "selected_model.rds"))

  png(file.path(run_dir, "plots", "model_selection.png"), width = 1200, height = 650, res = 140)
  old_par <- par(mfrow = c(1, 2), mar = c(4, 12, 3, 1))
  colors <- ifelse(seq_len(nrow(comparison)) == 1L, "#3478A8", "grey75")
  barplot(comparison$waic - min(comparison$waic), names.arg = comparison$model,
          horiz = TRUE, las = 1, cex.names = 0.65, col = colors,
          xlab = "Delta WAIC from minimum", main = "Joint posterior fit")
  barplot(comparison$lcpo - min(comparison$lcpo), names.arg = comparison$model,
          horiz = TRUE, las = 1, cex.names = 0.65, col = colors,
          xlab = "Delta mean -log CPO", main = "Athlete-level LOO")
  par(old_par)
  dev.off()

  png(file.path(run_dir, "plots", "component_fit.png"), width = 1000, height = 650, res = 140)
  values <- rbind(component_summary$observed_mean, component_summary$posterior_mean)
  centers <- barplot(values, beside = TRUE, names.arg = component_summary$component,
                     las = 2, cex.names = 0.75, col = c("grey65", "#3478A8"),
                     ylab = "Mean composition", main = best_name)
  arrows(centers[2, ], component_summary$lower_95, centers[2, ], component_summary$upper_95,
         angle = 90, code = 3, length = 0.04)
  legend("topright", c("Observed", "Posterior mean"), fill = c("grey65", "#3478A8"), bty = "n")
  dev.off()
  list(comparison = comparison, best = best_name)
}

audit_dir <- file.path(out_root, "00_data_audit", "tables")
positive <- rowSums(parts > 0) == ncol(parts)
positive_composition <- close_composition(parts[positive, , drop = FALSE])
log_variance <- apply(log(positive_composition), 2L, var)
reference <- which.min(log_variance)
write_table(data.frame(
  component = colnames(parts), zero_count = colSums(parts == 0),
  mean_raw_score = colMeans(parts), log_variance_positive = log_variance,
  selected_reference = seq_along(log_variance) == reference
), file.path(audit_dir, "composition_support.csv"))
write_table(data.frame(
  source_rows = nrow(read.csv(source_file)), eligible_rows = sum(eligible),
  complete_rows = nrow(parts), primary_all_positive_rows = sum(positive),
  excluded_from_primary_for_zero = sum(!positive), source_file = source_file
), file.path(audit_dir, "sample_flow.csv"))

primary <- fit_model_set(
  positive_composition, covariates[positive, , drop = FALSE], reference,
  file.path(out_root, "01_primary_positive"), 1001L
)

pseudocounts <- c(0.25, 0.5, 1)
sensitivity <- vector("list", length(pseudocounts))
for (i in seq_along(pseudocounts)) {
  label <- gsub("\\.", "p", format(pseudocounts[i], trim = TRUE))
  result <- fit_model_set(
    replace_zeros(parts, pseudocounts[i]), covariates, reference,
    file.path(out_root, "02_zero_sensitivity", paste0("pseudocount_", label)), 2000L + i
  )
  sensitivity[[i]] <- data.frame(
    pseudocount = pseudocounts[i], selected_model = result$best,
    lcpo = result$comparison$lcpo[1], waic = result$comparison$waic[1]
  )
}
write_table(do.call(rbind, sensitivity),
            file.path(out_root, "02_zero_sensitivity", "tables", "selection_summary.csv"))

reference_check <- list()
counter <- 1L
for (ref in seq_len(ncol(positive_composition))) {
  y <- alr_transform(positive_composition, ref)
  for (i in seq_along(model_specs)) {
    x <- build_x(covariates[positive, , drop = FALSE], model_specs[[i]])
    loo <- exact_loo(x, y)
    reference_check[[counter]] <- data.frame(
      reference = colnames(positive_composition)[ref], model = names(model_specs)[i],
      elpd_loo = loo$elpd_loo, lcpo = loo$lcpo
    )
    counter <- counter + 1L
  }
}
reference_check <- do.call(rbind, reference_check)
write_table(reference_check,
            file.path(out_root, "03_reference_check", "tables", "model_comparison.csv"))

write_table(data.frame(
  mode = if (quick) "quick" else "full", posterior_draws = n_draws,
  seed = 20260719, primary_reference = colnames(parts)[reference],
  primary_selected_model = primary$best,
  zero_replacement = "zeros only; positive scores unchanged; row then closed",
  covariance = "unrestricted 3x3 logistic-normal residual covariance",
  validation = "exact leave-one-athlete-out joint predictive density"
), file.path(out_root, "run_config.csv"))

cat("CoDa experiments complete. Primary model:", primary$best, "\n")
