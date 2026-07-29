options(stringsAsFactors = FALSE)
source("experiments/multilevel/toolchain.R")

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
requested <- sub("^--exposure=", "", grep("^--exposure=", args, value = TRUE))
requested <- if (length(requested)) requested else c("share", "carbon_frequency")
stopifnot(all(requested %in% c("share", "carbon_frequency")))

seed <- 20260719L
set.seed(seed)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

root <- "experiments/multilevel"
out_root <- file.path(root, "runs", "06_exposure_comparison",
                      if (quick) "quick" else "final")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}

summarize_draws <- function(x) c(
  median = median(x),
  lower_95 = unname(quantile(x, 0.025)),
  upper_95 = unname(quantile(x, 0.975)),
  probability_positive = mean(x > 0)
)

z_score <- function(x) {
  s <- sd(x)
  stopifnot(is.finite(s), s > 0)
  as.numeric((x - mean(x)) / s)
}

raw_all <- read.csv(
  "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv",
  check.names = FALSE
)
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
required <- c(
  "consenso", "eleggibile_calc", "survey_atleta_it_complete", "infortunio_stagione",
  "inj_sing_periodo", "inj_multi1_periodo", "inj_multi2_periodo",
  "eta", "anni_atletica", "genere_gara", "statura_cm", "peso_kg",
  "infortunio_ultime_due_stagioni", paste0("disciplina___", 1:11),
  unlist(lapply(c("prep_carbon_", "prep_no_carbon_", "gara_carbon_", "gara_no_carbon_"),
                paste0, domains)),
  "prep_no_carbon_palestra", "gara_no_carbon_palestra"
)
stopifnot(all(required %in% names(raw_all)))

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
  hit <- rowSums((multiple == period_codes[1]) | (multiple == period_codes[2]),
                 na.rm = TRUE) > 0
  outcome[count == 2] <- as.integer(hit[count == 2])
  outcome
}

prep_injury <- derive_period_injury(c(1, 3))
competition_injury <- derive_period_injury(c(2, 3))
complete_timing <- raw$infortunio_stagione == 0 |
  (raw$infortunio_stagione == 1 & !is.na(raw$inj_sing_periodo)) |
  (raw$infortunio_stagione == 2 & !is.na(raw$inj_multi1_periodo) &
     !is.na(raw$inj_multi2_periodo))
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
X <- cbind(
  athlete_covariates[rep(seq_len(J), 2L), ],
  domain_totals,
  palestra_z = palestra
)
X <- as.matrix(X)
storage.mode(X) <- "double"

discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
active_codes <- which(discipline_counts > 0)
discipline_athlete <- as.matrix(raw[paste0("disciplina___", active_codes)])
colnames(discipline_athlete) <- paste0("discipline_", active_codes)
discipline <- discipline_athlete[rep(seq_len(J), 2L), , drop = FALSE]
storage.mode(discipline) <- "double"

stopifnot(
  J == 352L,
  length(y) == 2L * J,
  all(y %in% 0:1),
  all(is.finite(X)),
  all(is.finite(discipline)),
  all(share >= 0 & share <= 1)
)

standardize_users <- function(values) {
  center <- mean(values[any_carbon == 1])
  scale <- sd(values[any_carbon == 1])
  stopifnot(is.finite(scale), scale > 0)
  list(values = ifelse(any_carbon == 1, (values - center) / scale, 0),
       center = center, scale = scale)
}

exposures <- list(
  share = standardize_users(share),
  carbon_frequency = standardize_users(carbon)
)
exposure_labels <- c(
  share = "Carbon share",
  carbon_frequency = "Absolute carbon-frequency score"
)
write_table(do.call(rbind, lapply(names(exposures), function(name) data.frame(
  exposure = name,
  unit = if (name == "share") "proportion" else "summed ordinal frequency score",
  user_mean = exposures[[name]]$center,
  user_sd = exposures[[name]]$scale,
  user_min = min(if (name == "share") share[any_carbon == 1] else carbon[any_carbon == 1]),
  user_max = max(if (name == "share") share[any_carbon == 1] else carbon[any_carbon == 1])
))), file.path(out_root, "exposure_scales.csv"))

stan_data_for <- function(exposure) list(
  N = length(y),
  J = J,
  P = ncol(X),
  D = ncol(discipline),
  K = 1L,
  y = y,
  athlete = athlete,
  period = period,
  any_carbon = any_carbon,
  share_basis = matrix(exposure$values, ncol = 1L),
  X = X,
  discipline = discipline,
  include_exposure = 1L,
  prior_exposure_sd = 0.7,
  prior_covariate_sd = 0.5
)

sampling_config <- if (quick) {
  list(chains = 2L, iter = 600L, warmup = 300L)
} else {
  list(chains = 4L, iter = 2000L, warmup = 1000L)
}
model <- rstan::stan_model(file = file.path(root, "model.stan"))

x_reference <- apply(X, 2L, median)
modal_code <- active_codes[which.max(discipline_counts[active_codes])]
discipline_reference <- as.numeric(active_codes == modal_code)
write_table(rbind(
  data.frame(group = "covariate", term = names(x_reference), value = x_reference),
  data.frame(group = "discipline", term = paste0("discipline_", active_codes),
             value = discipline_reference)
), file.path(out_root, "reference_profile.csv"))

fit_one <- function(name, offset) {
  message("Fitting ", name)
  run_dir <- file.path(out_root, name)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  fit <- rstan::sampling(
    model,
    data = stan_data_for(exposures[[name]]),
    chains = sampling_config$chains,
    iter = sampling_config$iter,
    warmup = sampling_config$warmup,
    seed = seed + offset,
    refresh = if (quick) 0 else 200,
    control = list(adapt_delta = 0.98, max_treedepth = 12)
  )

  pars <- c("alpha", "beta_any", "beta_share", "beta_common",
            "beta_discipline", "tau_discipline", "sigma_athlete")
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)
  saveRDS(draws, file.path(run_dir, "posterior_draws.rds"))

  fit_summary <- rstan::summary(fit, pars = pars)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  diagnostics <- data.frame(
    exposure = name,
    max_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
    min_neff = min(fit_summary[, "n_eff"], na.rm = TRUE),
    divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
    max_treedepth_hits = sum(vapply(
      sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1)
    )),
    chains = sampling_config$chains,
    iterations = sampling_config$iter,
    warmup = sampling_config$warmup
  )
  write_table(diagnostics, file.path(run_dir, "sampling_diagnostics.csv"))

  coefficient_rows <- do.call(rbind, lapply(1:2, function(p) {
    values <- draws$beta_share[, p, 1]
    data.frame(
      exposure = name,
      period = c("preparation", "competition")[p],
      estimand = "log_odds_per_1_user_SD",
      t(summarize_draws(values)),
      odds_ratio_median = median(exp(values)),
      odds_ratio_lower_95 = unname(quantile(exp(values), 0.025)),
      odds_ratio_upper_95 = unname(quantile(exp(values), 0.975)),
      check.names = FALSE
    )
  }))
  write_table(coefficient_rows, file.path(run_dir, "standardized_coefficients.csv"))

  risk_rows <- lapply(1:2, function(p) {
    base <- draws$alpha[, p] +
      as.numeric(draws$beta_common %*% x_reference) +
      as.numeric(draws$beta_discipline %*% discipline_reference) +
      draws$beta_any[, p]
    beta <- draws$beta_share[, p, 1]
    low <- plogis(base - 0.5 * beta)
    high <- plogis(base + 0.5 * beta)
    difference <- high - low
    data.frame(
      exposure = name,
      period = c("preparation", "competition")[p],
      contrast = "one_SD_increase_among_users",
      raw_from = exposures[[name]]$center - 0.5 * exposures[[name]]$scale,
      raw_to = exposures[[name]]$center + 0.5 * exposures[[name]]$scale,
      risk_from_median = median(low),
      risk_to_median = median(high),
      t(summarize_draws(difference)),
      check.names = FALSE
    )
  })
  risk_draws <- lapply(1:2, function(p) {
    base <- draws$alpha[, p] +
      as.numeric(draws$beta_common %*% x_reference) +
      as.numeric(draws$beta_discipline %*% discipline_reference) +
      draws$beta_any[, p]
    beta <- draws$beta_share[, p, 1]
    plogis(base + 0.5 * beta) - plogis(base - 0.5 * beta)
  })
  names(risk_draws) <- c("preparation", "competition")
  write_table(do.call(rbind, risk_rows), file.path(run_dir, "standardized_risk_differences.csv"))

  rm(fit)
  gc()
  list(diagnostics = diagnostics, coefficients = coefficient_rows,
       risk = do.call(rbind, risk_rows), risk_draws = risk_draws)
}

fits <- setNames(vector("list", length(requested)), requested)
for (i in seq_along(requested)) {
  fits[[i]] <- fit_one(requested[i], i * 100L)
}

write_table(do.call(rbind, lapply(fits, `[[`, "diagnostics")),
            file.path(out_root, "sampling_diagnostics.csv"))
write_table(do.call(rbind, lapply(fits, `[[`, "coefficients")),
            file.path(out_root, "standardized_coefficients.csv"))
all_risks <- do.call(rbind, lapply(fits, `[[`, "risk"))
write_table(all_risks, file.path(out_root, "standardized_risk_differences.csv"))

if (all(c("share", "carbon_frequency") %in% names(fits))) {
  direct <- do.call(rbind, lapply(c("preparation", "competition"), function(p) {
    difference <- fits$share$risk_draws[[p]] - fits$carbon_frequency$risk_draws[[p]]
    data.frame(
      period = p,
      contrast = "share_minus_absolute_carbon_frequency",
      t(summarize_draws(difference)),
      check.names = FALSE
    )
  }))
  write_table(direct, file.path(out_root, "direct_standardized_comparison.csv"))

  png(file.path(out_root, "standardized_risk_differences.png"),
      width = 1600, height = 1000, res = 180)
  old_par <- par(mar = c(8, 5, 2, 1))
  plot(c(0.9, 2.1), range(c(all_risks$lower_95, all_risks$upper_95, 0)),
       type = "n", xaxt = "n", xlab = "", ylab = "Risk difference for +1 SD",
       bty = "l")
  abline(h = 0, lty = 2, col = "grey50")
  positions <- c(0.9, 1.1, 1.9, 2.1)
  colors <- rep(c("#2C7FB8", "#F28E2B"), 2)
  points(positions, all_risks$median, pch = 19, col = colors)
  arrows(positions, all_risks$lower_95, positions, all_risks$upper_95,
         angle = 90, code = 3, length = 0.05, col = colors)
  axis(1, at = positions,
       labels = paste(exposure_labels[all_risks$exposure], all_risks$period, sep = "\n"),
       las = 2, cex.axis = 0.8)
  par(old_par)
  dev.off()
}

writeLines(capture.output(sessionInfo()), file.path(out_root, "session_info.txt"))
message("Done: ", out_root)
