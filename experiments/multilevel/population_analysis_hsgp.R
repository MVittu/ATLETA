options(stringsAsFactors = FALSE)
source("experiments/multilevel/toolchain.R")

# HSGP counterpart of population_analysis.R: poststratification and sex/age
# conditional estimands for the primary HSGP curve (runs/07_hsgp) rather than
# the natural-spline curve. Self-contained (data loading duplicated, matching
# the convention of fit_exposure_comparison.R / analysis_hsgp.R) and does not
# touch population_analysis.R or its runs/03-04 outputs, which remain the
# spline-era reference.

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
seed <- 20260719L
set.seed(seed)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

root <- "experiments/multilevel"
hsgp_root <- file.path(root, "runs", "07_hsgp")
weighted_root <- file.path(root, "runs/08_hsgp_population_weighting")
subgroup_root <- file.path(root, "runs/09_hsgp_subgroups")
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
population_file <- "tables/recall/population_history_by_age_sex.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE)
}
z_score <- function(x) as.numeric((x - mean(x)) / sd(x))
weighted_quantile <- function(x, w, probs) {
  order_x <- order(x)
  x <- x[order_x]
  cumulative <- cumsum(w[order_x]) / sum(w)
  vapply(probs, function(prob) x[which(cumulative >= prob)[1]], numeric(1))
}
summary_row <- function(x, measure) {
  base <- c(mean = mean(x), median = median(x), sd = sd(x),
            lower_95 = unname(quantile(x, 0.025)), upper_95 = unname(quantile(x, 0.975)))
  probability <- switch(
    measure,
    prevalence_difference = c(
      probability_above_null = mean(x > 0), probability_above_threshold = mean(x > 0.05),
      probability_in_rope = mean(abs(x) < 0.02)
    ),
    prevalence_ratio = c(
      probability_above_null = mean(x > 1), probability_above_threshold = mean(x > 1.10),
      probability_in_rope = mean(x > 0.95 & x < 1.05)
    ),
    c(probability_above_null = NA, probability_above_threshold = NA,
      probability_in_rope = NA)
  )
  c(base, probability)
}

raw_all <- read.csv(source_file, check.names = FALSE)
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
  multiple_hit <- rowSums((multiple == period_codes[1]) | (multiple == period_codes[2]),
                          na.rm = TRUE) > 0
  outcome[count == 2] <- as.integer(multiple_hit[count == 2])
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
sex <- factor(ifelse(raw$genere_gara == 1, "Men", "Women"), c("Men", "Women"))
age_group <- cut(as.numeric(raw$eta), c(-Inf, 19, 22, Inf),
                 labels = c("Under 20", "Under 23", "Senior"))

# Post-stratify to the championship frame used by the existing recall plan.
population_history <- read.csv(population_file)
population_history$age_group <- factor(
  population_history$age_group, c("under20", "under23", "assoluta")
)
target <- do.call(rbind, lapply(split(
  population_history, list(population_history$age_group, population_history$sex), drop = TRUE
), function(cell) {
  current <- cell$n[cell$year == 2026]
  estimate <- if (length(current)) current else round(mean(cell$n[cell$year < 2026]))
  data.frame(age_group = as.character(cell$age_group[1]), sex = cell$sex[1], target_n = estimate)
}))
target$age_group <- c(under20 = "Under 20", under23 = "Under 23", assoluta = "Senior")[target$age_group]
target$sex <- c(uomini = "Men", donne = "Women")[target$sex]
sample_cells <- as.data.frame(table(age_group = age_group, sex = sex), responseName = "sample_n")
weight_cells <- merge(target, sample_cells, by = c("age_group", "sex"), sort = FALSE)
stopifnot(nrow(weight_cells) == 6L, all(weight_cells$sample_n > 0), sum(weight_cells$target_n) == 1328L)
weight_cells$raw_weight <- weight_cells$target_n / weight_cells$sample_n
weight_cells$model_weight <- weight_cells$raw_weight / weighted.mean(
  weight_cells$raw_weight, weight_cells$sample_n
)
cell_key <- paste(age_group, sex)
weight_key <- paste(weight_cells$age_group, weight_cells$sex)
raw_weight <- weight_cells$raw_weight[match(cell_key, weight_key)]
model_weight <- weight_cells$model_weight[match(cell_key, weight_key)]
weight_cells$sample_share <- weight_cells$sample_n / sum(weight_cells$sample_n)
weight_cells$target_share <- weight_cells$target_n / sum(weight_cells$target_n)
weight_cells$kish_effective_n <- (sum(model_weight)^2) / sum(model_weight^2)
write_table(weight_cells, file.path(weighted_root, "tables", "weight_diagnostics.csv"))

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
X <- as.matrix(cbind(
  athlete_covariates[rep(seq_len(J), 2L), ], domain_totals, palestra_z = palestra
))
discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
discipline <- as.matrix(raw[paste0("disciplina___", which(discipline_counts > 0))])
discipline <- discipline[rep(seq_len(J), 2L), , drop = FALSE]

# HSGP basis: reuse the exact configuration fit for the primary HSGP model
# (boundary L, number of basis functions M, share_center, lengthscale prior)
# rather than recomputing it, so the two analyses share one basis exactly.
hsgp_config <- read.csv(file.path(hsgp_root, "hsgp_configuration.csv"))
hsgp_L <- hsgp_config$boundary_L
hsgp_M <- hsgp_config$num_basis_functions
share_center <- hsgp_config$share_center
matern_nu <- hsgp_config$matern_nu
gp_rho_shape <- hsgp_config$gp_rho_prior_shape
gp_rho_rate <- hsgp_config$gp_rho_prior_rate
sqrt_lambda <- seq_len(hsgp_M) * pi / (2 * hsgp_L)
basis_gp <- function(values) sin(outer(values - share_center + hsgp_L, sqrt_lambda)) / sqrt(hsgp_L)
gp_basis <- matrix(0, length(share), hsgp_M)
gp_basis[any_carbon == 1, ] <- basis_gp(share[any_carbon == 1])

stopifnot(J > 0L, all(is.finite(model_weight)), all(model_weight > 0),
          abs(mean(model_weight) - 1) < 1e-12)
sampling_config <- if (quick) list(chains = 2L, iter = 600L, warmup = 300L) else
  list(chains = 4L, iter = 2000L, warmup = 1000L)
model <- rstan::stan_model(file = file.path(root, "weighted_model_hsgp.stan"))
fit <- rstan::sampling(
  model,
  data = list(
    N = length(y), J = J, P = ncol(X), D = ncol(discipline), M = hsgp_M,
    y = y, athlete = athlete, period = period,
    survey_weight = rep(model_weight, 2L), any_carbon = any_carbon,
    gp_basis = gp_basis, sqrt_lambda = sqrt_lambda, matern_nu = matern_nu,
    X = X, discipline = discipline,
    prior_exposure_sd = 0.7, prior_covariate_sd = 0.5, prior_gp_alpha_sd = 0.7,
    gp_rho_shape = gp_rho_shape, gp_rho_rate = gp_rho_rate
  ),
  chains = sampling_config$chains, iter = sampling_config$iter,
  warmup = sampling_config$warmup, seed = seed + 500L,
  refresh = if (quick) 0 else 200,
  control = list(adapt_delta = 0.98, max_treedepth = 12)
)

monitored <- c("alpha", "beta_any", "beta_share", "beta_common", "beta_discipline",
               "tau_discipline", "sigma_athlete", "gp_alpha", "gp_rho")
summary_matrix <- rstan::summary(fit, pars = monitored)$summary
sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
ebfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2) / var(x[, "energy__"]), numeric(1))
diagnostics <- data.frame(
  model = "age_sex_poststratified_hsgp", max_rhat = max(summary_matrix[, "Rhat"], na.rm = TRUE),
  min_bulk_proxy_neff = min(summary_matrix[, "n_eff"], na.rm = TRUE),
  divergences = sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1))),
  max_treedepth_hits = sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= 12), numeric(1))),
  min_ebfmi = min(ebfmi), chains = sampling_config$chains,
  iterations = sampling_config$iter, warmup = sampling_config$warmup
)
write_table(diagnostics, file.path(weighted_root, "tables", "sampling_diagnostics.csv"))
draws <- rstan::extract(fit, pars = monitored, permuted = TRUE)
dir.create(file.path(weighted_root, "models"), recursive = TRUE, showWarnings = FALSE)
saveRDS(draws, file.path(weighted_root, "models", "posterior_draws.rds"))
rm(fit)
gc()

gauss_hermite_normal <- function(n = 15L) {
  jacobi <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  list(nodes = sqrt(2) * eig$values, weights = eig$vectors[1, ]^2)
}
gh <- gauss_hermite_normal()
marginal_prevalence <- function(eta, sigma) {
  as.numeric(plogis(outer(eta, sigma * gh$nodes, "+")) %*% gh$weights)
}

groups <- list(
  Overall = seq_len(J), Men = which(sex == "Men"), Women = which(sex == "Women"),
  `Under 20` = which(age_group == "Under 20"),
  `Under 23` = which(age_group == "Under 23"), Senior = which(age_group == "Senior")
)
group_type <- c(Overall = "overall", Men = "sex", Women = "sex",
                `Under 20` = "age", `Under 23` = "age", Senior = "age")
samples <- nrow(draws$alpha)
summary_output <- list()
contrast_output <- list()
counter <- 1L

for (p in 1:2) {
  record_index <- seq_len(J) + (p - 1L) * J
  user_index <- record_index[any_carbon[record_index] == 1]
  values <- weighted_quantile(
    share[user_index], raw_weight[user_index - (p - 1L) * J], c(0.25, 0.5, 0.75)
  )
  names(values) <- c("low", "median", "high")
  basis_values <- basis_gp(values)
  rownames(basis_values) <- names(values)
  scenario <- array(
    NA_real_, dim = c(samples, length(groups), 4L),
    dimnames = list(NULL, names(groups), c("nonuser", "median_user", "low_user", "high_user"))
  )
  for (s in seq_len(samples)) {
    base <- draws$alpha[s, p] +
      as.numeric(X[record_index, , drop = FALSE] %*% draws$beta_common[s, ]) +
      as.numeric(discipline[record_index, , drop = FALSE] %*% draws$beta_discipline[s, ])
    probabilities <- matrix(
      NA_real_, J, 4L,
      dimnames = list(NULL, c("nonuser", "median_user", "low_user", "high_user"))
    )
    probabilities[, "nonuser"] <- marginal_prevalence(base, draws$sigma_athlete[s])
    for (name in c("median", "low", "high")) {
      exposure <- draws$beta_any[s, p] +
        sum(basis_values[name, ] * draws$beta_share[s, p, ])
      probabilities[, paste0(name, "_user")] <-
        marginal_prevalence(base + exposure, draws$sigma_athlete[s])
    }
    for (g in names(groups)) {
      index <- groups[[g]]
      scenario[s, g, ] <- colSums(probabilities[index, , drop = FALSE] * raw_weight[index]) /
        sum(raw_weight[index])
    }
  }
  period_name <- c("preparation", "competition")[p]
  for (g in names(groups)) {
    for (name in dimnames(scenario)[[3]]) {
      summary_output[[counter]] <- data.frame(
        group = g, group_type = group_type[g], period = period_name,
        estimand = paste0("prevalence_", name), measure = "prevalence",
        share_from = NA_real_,
        share_to = if (name == "nonuser") 0 else values[sub("_user", "", name)],
        t(summary_row(scenario[, g, name], "prevalence")), check.names = FALSE
      )
      counter <- counter + 1L
    }
    contrasts <- list(
      adoption_difference = scenario[, g, "median_user"] - scenario[, g, "nonuser"],
      adoption_ratio = scenario[, g, "median_user"] / scenario[, g, "nonuser"],
      dose_difference = scenario[, g, "high_user"] - scenario[, g, "low_user"],
      dose_ratio = scenario[, g, "high_user"] / scenario[, g, "low_user"]
    )
    for (name in names(contrasts)) {
      measure <- if (grepl("difference$", name)) "prevalence_difference" else "prevalence_ratio"
      dose <- grepl("^dose", name)
      summary_output[[counter]] <- data.frame(
        group = g, group_type = group_type[g], period = period_name,
        estimand = name, measure = measure,
        share_from = if (dose) values["low"] else 0,
        share_to = if (dose) values["high"] else values["median"],
        t(summary_row(contrasts[[name]], measure)), check.names = FALSE
      )
      contrast_output[[length(contrast_output) + 1L]] <- data.frame(
        draw = seq_len(samples), group = g, group_type = group_type[g],
        period = period_name, estimand = name, value = contrasts[[name]]
      )
      counter <- counter + 1L
    }
  }
}

estimands <- do.call(rbind, summary_output)
contrast_draws <- do.call(rbind, contrast_output)
write_table(estimands[estimands$group == "Overall", ],
            file.path(weighted_root, "tables", "weighted_association_estimands.csv"))
write_table(estimands[estimands$group_type == "sex", ],
            file.path(subgroup_root, "tables", "sex_conditional_estimands.csv"))
write_table(estimands[estimands$group_type == "age", ],
            file.path(subgroup_root, "tables", "age_group_conditional_estimands.csv"))
saveRDS(contrast_draws, file.path(subgroup_root, "tables", "conditional_contrast_draws.rds"))

unweighted <- read.csv(file.path(
  hsgp_root, "01_primary", "association", "tables", "association_estimands.csv"
))
unweighted <- unweighted[unweighted$measure %in% c("prevalence_difference", "prevalence_ratio"), ]
weighted <- estimands[estimands$group == "Overall" &
                         estimands$measure %in% c("prevalence_difference", "prevalence_ratio"), ]
comparison <- merge(
  unweighted[c("period", "estimand", "mean", "lower_95", "upper_95")],
  weighted[c("period", "estimand", "mean", "lower_95", "upper_95")],
  by = c("period", "estimand"), suffixes = c("_unweighted", "_weighted")
)
write_table(comparison, file.path(weighted_root, "tables", "weighted_vs_unweighted.csv"))

# Target-weighted fit check under observed exposures.
fit_check <- do.call(rbind, lapply(1:2, function(p) {
  record_index <- seq_len(J) + (p - 1L) * J
  predicted <- vapply(seq_len(samples), function(s) {
    eta <- draws$alpha[s, p] +
      as.numeric(X[record_index, , drop = FALSE] %*% draws$beta_common[s, ]) +
      as.numeric(discipline[record_index, , drop = FALSE] %*% draws$beta_discipline[s, ]) +
      any_carbon[record_index] * draws$beta_any[s, p] +
      rowSums(sweep(gp_basis[record_index, , drop = FALSE], 2L,
                    draws$beta_share[s, p, ], "*"))
    weighted.mean(marginal_prevalence(eta, draws$sigma_athlete[s]), raw_weight)
  }, numeric(1))
  data.frame(
    period = c("preparation", "competition")[p],
    observed_weighted_prevalence = weighted.mean(y[record_index], raw_weight),
    posterior_mean = mean(predicted), lower_95 = quantile(predicted, 0.025),
    upper_95 = quantile(predicted, 0.975)
  )
}))
write_table(fit_check, file.path(weighted_root, "tables", "weighted_fit_check.csv"))

# One-sided posterior densities for the requested conditional analyses.
plot_conditional <- function(type, filenames, title) {
  plot_data <- contrast_draws[
    contrast_draws$group_type == type & contrast_draws$period == "competition" &
      contrast_draws$estimand == "dose_difference",
  ]
  medians <- aggregate(value ~ group, plot_data, median)
  plot_data$group <- factor(plot_data$group, unique(plot_data$group))
  medians$group <- factor(medians$group, levels(plot_data$group))
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(value, fill = group)) +
    ggplot2::annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf,
                      fill = "grey70", alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey35") +
    ggplot2::geom_density(alpha = 0.82, color = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = medians, ggplot2::aes(value, 0), inherit.aes = FALSE,
                        shape = 21, fill = "white", size = 2.6, stroke = 0.7) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", switch = "y") +
    ggplot2::scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL,
                                expand = ggplot2::expansion(mult = c(0.04, 0.08))) +
    ggplot2::labs(
      title = title, subtitle = "Target-standardized Q75 vs Q25 carbon share among users (HSGP)",
      x = "Adjusted prevalence difference", y = NULL,
      caption = "Density is above baseline; dot is posterior median. Common coefficients: descriptive conditional estimates."
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "none", plot.title = ggplot2::element_text(face = "bold"),
      strip.placement = "outside", strip.background = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0)
    )
  for (path in filenames) {
    ggplot2::ggsave(path, plot, width = 9.5, height = if (type == "sex") 4.8 else 6,
                    dpi = 300, bg = "white")
  }
}
plot_conditional(
  "sex",
  c(file.path(root, "plots", "hsgp_weighted_competition_dose_by_sex.png"),
    file.path("paper", "figures", "fig_weighted_dose_by_sex.png")),
  "Competition association conditional on sex"
)
plot_conditional(
  "age",
  c(file.path(root, "plots", "hsgp_weighted_competition_dose_by_age.png"),
    file.path("paper", "figures", "fig_weighted_dose_by_age.png")),
  "Competition association conditional on age group"
)

write_table(data.frame(
  mode = if (quick) "quick" else "full", seed = seed,
  target_population = "2026 Italian championship frame estimate",
  poststrata = "sex x age category", target_n = sum(weight_cells$target_n),
  sample_n = J, kish_effective_n = unique(weight_cells$kish_effective_n),
  estimand_note = "Target-standardized descriptive associations; not causal or interaction estimates"
), file.path(weighted_root, "run_config.csv"))

if (!quick) stopifnot(diagnostics$max_rhat < 1.01, diagnostics$divergences == 0)
cat("HSGP weighted and conditional multilevel analyses complete.\n")
