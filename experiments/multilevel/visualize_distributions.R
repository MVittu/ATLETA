library(ggplot2)

root <- "experiments/multilevel"
plot_dir <- file.path(root, "plots")
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

z_score <- function(x) as.numeric((x - mean(x)) / sd(x))
raw_all <- read.csv(source_file, check.names = FALSE)
eligible <- raw_all$consenso == 1 & raw_all$eleggibile_calc == 1 &
  raw_all$survey_atleta_it_complete == 2
raw <- raw_all[eligible, , drop = FALSE]
sum_named <- function(prefix) rowSums(raw[paste0(prefix, domains)])
prep_carbon <- sum_named("prep_carbon_")
prep_noncarbon <- sum_named("prep_no_carbon_")
competition_carbon <- sum_named("gara_carbon_")
competition_noncarbon <- sum_named("gara_no_carbon_")
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

J <- nrow(raw)
period <- c(rep(1L, J), rep(2L, J))
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
X <- as.matrix(cbind(
  athlete_covariates[rep(seq_len(J), 2L), ], domain_totals, palestra_z = palestra
))

discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
discipline_athlete <- as.matrix(raw[paste0("disciplina___", which(discipline_counts > 0))])
discipline_profile <- apply(discipline_athlete, 1L, paste0, collapse = "")
reference_discipline <- discipline_athlete[
  match(names(which.max(table(discipline_profile))), discipline_profile), ,
  drop = FALSE
]
discipline <- discipline_athlete[rep(seq_len(J), 2L), , drop = FALSE]

share_center <- median(share[any_carbon == 1])
ns_fit <- splines::ns(share[any_carbon == 1], df = 3, Boundary.knots = c(0, 1))
ns_center <- as.numeric(predict(ns_fit, share_center))
basis_spline <- function(values) sweep(predict(ns_fit, values), 2L, ns_center, "-")
basis_linear <- function(values) matrix(values - share_center, ncol = 1L)

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

derive_contrasts <- function(draws, basis_function, model_name) {
  samples <- nrow(draws$alpha)
  output <- vector("list", 2L)
  for (p in 1:2) {
    index <- which(period == p)
    user_share <- share[index][any_carbon[index] == 1]
    values <- c(low = unname(quantile(user_share, 0.25)),
                median = median(user_share), high = unname(quantile(user_share, 0.75)))
    basis <- basis_function(values)
    rownames(basis) <- names(values)
    scenario <- matrix(NA_real_, samples, 4L,
                       dimnames = list(NULL, c("nonuser", "median", "low", "high")))
    for (s in seq_len(samples)) {
      base <- draws$alpha[s, p] +
        as.numeric(X[index, , drop = FALSE] %*% draws$beta_common[s, ]) +
        as.numeric(discipline[index, , drop = FALSE] %*% draws$beta_discipline[s, ])
      scenario[s, "nonuser"] <- mean(marginal_prevalence(base, draws$sigma_athlete[s]))
      for (name in c("median", "low", "high")) {
        exposure <- draws$beta_any[s, p] + sum(basis[name, ] * draws$beta_share[s, p, ])
        scenario[s, name] <- mean(marginal_prevalence(base + exposure, draws$sigma_athlete[s]))
      }
    }
    output[[p]] <- rbind(
      data.frame(model = model_name, draw = seq_len(samples), period = p,
                 estimand = "adoption", value = scenario[, "median"] - scenario[, "nonuser"]),
      data.frame(model = model_name, draw = seq_len(samples), period = p,
                 estimand = "higher share", value = scenario[, "high"] - scenario[, "low"])
    )
  }
  do.call(rbind, output)
}

model_files <- c(
  primary = "runs/01_primary/association/models/posterior_draws.rds",
  linear = "runs/02_sensitivity/linear_share/models/posterior_draws.rds",
  wide = "runs/02_sensitivity/wide_prior/models/posterior_draws.rds"
)
contrast_draws <- do.call(rbind, lapply(names(model_files), function(name) {
  derive_contrasts(
    readRDS(file.path(root, model_files[[name]])),
    if (name == "linear") basis_linear else basis_spline,
    name
  )
}))
stopifnot(nrow(contrast_draws) == 3L * 2L * 2L * nrow(
  readRDS(file.path(root, model_files[["primary"]]))$alpha
))

period_colors <- c(Preparation = "#0072B2", Competition = "#D55E00")
theme_results <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"), legend.position = "top"
  )

# Profile-specific predictions from the existing primary model.
draws <- readRDS(file.path(root, model_files[["primary"]]))
reference_x <- apply(X, 2L, median)
reference_base <- sapply(1:2, function(p) {
  draws$alpha[, p] +
    as.numeric(draws$beta_common %*% reference_x) +
    as.numeric(draws$beta_discipline %*% as.numeric(reference_discipline))
})
adoption_draws <- do.call(rbind, lapply(1:2, function(p) {
  nonuser <- plogis(reference_base[, p])
  user <- plogis(reference_base[, p] + draws$beta_any[, p])
  data.frame(
    Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
    nonuser = nonuser, user = user, difference = user - nonuser
  )
}))
adoption_summary <- do.call(rbind, lapply(split(adoption_draws, adoption_draws$Period), function(x) {
  do.call(rbind, lapply(c("nonuser", "user", "difference"), function(name) data.frame(
    Period = x$Period[1], estimand = name, median = median(x[[name]]),
    lower = unname(quantile(x[[name]], 0.025)),
    upper = unname(quantile(x[[name]], 0.975))
  )))
}))
adoption_difference <- adoption_summary[adoption_summary$estimand == "difference", ]

p_posterior <- ggplot(adoption_draws, aes(Period, difference, fill = Period)) +
  geom_hline(yintercept = 0, color = "grey35", linetype = 2) +
  geom_violin(trim = FALSE, alpha = 0.78, color = "white", linewidth = 0.4) +
  geom_pointrange(
    data = adoption_difference,
    aes(y = median, ymin = lower, ymax = upper),
    color = "black", fill = "white", shape = 21, linewidth = 0.6
  ) +
  scale_fill_manual(values = period_colors, guide = "none") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Adoption and profile-specific injury risk",
    subtitle = sprintf(
      "User at %.0f%% carbon share minus nonuser, with all other predictors fixed",
      100 * share_center
    ),
    x = NULL, y = "Difference in predicted injury risk",
    caption = "Violin: posterior distribution; dot and line: median and 95% credible interval.\nProfile: cohort medians, modal discipline indicators, athlete random effect = 0."
  ) +
  theme_results
ggsave(file.path(plot_dir, "posterior_association_distributions.png"), p_posterior,
       width = 8, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path("paper", "figures", "fig_posterior_association_distributions.png"), p_posterior,
       width = 8, height = 5.5, dpi = 300, bg = "white")

# Exact posterior distribution of the key result under all specifications.
robust <- contrast_draws[
  contrast_draws$period == 2 & contrast_draws$estimand == "higher share",
]
robust$Model <- factor(
  c(primary = "Primary nonlinear share", linear = "Linear share", wide = "Wider priors")[robust$model],
  c("Primary nonlinear share", "Linear share", "Wider priors")
)
robust_medians <- aggregate(value ~ Model, robust, median)
p_sensitivity <- ggplot(robust, aes(value, fill = Model)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  geom_vline(xintercept = 0, color = "grey35", linetype = 2) +
  geom_density(alpha = 0.82, color = "white", linewidth = 0.3) +
  geom_point(data = robust_medians, aes(value, 0), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2.6, stroke = 0.7) +
  facet_grid(Model ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(
    "Primary nonlinear share" = "#0072B2", "Linear share" = "#E69F00",
    "Wider priors" = "#6A3D9A"
  ), guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Full posterior distribution of the competition association",
    subtitle = "Q75 vs Q25 carbon share among competition-period users",
    x = "Adjusted prevalence difference", y = NULL,
    caption = "Density is shown above the baseline; dot is the posterior median."
  ) +
  theme_results +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0), panel.spacing.y = grid::unit(0.3, "lines"))
ggsave(file.path(plot_dir, "posterior_competition_sensitivity.png"), p_sensitivity,
       width = 9.5, height = 5.4, dpi = 300, bg = "white")

# Empirical exposure support, including the mass of nonusers at zero.
observed <- data.frame(
  Period = factor(rep(c("Preparation", "Competition"), each = J),
                  c("Preparation", "Competition")),
  share = share
)
p_observed <- ggplot(observed, aes(share, fill = Period)) +
  geom_histogram(binwidth = 0.05, boundary = 0, closed = "left",
                 color = "white", linewidth = 0.2) +
  facet_wrap(~Period, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = period_colors, guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Observed carbon-frequency share",
    subtitle = sprintf(
      "All %d athletes; zero marks nonuse and one means all reported type-frequency was carbon", J
    ),
    x = "Carbon share of reported weekly type-frequency counts", y = "Athletes",
    caption = "Frequency items are 0, 1, 2, 3, 4, or 5+ sessions/week; the dataset encodes 5+ as 5."
  ) +
  theme_results
ggsave(file.path(plot_dir, "observed_carbon_share_distribution.png"), p_observed,
       width = 9.5, height = 7, dpi = 300, bg = "white")

# Absolute fitted risk among users with both periods shown on the same reference profile.
user_ranges <- tapply(share[any_carbon == 1], period[any_carbon == 1], range)
share_grid <- seq(max(vapply(user_ranges, `[`, numeric(1), 1L)),
                  min(vapply(user_ranges, `[`, numeric(1), 2L)), length.out = 121)
profile_risk_draws <- function(values, p) {
  eta <- sweep(basis_spline(values) %*% t(draws$beta_share[, p, ]), 2L,
               reference_base[, p] + draws$beta_any[, p], "+")
  plogis(eta)
}
profile_risk <- function(values, p) {
  risk <- profile_risk_draws(values, p)
  data.frame(
    Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
    share = values,
    median = apply(risk, 1, median),
    lower = apply(risk, 1, quantile, 0.025),
    upper = apply(risk, 1, quantile, 0.975)
  )
}
spline_curve <- do.call(rbind, lapply(1:2, function(p) profile_risk(share_grid, p)))
selected_risks <- do.call(rbind, lapply(1:2, function(p) profile_risk(c(0.25, 0.5, 0.75), p)))
share_contrasts <- do.call(rbind, lapply(1:2, function(p) {
  risk <- profile_risk_draws(c(0.25, 0.75), p)
  difference <- risk[2, ] - risk[1, ]
  data.frame(
    contrast = "75% minus 25%", Period = c("Preparation", "Competition")[p],
    median = median(difference), lower = unname(quantile(difference, 0.025)),
    upper = unname(quantile(difference, 0.975))
  )
}))
phase_contrasts <- do.call(rbind, lapply(c(0.25, 0.5, 0.75), function(value) {
  difference <- profile_risk_draws(value, 2) - profile_risk_draws(value, 1)
  data.frame(
    contrast = "competition minus preparation", share = value,
    median = median(difference), lower = unname(quantile(difference, 0.025)),
    upper = unname(quantile(difference, 0.975))
  )
}))
user_support <- observed[observed$share > 0, ]
stopifnot(all(spline_curve$lower <= spline_curve$median),
          all(spline_curve$median <= spline_curve$upper))

profile_dir <- file.path(root, "runs", "05_profile_predictions", "tables")
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(adoption_summary, file.path(profile_dir, "adoption_risk.csv"), row.names = FALSE)
write.csv(selected_risks, file.path(profile_dir, "share_risk_selected.csv"), row.names = FALSE)
write.csv(share_contrasts, file.path(profile_dir, "share_risk_contrast.csv"), row.names = FALSE)
write.csv(phase_contrasts, file.path(profile_dir, "period_risk_contrast.csv"), row.names = FALSE)
write.csv(spline_curve, file.path(profile_dir, "share_risk_curve.csv"), row.names = FALSE)
write.csv(rbind(
  data.frame(type = "covariate", term = colnames(X), value = reference_x),
  data.frame(type = "discipline_indicator", term = colnames(discipline_athlete),
             value = as.numeric(reference_discipline)),
  data.frame(type = "random_effect", term = "athlete", value = 0)
), file.path(profile_dir, "reference_profile.csv"), row.names = FALSE)

p_spline <- ggplot(spline_curve, aes(share, median, color = Period, fill = Period)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = user_support, aes(x = share, color = Period), inherit.aes = FALSE,
           sides = "b", alpha = 0.16, length = grid::unit(0.035, "npc")) +
  scale_color_manual(values = period_colors) +
  scale_fill_manual(values = period_colors) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Profile-specific injury risk across carbon share",
    subtitle = "Preparation and competition predictions for the same reference profile",
    x = "Carbon share of reported weekly type-frequency counts", y = "Predicted injury risk",
    caption = "Natural-spline posterior median and 95% credible interval; rug: observed shares.\nProfile: cohort medians, modal discipline indicators, athlete random effect = 0."
  ) +
  theme_results
ggsave(file.path(plot_dir, "spline_dose_response.png"), p_spline,
       width = 9.5, height = 7, dpi = 300, bg = "white")
ggsave(file.path("paper", "figures", "fig_spline_dose_response.png"), p_spline,
       width = 9.5, height = 7, dpi = 300, bg = "white")

expected <- file.path(plot_dir, c(
  "posterior_association_distributions.png",
  "posterior_competition_sensitivity.png",
  "observed_carbon_share_distribution.png",
  "spline_dose_response.png"
))
stopifnot(all(file.exists(expected)), all(file.info(expected)$size > 1000))
cat("Saved 4 full-distribution figures to", plot_dir, "\n")
