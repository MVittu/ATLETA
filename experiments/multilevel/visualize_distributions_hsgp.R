options(stringsAsFactors = FALSE)
library(ggplot2)

# GP-based counterparts of visualize_distributions.R / plot_adoption_posterior_density.R /
# plot_share_interpretive_anchors.R / visualize_ppc_density.R, reading from runs/07_hsgp
# instead of runs/01_primary and runs/02_sensitivity. Overwrites the same paper/figures
# filenames (fig_posterior_association_distributions.png, fig_dose_response.png --
# renamed from fig_spline_dose_response.png since the curve is no longer a spline --
# fig_competition_dose_sensitivity.png, fig_posterior_predictive.png) so the manuscript
# figures reflect the HSGP primary analysis. Does not touch the original spline-based
# scripts or their run outputs.

set.seed(20260719)
root <- "experiments/multilevel"
plot_dir <- file.path(root, "plots")
run_root <- file.path(root, "runs", "07_hsgp")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")

z_score <- function(x) as.numeric((x - mean(x)) / sd(x))

# ---- Data (mirrors analysis_hsgp.R) ------------------------------------------------

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
prep_carbon <- prep_carbon[valid]; prep_noncarbon <- prep_noncarbon[valid]
competition_carbon <- competition_carbon[valid]; competition_noncarbon <- competition_noncarbon[valid]
prep_injury <- prep_injury[valid]; competition_injury <- competition_injury[valid]

J <- nrow(raw)
stopifnot(J == 352L)
period <- c(rep(1L, J), rep(2L, J))
athlete <- rep(seq_len(J), 2L)
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
X <- as.matrix(cbind(athlete_covariates[rep(seq_len(J), 2L), ], domain_totals, palestra_z = palestra))
storage.mode(X) <- "double"

discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
active_discipline_codes <- which(discipline_counts > 0)
discipline_athlete <- as.matrix(raw[paste0("disciplina___", active_discipline_codes)])
discipline_profile <- apply(discipline_athlete, 1L, paste0, collapse = "")
reference_discipline <- as.numeric(discipline_athlete[
  match(names(which.max(table(discipline_profile))), discipline_profile), , drop = FALSE
])
discipline <- discipline_athlete[rep(seq_len(J), 2L), , drop = FALSE]
storage.mode(discipline) <- "double"
reference_x <- apply(X, 2L, median)

# ---- HSGP basis (must match analysis_hsgp.R exactly; read back the fitted config) ----

hsgp_config <- read.csv(file.path(run_root, "hsgp_configuration.csv"))
hsgp_L <- hsgp_config$boundary_L
hsgp_M <- hsgp_config$num_basis_functions
share_center <- hsgp_config$share_center
sqrt_lambda <- seq_len(hsgp_M) * pi / (2 * hsgp_L)

basis_gp <- function(values, used = rep(1, length(values))) {
  out <- matrix(0, length(values), hsgp_M)
  index <- used == 1
  if (any(index)) out[index, ] <- sin(outer(values[index] - share_center + hsgp_L, sqrt_lambda)) / sqrt(hsgp_L)
  out
}
basis_linear <- function(values, used = rep(1, length(values))) {
  out <- matrix(0, length(values), 1L)
  out[used == 1, 1] <- values[used == 1] - share_center
  out
}

gauss_hermite_normal <- function(n = 15L) {
  jacobi <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1L) / 2)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off
  eig <- eigen(jacobi, symmetric = TRUE)
  list(nodes = sqrt(2) * eig$values, weights = eig$vectors[1, ]^2)
}
gh <- gauss_hermite_normal()
marginal_prevalence <- function(eta, sigma) as.numeric(plogis(outer(eta, sigma * gh$nodes, "+")) %*% gh$weights)

period_colors <- c(Preparation = "#0072B2", Competition = "#D55E00")
theme_results <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"), legend.position = "top"
  )

# ---- Figure 1: posterior distribution of the adoption contrast (fixed profile) -------

draws <- readRDS(file.path(run_root, "01_primary", "association", "models", "posterior_draws.rds"))
reference_base <- sapply(1:2, function(p) {
  draws$alpha[, p] + as.numeric(draws$beta_common %*% reference_x) +
    as.numeric(draws$beta_discipline %*% reference_discipline)
})
adoption_draws <- do.call(rbind, lapply(1:2, function(p) {
  data.frame(
    Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
    difference = plogis(reference_base[, p] + draws$beta_any[, p]) - plogis(reference_base[, p])
  )
}))
adoption_summary <- do.call(rbind, lapply(split(adoption_draws$difference, adoption_draws$Period), function(x) {
  data.frame(median = median(x), lower = unname(quantile(x, 0.025)), upper = unname(quantile(x, 0.975)))
}))
adoption_summary$Period <- factor(rownames(adoption_summary), c("Preparation", "Competition"))

p_posterior <- ggplot(adoption_draws, aes(difference, fill = Period)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.22) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45, linetype = 2) +
  geom_density(alpha = 0.82, colour = "white", linewidth = 0.35) +
  geom_segment(data = adoption_summary, aes(x = lower, xend = upper, y = 0, yend = 0),
               inherit.aes = FALSE, linewidth = 0.9) +
  geom_point(data = adoption_summary, aes(x = median, y = 0), inherit.aes = FALSE,
             shape = 21, fill = "white", size = 2.6, stroke = 0.75) +
  facet_grid(Period ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = period_colors, guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Posterior distribution of the adoption contrast",
    subtitle = sprintf("User at exactly %.0f%% carbon share minus non-user, for the same fixed profile", 100 * share_center),
    x = "Difference in predicted injury risk", y = NULL,
    caption = "Density: posterior distribution; dot and segment: median and 95% credible interval; grey band: ±2-pp ROPE."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    strip.placement = "outside", strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "#444444"),
    plot.caption = element_text(hjust = 0), panel.spacing.y = grid::unit(0.35, "lines")
  )
for (path in c(file.path(plot_dir, "hsgp_posterior_association_distributions.png"),
               file.path("paper", "figures", "fig_posterior_association_distributions.png"))) {
  ggsave(path, p_posterior, width = 8, height = 4.8, dpi = 320, bg = "white")
}

# ---- Figure 2: population-standardized curve (with posterior draws) and ------------
# ---- fixed-profile curve, stacked so the primary, population-level estimand --------
# ---- is the headline panel and the single-athlete reading is secondary. -----------

user_ranges <- tapply(share[any_carbon == 1], period[any_carbon == 1], range)
share_grid <- seq(max(vapply(user_ranges, `[`, numeric(1), 1L)),
                  min(vapply(user_ranges, `[`, numeric(1), 2L)), length.out = 121)
pop_share_grid <- seq(max(vapply(user_ranges, `[`, numeric(1), 1L)),
                      min(vapply(user_ranges, `[`, numeric(1), 2L)), length.out = 41)

profile_risk_draws <- function(values, p) {
  eta <- sweep(basis_gp(values) %*% t(draws$beta_share[, p, ]), 2L,
               reference_base[, p] + draws$beta_any[, p], "+")
  plogis(eta)
}
profile_risk <- function(values, p) {
  risk <- profile_risk_draws(values, p)
  data.frame(
    Type = "Fixed profile",
    Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
    share = values, median = apply(risk, 1, median),
    lower = apply(risk, 1, quantile, 0.025), upper = apply(risk, 1, quantile, 0.975)
  )
}
profile_curve <- do.call(rbind, lapply(1:2, function(p) profile_risk(share_grid, p)))
profile_anchors <- do.call(rbind, lapply(1:2, function(p) profile_risk(c(0.25, 0.5, 0.75), p)))

# Population-standardized: for a grid of counterfactual share values, average
# predicted prevalence over every athlete's observed covariates/discipline in
# that period (same estimand as the P25/P50/P75 table rows, read continuously).
# A 1000-draw subsample keeps this affordable; a further 60-draw subsample of
# individual curves ("spaghetti") shows how much curve *shape*, not just
# pointwise level, the posterior leaves open -- the ribbon alone cannot show
# that the same draw is smooth or steep across its whole length.
S <- nrow(draws$alpha)
ribbon_sub <- sample(seq_len(S), min(1000L, S))
spaghetti_sub <- sample(ribbon_sub, 60L)
spaghetti_rows <- match(spaghetti_sub, ribbon_sub)

population_curve_draws <- function(values, p, sub) {
  index <- which(period == p)
  Xp <- X[index, , drop = FALSE]
  Dp <- discipline[index, , drop = FALSE]
  basis_grid <- basis_gp(values)
  out <- matrix(NA_real_, length(sub), length(values))
  for (row in seq_along(sub)) {
    s <- sub[row]
    base <- draws$alpha[s, p] + as.numeric(Xp %*% draws$beta_common[s, ]) +
      as.numeric(Dp %*% draws$beta_discipline[s, ])
    exposure <- draws$beta_any[s, p] + as.numeric(basis_grid %*% draws$beta_share[s, p, ])
    for (i in seq_along(values)) {
      out[row, i] <- mean(marginal_prevalence(base + exposure[i], draws$sigma_athlete[s]))
    }
  }
  out
}
population_summary <- function(risk, values, p) data.frame(
  Type = "Population-standardized",
  Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
  share = values, median = apply(risk, 2, median),
  lower = apply(risk, 2, quantile, 0.025), upper = apply(risk, 2, quantile, 0.975)
)

pop_curve <- list(); pop_anchor_list <- list(); spaghetti_list <- list()
for (p in 1:2) {
  risk_grid <- population_curve_draws(pop_share_grid, p, ribbon_sub)
  pop_curve[[p]] <- population_summary(risk_grid, pop_share_grid, p)
  risk_anchor <- population_curve_draws(c(0.25, 0.5, 0.75), p, ribbon_sub)
  pop_anchor_list[[p]] <- population_summary(risk_anchor, c(0.25, 0.5, 0.75), p)
  spaghetti_list[[p]] <- do.call(rbind, lapply(seq_along(spaghetti_sub), function(row) data.frame(
    Type = "Population-standardized",
    Period = factor(c("Preparation", "Competition")[p], c("Preparation", "Competition")),
    draw = spaghetti_sub[row], share = pop_share_grid,
    value = risk_grid[spaghetti_rows[row], ]
  )))
}
pop_curve <- do.call(rbind, pop_curve)
pop_anchors <- do.call(rbind, pop_anchor_list)
spaghetti <- do.call(rbind, spaghetti_list)

curve_data <- rbind(pop_curve, profile_curve)
curve_data$Type <- factor(curve_data$Type, c("Population-standardized", "Fixed profile"))
anchor_data <- rbind(pop_anchors, profile_anchors)
anchor_data$Type <- factor(anchor_data$Type, levels(curve_data$Type))
anchor_data$label <- sprintf("%.1f%%", 100 * anchor_data$median)
spaghetti$Type <- factor(spaghetti$Type, levels(curve_data$Type))
stopifnot(all(curve_data$lower <= curve_data$median), all(curve_data$median <= curve_data$upper))

rug <- data.frame(
  Period = factor(rep(c("Preparation", "Competition"), each = J), c("Preparation", "Competition")),
  share = share
)
rug <- rug[rug$share > 0, ]

p_gp <- ggplot(curve_data, aes(share, median, colour = Period, fill = Period)) +
  geom_line(data = spaghetti, aes(share, value, group = interaction(Period, draw)),
            inherit.aes = FALSE, colour = "grey40", alpha = 0.07, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_vline(xintercept = c(0.25, 0.50, 0.75), colour = "grey55", linewidth = 0.4, linetype = 3) +
  geom_line(linewidth = 1.1) +
  geom_errorbar(data = anchor_data, aes(ymin = lower, ymax = upper), width = 0.012, linewidth = 0.65) +
  geom_point(data = anchor_data, shape = 21, fill = "white", size = 2.5, stroke = 0.8) +
  geom_label(data = anchor_data, aes(label = label), nudge_x = -0.035, nudge_y = 0.035,
             fill = "white", linewidth = 0, label.padding = grid::unit(0.08, "lines"),
             show.legend = FALSE, size = 3.1) +
  geom_rug(data = rug, aes(x = share, colour = Period), inherit.aes = FALSE,
           sides = "b", alpha = 0.16, length = grid::unit(0.035, "npc")) +
  facet_grid(Type ~ Period) +
  scale_colour_manual(values = period_colors, guide = "none") +
  scale_fill_manual(values = period_colors, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1), labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Injury risk across carbon share: population-standardized (top) and fixed-profile (bottom)",
    subtitle = "Thin grey lines: 60 posterior draws of the top curve, showing shape uncertainty the ribbon alone hides",
    x = "Carbon share of reported training-frequency score", y = "Predicted injury risk",
    caption = paste(
      "25%, 50%, and 75% are exact carbon-share values, not percentiles, bins, or clinical thresholds.",
      "Bars are 95% credible intervals; rug marks observed user shares.",
      sep = "\n"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(), strip.text = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "#444444"), plot.caption = element_text(hjust = 0),
    panel.spacing = grid::unit(1.1, "lines")
  )
for (path in c(file.path(plot_dir, "hsgp_dose_response.png"),
               file.path("paper", "figures", "fig_dose_response.png"))) {
  ggsave(path, p_gp, width = 8.5, height = 9, dpi = 320, bg = "white")
}

# ---- Figure 3: full posterior density of the competition dose contrast, 3 specs ------

derive_competition_dose <- function(draws, basis_function) {
  p <- 2L
  index <- which(period == p)
  user_share <- share[index][any_carbon[index] == 1]
  values <- c(low = unname(quantile(user_share, 0.25)), high = unname(quantile(user_share, 0.75)))
  basis <- basis_function(values)
  rownames(basis) <- names(values)
  samples <- nrow(draws$alpha)
  out <- numeric(samples)
  for (s in seq_len(samples)) {
    base <- draws$alpha[s, p] + as.numeric(X[index, , drop = FALSE] %*% draws$beta_common[s, ]) +
      as.numeric(discipline[index, , drop = FALSE] %*% draws$beta_discipline[s, ])
    low <- mean(marginal_prevalence(base + draws$beta_any[s, p] + sum(basis["low", ] * draws$beta_share[s, p, ]), draws$sigma_athlete[s]))
    high <- mean(marginal_prevalence(base + draws$beta_any[s, p] + sum(basis["high", ] * draws$beta_share[s, p, ]), draws$sigma_athlete[s]))
    out[s] <- high - low
  }
  out
}

draws_linear <- readRDS(file.path(run_root, "02_sensitivity", "linear_share", "models", "posterior_draws.rds"))
draws_wide <- readRDS(file.path(run_root, "02_sensitivity", "wide_prior", "models", "posterior_draws.rds"))
robust <- rbind(
  data.frame(Model = "Primary (HSGP share)", value = derive_competition_dose(draws, basis_gp)),
  data.frame(Model = "Linear share", value = derive_competition_dose(draws_linear, basis_linear)),
  data.frame(Model = "Wider GP priors", value = derive_competition_dose(draws_wide, basis_gp))
)
robust$Model <- factor(robust$Model, c("Primary (HSGP share)", "Linear share", "Wider GP priors"))
robust_medians <- aggregate(value ~ Model, robust, median)

p_sensitivity <- ggplot(robust, aes(value, fill = Model)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.2) +
  geom_vline(xintercept = 0, color = "grey35", linetype = 2) +
  geom_density(alpha = 0.82, color = "white", linewidth = 0.3) +
  geom_point(data = robust_medians, aes(value, 0), inherit.aes = FALSE,
             shape = 21, fill = "white", color = "black", size = 2.6, stroke = 0.7) +
  facet_grid(Model ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(
    "Primary (HSGP share)" = "#0072B2", "Linear share" = "#E69F00", "Wider GP priors" = "#6A3D9A"
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
for (path in c(file.path(plot_dir, "hsgp_competition_dose_sensitivity.png"),
               file.path("paper", "figures", "fig_competition_dose_sensitivity.png"))) {
  ggsave(path, p_sensitivity, width = 9.5, height = 5.4, dpi = 320, bg = "white")
}

# ---- Figure 4: posterior predictive density of injured-athlete counts ----------------

gp_basis_full <- basis_gp(share, any_carbon)
S <- nrow(draws$alpha)
prep_index <- period == 1L
comp_index <- period == 2L
replicated <- matrix(NA_real_, S, 2L, dimnames = list(NULL, c("Preparation", "Competition")))
for (s in seq_len(S)) {
  fixed <- draws$alpha[s, period] + as.numeric(X %*% draws$beta_common[s, ]) +
    as.numeric(discipline %*% draws$beta_discipline[s, ]) +
    any_carbon * draws$beta_any[s, period] +
    rowSums(gp_basis_full * draws$beta_share[s, period, ])
  athlete_effect <- rnorm(J, 0, draws$sigma_athlete[s])
  eta <- fixed + athlete_effect[athlete]
  y_rep <- rbinom(length(eta), 1L, plogis(eta))
  replicated[s, "Preparation"] <- sum(y_rep[prep_index])
  replicated[s, "Competition"] <- sum(y_rep[comp_index])
}
observed <- c(Preparation = sum(prep_injury), Competition = sum(competition_injury))
ppc_data <- data.frame(
  Period = factor(rep(colnames(replicated), each = S), c("Preparation", "Competition")),
  Count = as.vector(replicated)
)
observed_data <- data.frame(Period = factor(names(observed), c("Preparation", "Competition")), Count = observed)

p_ppc <- ggplot(ppc_data, aes(Count, fill = Period)) +
  geom_density(alpha = 0.55, color = "white", linewidth = 0.3) +
  geom_vline(data = observed_data, aes(xintercept = Count), linetype = 2, color = "grey15", linewidth = 0.8) +
  geom_text(data = observed_data, aes(x = Count, y = 0, label = paste("Observed:", Count)),
           vjust = -0.6, hjust = -0.08, color = "grey15", size = 3.6) +
  facet_wrap(~Period, scales = "free_x") +
  scale_fill_manual(values = period_colors, guide = "none") +
  labs(
    title = "Posterior predictive distribution of injured-athlete counts",
    subtitle = sprintf(
      "Replicated counts from the primary HSGP model (%d posterior draws) versus the observed count, by period", S
    ),
    x = "Number of athletes with injury", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )
for (path in c(file.path(plot_dir, "hsgp_posterior_predictive_density.png"),
               file.path("paper", "figures", "fig_posterior_predictive.png"))) {
  ggsave(path, p_ppc, width = 9.5, height = 5.2, dpi = 300, bg = "white")
}

cat("Bayesian p (P(rep >= obs)): prep", mean(replicated[, "Preparation"] >= observed["Preparation"]),
   "comp", mean(replicated[, "Competition"] >= observed["Competition"]), "\n")

expected <- file.path("paper", "figures", c(
  "fig_posterior_association_distributions.png", "fig_dose_response.png",
  "fig_competition_dose_sensitivity.png", "fig_posterior_predictive.png"
))
stopifnot(all(file.exists(expected)), all(file.info(expected)$size > 1000))
cat("Saved 4 HSGP figures to paper/figures\n")
