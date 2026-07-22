options(stringsAsFactors = FALSE)
library(ggplot2)

set.seed(20260719)
root <- "experiments/multilevel"
plot_dir <- file.path(root, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")

z_score <- function(x) {
  s <- sd(x)
  as.numeric((x - mean(x)) / s)
}

# Reconstruct the exact primary-model design matrices (mirrors analysis.R).
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
stopifnot(J > 0L)
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
domains_out <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
domain_totals <- do.call(cbind, lapply(domains_out, function(domain) c(
  raw[[paste0("prep_carbon_", domain)]] + raw[[paste0("prep_no_carbon_", domain)]],
  raw[[paste0("gara_carbon_", domain)]] + raw[[paste0("gara_no_carbon_", domain)]]
)))
colnames(domain_totals) <- paste0(domains_out, "_total_z")
domain_totals <- apply(domain_totals, 2L, z_score)
palestra <- z_score(c(raw$prep_no_carbon_palestra, raw$gara_no_carbon_palestra))
X <- cbind(athlete_covariates[rep(seq_len(J), 2L), ], domain_totals, palestra_z = palestra)
X <- as.matrix(X); storage.mode(X) <- "double"

discipline_counts <- colSums(raw[paste0("disciplina___", 1:11)])
active_discipline_codes <- which(discipline_counts > 0)
discipline <- as.matrix(raw[paste0("disciplina___", active_discipline_codes)])
discipline <- discipline[rep(seq_len(J), 2L), , drop = FALSE]
storage.mode(discipline) <- "double"

share_center <- median(share[any_carbon == 1])
ns_fit <- splines::ns(share[any_carbon == 1], df = 3, Boundary.knots = c(0, 1))
ns_center <- as.numeric(predict(ns_fit, share_center))
share_basis <- matrix(0, length(share), 3L)
idx <- any_carbon == 1
share_basis[idx, ] <- sweep(predict(ns_fit, share[idx]), 2L, ns_center, "-")

# Posterior draws from the fitted primary model.
draws <- readRDS(file.path(root, "runs/01_primary/association/models/posterior_draws.rds"))
S <- nrow(draws$alpha)

prep_index <- period == 1L
comp_index <- period == 2L
replicated <- matrix(NA_real_, S, 2L, dimnames = list(NULL, c("Preparation", "Competition")))
for (s in seq_len(S)) {
  fixed <- draws$alpha[s, period] + as.numeric(X %*% draws$beta_common[s, ]) +
    as.numeric(discipline %*% draws$beta_discipline[s, ]) +
    any_carbon * draws$beta_any[s, period] +
    rowSums(share_basis * draws$beta_share[s, period, ])
  athlete_effect <- rnorm(J, 0, draws$sigma_athlete[s])
  eta <- fixed + athlete_effect[athlete]
  y_rep <- rbinom(length(eta), 1L, plogis(eta))
  replicated[s, "Preparation"] <- sum(y_rep[prep_index])
  replicated[s, "Competition"] <- sum(y_rep[comp_index])
}

observed <- c(Preparation = sum(prep_injury), Competition = sum(competition_injury))
plot_data <- data.frame(
  Period = factor(rep(colnames(replicated), each = S), c("Preparation", "Competition")),
  Count = as.vector(replicated)
)
observed_data <- data.frame(
  Period = factor(names(observed), c("Preparation", "Competition")), Count = observed
)

period_colors <- c(Preparation = "#0072B2", Competition = "#D55E00")
p <- ggplot(plot_data, aes(Count, fill = Period)) +
  geom_density(alpha = 0.55, color = "white", linewidth = 0.3) +
  geom_vline(data = observed_data, aes(xintercept = Count), linetype = 2, color = "grey15", linewidth = 0.8) +
  geom_text(data = observed_data, aes(x = Count, y = 0, label = paste("Observed:", Count)),
           vjust = -0.6, hjust = -0.08, color = "grey15", size = 3.6) +
  facet_wrap(~Period, scales = "free_x") +
  scale_fill_manual(values = period_colors, guide = "none") +
  labs(
    title = "Posterior predictive distribution of injured-athlete counts",
    subtitle = sprintf(
      "Replicated counts from the primary model (%d posterior draws) versus the observed count, by period", S
    ),
    x = "Number of athletes with injury", y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30")
  )

ggsave(file.path(plot_dir, "posterior_predictive_density.png"), p,
       width = 9.5, height = 5.2, dpi = 300, bg = "white")
cat("Saved posterior_predictive_density.png\n")
cat("Bayesian p (P(rep >= obs)): prep", mean(replicated[, "Preparation"] >= observed["Preparation"]),
   "comp", mean(replicated[, "Competition"] >= observed["Competition"]), "\n")
