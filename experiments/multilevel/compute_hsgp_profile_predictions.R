options(stringsAsFactors = FALSE)

# One-off companion to visualize_distributions_hsgp.R: saves the same
# profile-prediction tables the original pipeline wrote to
# runs/05_profile_predictions/tables, but for the HSGP primary model, so the
# manuscript text can cite precise, reproducible numbers instead of reading
# them off a plot.

root <- "experiments/multilevel"
run_root <- file.path(root, "runs", "07_hsgp")
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
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

J <- nrow(raw)
period <- c(rep(1L, J), rep(2L, J))
carbon <- c(prep_carbon, competition_carbon)
noncarbon <- c(prep_noncarbon, competition_noncarbon)
share <- carbon / (carbon + noncarbon)
any_carbon <- as.numeric(carbon > 0)

height <- as.numeric(raw$statura_cm)
height[height > 0 & height < 3] <- height[height > 0 & height < 3] * 100
athlete_covariates <- data.frame(
  age_z = z_score(as.numeric(raw$eta)), sex_female = as.numeric(raw$genere_gara == 2),
  height_z = z_score(height), weight_z = z_score(as.numeric(raw$peso_kg)),
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
reference_x <- apply(X, 2L, median)

hsgp_config <- read.csv(file.path(run_root, "hsgp_configuration.csv"))
hsgp_L <- hsgp_config$boundary_L
hsgp_M <- hsgp_config$num_basis_functions
share_center <- hsgp_config$share_center
sqrt_lambda <- seq_len(hsgp_M) * pi / (2 * hsgp_L)
basis_gp <- function(values) sin(outer(values - share_center + hsgp_L, sqrt_lambda)) / sqrt(hsgp_L)

draws <- readRDS(file.path(run_root, "01_primary", "association", "models", "posterior_draws.rds"))
reference_base <- sapply(1:2, function(p) {
  draws$alpha[, p] + as.numeric(draws$beta_common %*% reference_x) +
    as.numeric(draws$beta_discipline %*% reference_discipline)
})

summarize <- function(x) c(median = median(x), lower = unname(quantile(x, 0.025)), upper = unname(quantile(x, 0.975)), p_dir = NA_real_)
summarize_contrast <- function(x) c(
  median = median(x), lower = unname(quantile(x, 0.025)), upper = unname(quantile(x, 0.975)),
  p_dir = max(mean(x > 0), mean(x < 0))
)

# Adoption at exactly 50% share (fixed profile).
adoption_rows <- do.call(rbind, lapply(1:2, function(p) {
  nonuser <- plogis(reference_base[, p])
  user <- plogis(reference_base[, p] + draws$beta_any[, p])
  rbind(
    data.frame(period = c("preparation", "competition")[p], estimand = "nonuser", t(summarize(nonuser))),
    data.frame(period = c("preparation", "competition")[p], estimand = "user", t(summarize(user))),
    data.frame(period = c("preparation", "competition")[p], estimand = "difference", t(summarize_contrast(user - nonuser)))
  )
}))
rownames(adoption_rows) <- NULL

# Continuous curve read at exactly 25/50/75% share (fixed profile).
profile_risk_draws <- function(values, p) {
  eta <- sweep(basis_gp(values) %*% t(draws$beta_share[, p, ]), 2L,
               reference_base[, p] + draws$beta_any[, p], "+")
  plogis(eta)
}
selected_rows <- do.call(rbind, lapply(1:2, function(p) {
  risk <- profile_risk_draws(c(0.25, 0.5, 0.75), p)
  do.call(rbind, lapply(1:3, function(i) data.frame(
    period = c("preparation", "competition")[p], share = c(0.25, 0.5, 0.75)[i],
    t(summarize(risk[i, ]))
  )))
}))
rownames(selected_rows) <- NULL

# 75%-minus-25% contrast, fixed profile.
share_contrast_rows <- do.call(rbind, lapply(1:2, function(p) {
  risk <- profile_risk_draws(c(0.25, 0.75), p)
  data.frame(period = c("preparation", "competition")[p], t(summarize_contrast(risk[2, ] - risk[1, ])))
}))
rownames(share_contrast_rows) <- NULL

# Competition-minus-preparation, at each of 25/50/75% share.
phase_rows <- do.call(rbind, lapply(c(0.25, 0.5, 0.75), function(value) {
  difference <- profile_risk_draws(value, 2) - profile_risk_draws(value, 1)
  data.frame(share = value, t(summarize_contrast(as.numeric(difference))))
}))
rownames(phase_rows) <- NULL

out_dir <- file.path(run_root, "05_profile_predictions", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(adoption_rows, file.path(out_dir, "adoption_risk.csv"), row.names = FALSE)
write.csv(selected_rows, file.path(out_dir, "share_risk_selected.csv"), row.names = FALSE)
write.csv(share_contrast_rows, file.path(out_dir, "share_risk_contrast.csv"), row.names = FALSE)
write.csv(phase_rows, file.path(out_dir, "period_risk_contrast.csv"), row.names = FALSE)

cat("=== adoption (fixed profile) ===\n"); print(adoption_rows, digits = 4)
cat("\n=== selected share risk (25/50/75, fixed profile) ===\n"); print(selected_rows, digits = 4)
cat("\n=== 75-25 share contrast (fixed profile) ===\n"); print(share_contrast_rows, digits = 4)
cat("\n=== competition-minus-preparation, by share ===\n"); print(phase_rows, digits = 4)
