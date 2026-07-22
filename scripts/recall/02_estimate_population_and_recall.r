# Estimate target population and recall counts by stratum.

library(readxl)

population_path <- "data/raw/recall/raccolta_dati_campione.xlsx"
membership_path <- "data/processed/recall/survey_strata_membership.csv"
private_out_dir <- "data/processed/recall"
public_out_dir <- "tables/recall"

prior_a <- 1
prior_b <- 1
credible_level <- 0.95
target_half_width <- 0.06

read_num <- function(x) suppressWarnings(as.numeric(x))

bayesian_target_n <- function(population_n) {
  alpha <- (1 - credible_level) / 2
  for (n in seq_len(population_n)) {
    successes <- floor(n / 2)
    lower <- qbeta(alpha, prior_a + successes, prior_b + n - successes)
    upper <- qbeta(1 - alpha, prior_a + successes, prior_b + n - successes)
    if ((upper - lower) / 2 <= target_half_width) return(n)
  }
  population_n
}

summary_sheet <- as.data.frame(suppressMessages(read_excel(population_path, sheet = "prospetto finale",
                                                            col_names = FALSE)))
membership <- read.csv(membership_path, stringsAsFactors = FALSE)

category_rows <- list(under20 = 3, under23 = 4, assoluta = 5)
years <- list(`2026` = c(2, 3), `2025` = c(5, 6), `2024` = c(8, 9), `2023` = c(11, 12))
sexes <- c("uomini", "donne")

sex_totals <- do.call(rbind, lapply(names(category_rows), function(age_group) {
  row <- category_rows[[age_group]]
  do.call(rbind, lapply(names(years), function(year) {
    data.frame(
      age_group = age_group,
      year = as.integer(year),
      sex = sexes,
      n = read_num(unlist(summary_sheet[row, years[[year]]])),
      stringsAsFactors = FALSE
    )
  }))
}))
sex_totals <- sex_totals[!is.na(sex_totals$n), ]

estimate_total <- do.call(rbind, lapply(split(sex_totals, list(sex_totals$age_group, sex_totals$sex), drop = TRUE), function(x) {
  current <- x$n[x$year == 2026]
  history <- x$n[x$year < 2026]
  estimate <- if (length(current) && !is.na(current)) current else round(mean(history))
  data.frame(
    age_group = x$age_group[1],
    sex = x$sex[1],
    estimate_n = estimate,
    lower_n = min(c(estimate, history), na.rm = TRUE),
    upper_n = max(c(estimate, history), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

section_rows <- list(under20 = 12:18, under23 = 25:31, assoluta = 38:44)
specialty_shares <- do.call(rbind, lapply(names(section_rows), function(age_group) {
  rows <- section_rows[[age_group]]
  data.frame(
    age_group = age_group,
    specialty = summary_sheet[[1]][rows],
    sex = rep(sexes, each = length(rows)),
    n_2025 = c(read_num(summary_sheet[[5]][rows]), read_num(summary_sheet[[7]][rows])),
    stringsAsFactors = FALSE
  )
}))
specialty_shares <- specialty_shares[specialty_shares$specialty != "totale", ]
specialty_shares$share <- ave(specialty_shares$n_2025, specialty_shares$age_group, specialty_shares$sex,
                              FUN = function(x) x / sum(x, na.rm = TRUE))

population <- merge(specialty_shares, estimate_total, by = c("age_group", "sex"))
population$population_estimate <- round(population$estimate_n * population$share)
population$population_lower <- round(population$lower_n * population$share)
population$population_upper <- round(population$upper_n * population$share)

total_population <- sum(population$population_estimate)
target_total <- bayesian_target_n(total_population)

population$target_raw <- target_total * population$population_estimate / total_population
population$target_n <- floor(population$target_raw)
missing_targets <- target_total - sum(population$target_n)
if (missing_targets > 0) {
  add_to <- order(population$target_raw - population$target_n, decreasing = TRUE)[seq_len(missing_targets)]
  population$target_n[add_to] <- population$target_n[add_to] + 1
}

current <- aggregate(weight ~ age_group + sex + specialty, membership, sum)
names(current)[names(current) == "weight"] <- "current_n"

recall <- merge(population, current, by = c("age_group", "sex", "specialty"), all.x = TRUE)
recall$current_n[is.na(recall$current_n)] <- 0
recall$recall_n <- ceiling(pmax(0, recall$target_n - recall$current_n))
recall <- recall[order(recall$age_group, recall$sex, -recall$recall_n), ]

dir.create(private_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(public_out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(sex_totals, file.path(public_out_dir, "population_history_by_age_sex.csv"), row.names = FALSE)
write.csv(population[, c("age_group", "sex", "specialty", "population_estimate",
                         "population_lower", "population_upper")],
          file.path(public_out_dir, "population_2026_estimate_by_stratum.csv"), row.names = FALSE)
write.csv(current, file.path(private_out_dir, "current_sample_by_stratum.csv"), row.names = FALSE)
write.csv(recall[, c("age_group", "sex", "specialty", "population_estimate",
                     "population_lower", "population_upper", "current_n", "target_n", "recall_n")],
          file.path(private_out_dir, "recall_targets_by_stratum.csv"), row.names = FALSE)
write.csv(data.frame(
  method = "Bayesian beta-binomial worst-case proportion",
  prior = sprintf("Beta(%s,%s)", prior_a, prior_b),
  credible_level = credible_level,
  target_half_width = target_half_width,
  population_estimate = total_population,
  target_sample = target_total
), file.path(public_out_dir, "sample_size_assumptions.csv"), row.names = FALSE)

cat("Estimated population:", total_population, "\n")
cat("Target sample:", target_total, "\n")
cat("Additional recalls:", sum(recall$recall_n), "\n")
