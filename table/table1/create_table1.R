options(stringsAsFactors = FALSE)

source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
population_file <- "tables/recall/population_history_by_age_sex.csv"
output_file <- "table/table1/table1_by_sex.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w)
  x <- x[keep]
  w <- w[keep]
  order_x <- order(x)
  cumulative <- cumsum(w[order_x]) / sum(w)
  vapply(probs, function(prob) x[order_x][which(cumulative >= prob)[1]], numeric(1))
}
weighted_var <- function(x, w) {
  keep <- is.finite(x) & is.finite(w)
  x <- x[keep]
  w <- w[keep]
  mean_x <- weighted.mean(x, w)
  sum(w * (x - mean_x)^2) / sum(w)
}
format_cell <- function(x, w, type) {
  keep <- is.finite(x) & is.finite(w)
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) return(NA_character_)
  if (type == "binary") {
    count <- sum(w * x)
    return(sprintf("%.1f (%.1f%%)", count, 100 * count / sum(w)))
  }
  if (type == "median") {
    q <- weighted_quantile(x, w, c(0.25, 0.5, 0.75))
    return(sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3]))
  }
  sprintf("%.1f (%.1f)", weighted.mean(x, w), sqrt(weighted_var(x, w)))
}
smd <- function(x, w, sex) {
  keep <- is.finite(x) & is.finite(w)
  x <- x[keep]
  w <- w[keep]
  sex <- sex[keep]
  women <- sex == "Women"
  men <- sex == "Men"
  denominator <- sqrt((weighted_var(x[women], w[women]) + weighted_var(x[men], w[men])) / 2)
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  (weighted.mean(x[women], w[women]) - weighted.mean(x[men], w[men])) / denominator
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
sex <- factor(ifelse(raw$genere_gara == 1, "Men", "Women"), c("Men", "Women"))
age_group <- cut(as.numeric(raw$eta), c(-Inf, 19, 22, Inf),
                 labels = c("Under 20", "Under 23", "Senior"))
population_history <- read.csv(population_file)
target <- do.call(rbind, lapply(split(
  population_history,
  list(population_history$age_group, population_history$sex), drop = TRUE
), function(cell) {
  current <- cell$n[cell$year == 2026]
  estimate <- if (length(current)) current else round(mean(cell$n[cell$year < 2026]))
  data.frame(age_group = cell$age_group[1], sex = cell$sex[1], target_n = estimate)
}))
target$age_group <- c(under20 = "Under 20", under23 = "Under 23", assoluta = "Senior")[target$age_group]
target$sex <- c(uomini = "Men", donne = "Women")[target$sex]
sample_cells <- as.data.frame(table(age_group = age_group, sex = sex), responseName = "sample_n")
weight_cells <- merge(target, sample_cells, by = c("age_group", "sex"))
weight_cells$weight <- weight_cells$target_n / weight_cells$sample_n
weights <- weight_cells$weight[match(paste(age_group, sex),
                                     paste(weight_cells$age_group, weight_cells$sex))]
stopifnot(J > 0L, sum(weight_cells$target_n) == 1328L, all(is.finite(weights)))

height <- as.numeric(raw$statura_cm)
height[height > 0 & height < 3] <- height[height > 0 & height < 3] * 100
prep_total <- prep_carbon + prep_noncarbon
competition_total <- competition_carbon + competition_noncarbon
prep_share <- prep_carbon / prep_total
competition_share <- competition_carbon / competition_total

rows <- list()
add_row <- function(section, characteristic, x, type = "mean", include = rep(TRUE, J)) {
  x[!include] <- NA_real_
  men <- sex == "Men"
  women <- sex == "Women"
  rows[[length(rows) + 1L]] <<- data.frame(
    section = section, characteristic = characteristic,
    display = c(mean = "Mean (SD)", median = "Median [Q1, Q3]", binary = "n (%)")[type],
    overall_sample = format_cell(x, rep(1, J), type),
    men_sample = format_cell(x[men], rep(1, sum(men)), type),
    women_sample = format_cell(x[women], rep(1, sum(women)), type),
    smd_sample_women_vs_men = smd(x, rep(1, J), sex),
    overall_target_weighted = format_cell(x, weights, type),
    men_target_weighted = format_cell(x[men], weights[men], type),
    women_target_weighted = format_cell(x[women], weights[women], type),
    smd_weighted_women_vs_men = smd(x, weights, sex),
    missing_n = sum(is.na(x) & include), stringsAsFactors = FALSE
  )
}

rows[[1]] <- data.frame(
  section = "Population", characteristic = "Athletes",
  display = "N", overall_sample = as.character(J),
  men_sample = as.character(sum(sex == "Men")), women_sample = as.character(sum(sex == "Women")),
  smd_sample_women_vs_men = NA_real_, overall_target_weighted = "1328",
  men_target_weighted = "743", women_target_weighted = "585",
  smd_weighted_women_vs_men = NA_real_, missing_n = 0L
)
add_row("Demographics", "Age, years", as.numeric(raw$eta))
for (level in levels(age_group))
  add_row("Demographics", paste0("Age group: ", level), as.numeric(age_group == level), "binary")
add_row("Anthropometrics", "Height, cm", height)
add_row("Anthropometrics", "Weight, kg", as.numeric(raw$peso_kg))
add_row("Athletics", "Athletics experience, years", as.numeric(raw$anni_atletica))
add_row("Injury", "Prior injury in previous two seasons", as.numeric(raw$infortunio_ultime_due_stagioni), "binary")
add_row("Injury", "Any current-season injury", as.numeric(raw$infortunio_stagione > 0), "binary")
add_row("Injury", "Preparation-period injury", prep_injury, "binary")
add_row("Injury", "Competition-period injury", competition_injury, "binary")
add_row("Carbon exposure", "Any preparation carbon use", as.numeric(prep_carbon > 0), "binary")
add_row("Carbon exposure", "Any competition carbon use", as.numeric(competition_carbon > 0), "binary")
add_row("Carbon exposure", "Preparation carbon share among users", prep_share, "median", prep_carbon > 0)
add_row("Carbon exposure", "Competition carbon share among users", competition_share, "median", competition_carbon > 0)
add_row("Training frequency", "Preparation total type-frequency", prep_total, "median")
add_row("Training frequency", "Competition total type-frequency", competition_total, "median")

domain_labels <- c(
  fondo_lungo = "Long endurance", fondo_medio = "Middle endurance",
  lattacido = "Lactate", max_velocity = "Maximum velocity", tecnica = "Technique"
)
for (domain in domains) {
  add_row("Preparation frequency", domain_labels[domain],
          raw[[paste0("prep_carbon_", domain)]] + raw[[paste0("prep_no_carbon_", domain)]],
          "median")
}
add_row("Preparation frequency", "Gym", as.numeric(raw$prep_no_carbon_palestra), "median")
for (domain in domains) {
  add_row("Competition frequency", domain_labels[domain],
          raw[[paste0("gara_carbon_", domain)]] + raw[[paste0("gara_no_carbon_", domain)]],
          "median")
}
add_row("Competition frequency", "Gym", as.numeric(raw$gara_no_carbon_palestra), "median")

discipline_labels <- c(
  "Speed", "Middle distance", "Long distance", "Hurdles",
  "Horizontal jumps", "Pole vault", "Combined events"
)
for (i in seq_along(discipline_labels))
  add_row("Discipline", discipline_labels[i], as.numeric(raw[[paste0("disciplina___", i)]]), "binary")

table1 <- do.call(rbind, rows)
table1$smd_sample_women_vs_men <- round(table1$smd_sample_women_vs_men, 3)
table1$smd_weighted_women_vs_men <- round(table1$smd_weighted_women_vs_men, 3)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(table1, output_file, row.names = FALSE)
stopifnot(nrow(table1) == 37L, file.info(output_file)$size > 1000)
cat("Saved", output_file, "\n")
