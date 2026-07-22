options(stringsAsFactors = FALSE)

source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
output_file <- "experiments/multilevel/runs/00_data_audit/tables/ostrc_by_injury_period.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
items <- c("diff_partecipazione", "modifica_allenamenti",
           "influenza_prestazione", "sintomi_infortunio")

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
  multiple <- cbind(raw$inj_multi1_periodo, raw$inj_multi2_periodo)
  outcome <- rep(NA_integer_, nrow(raw))
  outcome[raw$infortunio_stagione == 0] <- 0L
  one <- raw$infortunio_stagione == 1
  outcome[one] <- as.integer(raw$inj_sing_periodo[one] %in% period_codes)
  hit <- rowSums((multiple == period_codes[1]) | (multiple == period_codes[2]),
                 na.rm = TRUE) > 0
  outcome[raw$infortunio_stagione == 2] <- hit[raw$infortunio_stagione == 2]
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
prep_injury <- prep_injury[valid]
competition_injury <- competition_injury[valid]

item_codes <- sapply(raw[items], function(x) suppressWarnings(as.integer(x)))
item_codes[, 1] <- item_codes[, 1] - 1L
complete_ostrc <- raw$infortunio_stagione > 0 & complete.cases(item_codes)
stopifnot(nrow(raw) > 0L, sum(complete_ostrc) > 0L,
          all(item_codes[complete_ostrc, ] %in% 0:3))

score_map <- c(0, 8, 17, 25)
scores <- matrix(score_map[item_codes[complete_ostrc, ] + 1L],
                 nrow = sum(complete_ostrc), dimnames = list(NULL, c(
                   "Participation impact", "Training-volume impact",
                   "Performance impact", "Symptoms")))
scores <- cbind(scores, "Total severity" = rowSums(scores))
period_group <- ifelse(
  prep_injury[complete_ostrc] == 1 & competition_injury[complete_ostrc] == 1,
  "Both periods",
  ifelse(prep_injury[complete_ostrc] == 1, "Preparation only", "Competition only")
)
group_order <- c("Overall", "Preparation only", "Competition only", "Both periods")
reference <- scores[period_group == "Competition only", , drop = FALSE]

result <- do.call(rbind, lapply(seq_len(ncol(scores)), function(j) {
  do.call(rbind, lapply(group_order, function(group) {
    x <- if (group == "Overall") scores[, j] else scores[period_group == group, j]
    ref <- reference[, j]
    smd <- if (group == "Overall") NA_real_ else if (group == "Competition only") 0 else
      (mean(x) - mean(ref)) / sqrt((var(x) + var(ref)) / 2)
    data.frame(
      measure = colnames(scores)[j], group = group, n = length(x),
      mean = mean(x), sd = sd(x), median = median(x),
      q1 = unname(quantile(x, 0.25)), q3 = unname(quantile(x, 0.75)),
      smd_vs_competition_only = smd
    )
  }))
}))

result[c("mean", "sd", "median", "q1", "q3", "smd_vs_competition_only")] <-
  lapply(result[c("mean", "sd", "median", "q1", "q3", "smd_vs_competition_only")], round, 2)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_file, row.names = FALSE)
stopifnot(nrow(result) == 20L, file.info(output_file)$size > 1000)
cat("Saved", output_file, "\n")
