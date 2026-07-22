# Clean first-week survey data for recall planning.

raw_path <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
out_dir <- "data/processed/recall"

survey <- read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE)

disc_cols <- paste0("disciplina___", 1:11)
disciplines <- c(
  "velocità",
  "mezzofondo",
  "fondo",
  "ostacoli",
  "salti in estensione",
  "salto con l'asta",
  "prove multiple",
  "lanci",
  "salto in alto",
  "marcia",
  "altre discipline non track and field"
)

stopifnot(length(disc_cols) == length(disciplines), all(disc_cols %in% names(survey)))

carbon_cols <- grep("^(prep|gara)_carbon_", names(survey), value = TRUE)

clean <- data.frame(
  source_row = seq_len(nrow(survey)),
  record_id = survey[["record_id"]],
  eligible = survey[["eleggibile_calc"]] == 1,
  complete = survey[["survey_atleta_it_complete"]] == 2,
  sex = ifelse(survey[["genere_gara"]] == 1, "uomini",
               ifelse(survey[["genere_gara"]] == 2, "donne", NA)),
  age = suppressWarnings(as.numeric(survey[["eta"]])),
  age_group = cut(suppressWarnings(as.numeric(survey[["eta"]])),
                  breaks = c(-Inf, 19, 22, Inf),
                  labels = c("under20", "under23", "assoluta")),
  carbon_insole_usage = rowSums(survey[carbon_cols] > 0, na.rm = TRUE) > 0,
  injury_season = suppressWarnings(as.numeric(survey[["infortunio_stagione"]])) > 0,
  stringsAsFactors = FALSE
)

clean <- clean[clean$eligible & clean$complete & !is.na(clean$sex) &
                 !is.na(clean$age_group), ]

membership <- do.call(rbind, lapply(seq_len(nrow(clean)), function(i) {
  selected <- unique(disciplines[which(as.numeric(survey[clean$source_row[i], disc_cols]) == 1)])
  selected <- selected[!is.na(selected)]
  if (!length(selected)) return(NULL)

  data.frame(
    record_id = clean$record_id[i],
    age_group = clean$age_group[i],
    sex = clean$sex[i],
    specialty = selected,
    weight = 1 / length(selected),
    stringsAsFactors = FALSE
  )
}))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(clean, file.path(out_dir, "survey_clean.csv"), row.names = FALSE)
write.csv(membership, file.path(out_dir, "survey_strata_membership.csv"), row.names = FALSE)

cat("Clean survey rows:", nrow(clean), "\n")
cat("Strata memberships:", nrow(membership), "\n")
