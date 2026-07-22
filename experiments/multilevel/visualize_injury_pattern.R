options(stringsAsFactors = FALSE)
library(ggplot2)

root <- "experiments/multilevel"
plot_dir <- file.path(root, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
source_file <- "data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv"
domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")

raw_all <- read.csv(source_file, check.names = FALSE)
eligible <- raw_all$consenso == 1 & raw_all$eleggibile_calc == 1 &
  raw_all$survey_atleta_it_complete == 2
raw <- raw_all[eligible, , drop = FALSE]
sum_named <- function(prefix) rowSums(raw[paste0(prefix, domains)])
prep_carbon <- sum_named("prep_carbon_"); prep_noncarbon <- sum_named("prep_no_carbon_")
competition_carbon <- sum_named("gara_carbon_"); competition_noncarbon <- sum_named("gara_no_carbon_")
count <- raw$infortunio_stagione
complete_timing <- count == 0 |
  (count == 1 & !is.na(raw$inj_sing_periodo)) |
  (count == 2 & !is.na(raw$inj_multi1_periodo) & !is.na(raw$inj_multi2_periodo))
valid <- complete_timing & prep_carbon + prep_noncarbon > 0 &
  competition_carbon + competition_noncarbon > 0
raw <- raw[valid, , drop = FALSE]
stopifnot(nrow(raw) > 0L)

loc_labels <- c("Hip", "Groin", "Thigh", "Lower leg", "Ankle", "Foot", "Head/upper limb/trunk")
typ_labels <- c("Joint sprain", "Tendon rupture", "Ligament injury", "Muscle strain",
                "Bone stress injury\nor stress fracture", "Tendinopathy",
                "Plantar fasciopathy", "Other")

extract_injuries <- function(raw, loc_prefixes, typ_prefixes, which_records) {
  out <- list()
  for (i in seq_len(nrow(raw))) {
    row <- raw[i, ]
    c <- row$infortunio_stagione
    if (c == 1 && "single" %in% which_records) {
      out[[length(out) + 1]] <- list(
        loc = as.numeric(row[paste0("inj_sing_localizzazione___", 1:7)]),
        typ = as.numeric(row[paste0("inj_sing_tipologia___", 1:8)])
      )
    } else if (c == 2) {
      out[[length(out) + 1]] <- list(
        loc = as.numeric(row[paste0("inj_multi1_localizzazione___", 1:7)]),
        typ = as.numeric(row[paste0("inj_multi1_tipologia___", 1:8)])
      )
      out[[length(out) + 1]] <- list(
        loc = as.numeric(row[paste0("inj_multi2_localizzazione___", 1:7)]),
        typ = as.numeric(row[paste0("inj_multi2_tipologia___", 1:8)])
      )
    }
  }
  out
}
injuries <- extract_injuries(raw, which_records = "single")
n_injuries <- length(injuries)
loc_mat <- do.call(rbind, lapply(injuries, `[[`, "loc"))
typ_mat <- do.call(rbind, lapply(injuries, `[[`, "typ"))
loc_mat[is.na(loc_mat)] <- 0
typ_mat[is.na(typ_mat)] <- 0

is_upper <- loc_mat[, 7] == 1
n_lower <- sum(!is_upper)

loc_counts <- colSums(loc_mat)
typ_counts <- colSums(typ_mat[!is_upper, , drop = FALSE])

location_panel <- sprintf("Location (%% of %d injuries)", n_injuries)
type_panel <- sprintf("Type (%% of %d lower-limb injuries)", n_lower)
loc_df <- data.frame(
  Panel = location_panel, Category = loc_labels,
  Percent = 100 * loc_counts / n_injuries
)
typ_df <- data.frame(
  Panel = type_panel, Category = typ_labels,
  Percent = 100 * typ_counts / n_lower
)
plot_data <- rbind(loc_df, typ_df)
plot_data$Panel <- factor(plot_data$Panel, c(location_panel, type_panel))
plot_data$Category <- factor(plot_data$Category, levels = plot_data$Category[order(plot_data$Percent)])

p <- ggplot(plot_data, aes(Percent, Category, fill = Panel)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), hjust = -0.12, size = 3.4, color = "grey20") +
  facet_wrap(~Panel, scales = "free_y") +
  scale_fill_manual(values = c(
    setNames(c("#0072B2", "#D55E00"), c(location_panel, type_panel))
  )) +
  scale_x_continuous(limits = c(0, max(plot_data$Percent) * 1.22), expand = c(0, 0)) +
  labs(
    title = "Injury location and type in the balanced cohort",
    subtitle = "Multi-select items; percentages sum to more than 100% within each panel",
    x = "Percent", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"), plot.title.position = "plot",
    plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey30")
  )

ggsave(file.path(plot_dir, "injury_location_type.png"), p,
       width = 10.5, height = 5, dpi = 300, bg = "white")
cat("Saved injury_location_type.png\n")
cat("n_injuries", n_injuries, "n_lower", n_lower, "\n")
