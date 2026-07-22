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
prep_carbon <- sum_named("prep_carbon_")
prep_noncarbon <- sum_named("prep_no_carbon_")
competition_carbon <- sum_named("gara_carbon_")
competition_noncarbon <- sum_named("gara_no_carbon_")

count <- raw$infortunio_stagione
complete_timing <- count == 0 |
  (count == 1 & !is.na(raw$inj_sing_periodo)) |
  (count == 2 & !is.na(raw$inj_multi1_periodo) & !is.na(raw$inj_multi2_periodo))
valid <- complete_timing & prep_carbon + prep_noncarbon > 0 &
  competition_carbon + competition_noncarbon > 0
raw <- raw[valid, , drop = FALSE]
prep_carbon <- prep_carbon[valid]; prep_noncarbon <- prep_noncarbon[valid]
competition_carbon <- competition_carbon[valid]; competition_noncarbon <- competition_noncarbon[valid]

stopifnot(nrow(raw) > 0L)

sex <- factor(ifelse(raw$genere_gara == 1, "Men", "Women"), c("Men", "Women"))
age_group <- cut(as.numeric(raw$eta), c(-Inf, 19, 22, Inf),
                 labels = c("Under 20", "Under 23", "Senior"))

plot_data <- data.frame(
  Sex = rep(sex, 2), AgeGroup = rep(age_group, 2),
  Period = factor(rep(c("Preparation", "Competition"), each = nrow(raw)),
                  c("Preparation", "Competition")),
  Share = c(prep_carbon / (prep_carbon + prep_noncarbon),
           competition_carbon / (competition_carbon + competition_noncarbon))
)

period_colors <- c(Men = "#0072B2", Women = "#D55E00")
p <- ggplot(plot_data, aes(Share, fill = Sex, color = Sex)) +
  geom_density(alpha = 0.35, linewidth = 0.7, bounds = c(0, 1)) +
  facet_grid(rows = vars(Period), cols = vars(AgeGroup)) +
  scale_fill_manual(values = period_colors) +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Carbon-plate share of training is stable across sex and age group",
    subtitle = "Observed carbon share of reported training-type frequency, by period, age group, and sex",
    x = "Carbon share of training-type frequency", y = "Density",
    caption = sprintf(
      "Density estimated on the balanced analysis cohort (N = %d); Under 20 n = %d, Under 23 n = %d, Senior n = %d.",
      nrow(raw), sum(age_group == "Under 20"), sum(age_group == "Under 23"),
      sum(age_group == "Senior")
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), legend.position = "top",
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"), strip.text = element_text(face = "bold")
  )

ggsave(file.path(plot_dir, "descriptive_carbon_share_by_stratum.png"), p,
       width = 10.5, height = 6, dpi = 300, bg = "white")
cat("Saved descriptive_carbon_share_by_stratum.png\n")
