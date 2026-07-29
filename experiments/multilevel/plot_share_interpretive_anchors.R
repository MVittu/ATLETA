library(ggplot2)

root <- "experiments/multilevel"
profile_dir <- file.path(root, "runs", "05_profile_predictions", "tables")
curve <- read.csv(file.path(profile_dir, "share_risk_curve.csv"))
anchors <- read.csv(file.path(profile_dir, "share_risk_selected.csv"))
stopifnot(all(sort(unique(anchors$share)) == c(0.25, 0.50, 0.75)))

domains <- c("fondo_lungo", "fondo_medio", "lattacido", "max_velocity", "tecnica")
raw <- read.csv("data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv",
                check.names = FALSE)
raw <- raw[raw$consenso == 1 & raw$eleggibile_calc == 1 &
             raw$survey_atleta_it_complete == 2, , drop = FALSE]
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
rug <- data.frame(
  Period = rep(c("Preparation", "Competition"), each = sum(valid)),
  share = c(
    prep_carbon[valid] / (prep_carbon[valid] + prep_noncarbon[valid]),
    competition_carbon[valid] /
      (competition_carbon[valid] + competition_noncarbon[valid])
  )
)
rug <- rug[rug$share > 0, ]

levels <- c("Preparation", "Competition")
curve$Period <- factor(curve$Period, levels)
anchors$Period <- factor(anchors$Period, levels)
rug$Period <- factor(rug$Period, levels)
anchors$label <- sprintf("%.1f%%", 100 * anchors$median)
colors <- c(Preparation = "#0072B2", Competition = "#D55E00")

plot <- ggplot(curve, aes(share, median, colour = Period, fill = Period)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_vline(xintercept = c(0.25, 0.50, 0.75),
             colour = "grey55", linewidth = 0.4, linetype = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(
    data = anchors, aes(ymin = lower, ymax = upper),
    width = 0.012, linewidth = 0.65
  ) +
  geom_point(data = anchors, shape = 21, fill = "white", size = 2.5, stroke = 0.8) +
  geom_label(
    data = anchors, aes(label = label),
    nudge_x = -0.035, nudge_y = 0.035,
    fill = "white", linewidth = 0, label.padding = grid::unit(0.08, "lines"),
    show.legend = FALSE, size = 3.1
  ) +
  geom_rug(
    data = rug, aes(x = share, colour = Period), inherit.aes = FALSE,
    sides = "b", alpha = 0.18, length = grid::unit(0.035, "npc")
  ) +
  facet_wrap(~Period, nrow = 1) +
  scale_colour_manual(values = colors, guide = "none") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1), labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Continuous fitted risk with three interpretive anchors",
    subtitle = "25%, 50%, and 75% are exact carbon-share values—not percentiles, bins, or clinical thresholds",
    x = "Carbon share of reported training-frequency score",
    y = "Predicted injury risk",
    caption = paste(
      "Points read the same continuous natural-spline curve at three round values;",
      "bars are 95% credible intervals; rug marks observed user shares."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "#444444"),
    plot.caption = element_text(hjust = 0),
    panel.spacing = grid::unit(1.1, "lines")
  )

paths <- c(
  file.path(root, "plots", "spline_dose_response.png"),
  file.path("paper", "figures", "fig_spline_dose_response.png")
)
for (path in paths) {
  ggsave(path, plot, width = 8.5, height = 5.5, dpi = 320, bg = "white")
}
