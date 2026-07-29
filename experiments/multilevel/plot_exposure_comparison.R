library(ggplot2)

root <- "experiments/multilevel"
results <- read.csv(file.path(
  root, "runs", "06_exposure_comparison", "final",
  "standardized_risk_differences.csv"
))
stopifnot(nrow(results) == 4L, all(results$lower_95 <= results$median),
          all(results$median <= results$upper_95))

results$period <- factor(
  results$period,
  levels = c("preparation", "competition"),
  labels = c("Preparation", "Competition")
)
results$exposure <- factor(
  results$exposure,
  levels = c("share", "carbon_frequency"),
  labels = c("Carbon share", "Absolute carbon-frequency score")
)
results[c("median", "lower_95", "upper_95")] <-
  100 * results[c("median", "lower_95", "upper_95")]

plot <- ggplot(
  results,
  aes(x = period, y = median, colour = exposure, group = exposure)
) +
  geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.4, linetype = 2) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    position = position_dodge(width = 0.45), width = 0.08, linewidth = 0.7
  ) +
  geom_point(position = position_dodge(width = 0.45), size = 2.6) +
  scale_colour_manual(
    values = c("#2C7FB8", "#F28E2B"),
    labels = c(
      "Carbon share (+24.1 percentage points)",
      "Absolute score (+3.36 units)"
    ),
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_y_continuous(
    name = "Risk difference (pp)",
    breaks = seq(-5, 15, 5)
  ) +
  labs(x = NULL, colour = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(0, 0, 3, 0),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(5.5, 12, 5.5, 20)
  )

paths <- c(
  file.path(root, "plots", "exposure_scale_comparison.png"),
  file.path("paper", "figures", "fig_exposure_scale_comparison.png")
)
for (path in paths) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = 6.7, height = 3.8, dpi = 320, bg = "white")
}
