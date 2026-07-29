source("experiments/multilevel/visualize_descriptive_share.R")

sex_colors <- c(Men = "#16A6B6", Women = "#D84A73")
p <- p +
  scale_fill_manual(values = sex_colors) +
  scale_color_manual(values = sex_colors)

paths <- c(
  file.path(plot_dir, "descriptive_carbon_share_by_stratum.png"),
  file.path("paper", "figures", "fig_descriptive_carbon_share_by_stratum.png")
)
for (path in paths) {
  ggsave(path, p, width = 10.5, height = 6, dpi = 320, bg = "white")
}
