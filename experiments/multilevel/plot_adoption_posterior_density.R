library(ggplot2)

root <- "experiments/multilevel"
draws <- readRDS(file.path(
  root, "runs", "01_primary", "association", "models", "posterior_draws.rds"
))
profile <- read.csv(file.path(
  root, "runs", "05_profile_predictions", "tables", "reference_profile.csv"
))
reference_x <- profile$value[profile$type == "covariate"]
reference_discipline <- profile$value[profile$type == "discipline_indicator"]
stopifnot(
  length(reference_x) == ncol(draws$beta_common),
  length(reference_discipline) == ncol(draws$beta_discipline)
)

adoption <- do.call(rbind, lapply(1:2, function(period) {
  base <- draws$alpha[, period] +
    as.numeric(draws$beta_common %*% reference_x) +
    as.numeric(draws$beta_discipline %*% reference_discipline)
  data.frame(
    Period = factor(c("Preparation", "Competition")[period],
                    c("Preparation", "Competition")),
    difference = plogis(base + draws$beta_any[, period]) - plogis(base)
  )
}))
summary <- do.call(rbind, lapply(split(adoption$difference, adoption$Period), function(x) {
  data.frame(
    median = median(x),
    lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975))
  )
}))
summary$Period <- factor(rownames(summary), c("Preparation", "Competition"))
colors <- c(Preparation = "#0072B2", Competition = "#D55E00")

plot <- ggplot(adoption, aes(difference, fill = Period)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.22) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.45, linetype = 2) +
  geom_density(alpha = 0.82, colour = "white", linewidth = 0.35) +
  geom_segment(
    data = summary,
    aes(x = lower, xend = upper, y = 0, yend = 0),
    inherit.aes = FALSE, linewidth = 0.9
  ) +
  geom_point(
    data = summary, aes(x = median, y = 0),
    inherit.aes = FALSE, shape = 21, fill = "white", size = 2.6, stroke = 0.75
  ) +
  facet_grid(Period ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(NULL, breaks = NULL, expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Posterior distribution of the adoption contrast",
    subtitle = "User at exactly 50% carbon share minus non-user, for the same fixed profile",
    x = "Difference in predicted injury risk", y = NULL,
    caption = "Density: posterior distribution; dot and segment: median and 95% credible interval; grey band: ±2-pp ROPE."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "#444444"),
    plot.caption = element_text(hjust = 0),
    panel.spacing.y = grid::unit(0.35, "lines")
  )

paths <- c(
  file.path(root, "plots", "posterior_association_distributions.png"),
  file.path("paper", "figures", "fig_posterior_association_distributions.png")
)
for (path in paths) {
  ggsave(path, plot, width = 8, height = 4.8, dpi = 320, bg = "white")
}
