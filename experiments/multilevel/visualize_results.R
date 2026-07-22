library(ggplot2)

root <- "experiments/multilevel"
plot_dir <- file.path(root, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

primary <- read.csv(file.path(
  root, "runs/01_primary/association/tables/association_estimands.csv"
))
sensitivity <- read.csv(file.path(
  root, "runs/02_sensitivity/tables/all_association_estimands.csv"
))
stopifnot(
  all(c("period", "estimand", "mean", "lower_95", "upper_95") %in% names(primary)),
  all(c("model", "probability_above_null") %in% names(sensitivity))
)

period_colors <- c(Preparation = "#0072B2", Competition = "#D55E00")
theme_results <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    plot.title.position = "plot", plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey30"), legend.position = "top"
  )

# 1. Main inferential figure: all primary adjusted associations.
effects <- primary[primary$measure == "prevalence_difference", ]
effects$Period <- factor(tools::toTitleCase(effects$period), c("Preparation", "Competition"))
effects$Comparison <- ifelse(
  effects$estimand == "adoption_difference", "Adoption", "Higher share among users"
)
effects$label <- paste(effects$Period, effects$Comparison, sep = " - ")
effects$label <- factor(effects$label, rev(c(
  "Preparation - Adoption", "Preparation - Higher share among users",
  "Competition - Adoption", "Competition - Higher share among users"
)))

p_effects <- ggplot(effects, aes(mean, label, color = Period)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  geom_vline(xintercept = 0, color = "grey35", linetype = 2) +
  geom_segment(aes(x = lower_95, xend = upper_95, yend = label), linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Carbon use and adjusted injury prevalence",
    subtitle = "Adoption: median user vs nonuser; higher share: Q75 vs Q25 among users",
    x = "Adjusted prevalence difference", y = NULL,
    caption = "Points are posterior means; bars are 95% credible intervals. Grey band: +/-2 percentage points."
  ) +
  theme_results

ggsave(file.path(plot_dir, "primary_associations.png"), p_effects,
       width = 9.5, height = 5.8, dpi = 300, bg = "white")

# 2. Robustness figure: the key competition-dose result across specifications.
robust <- sensitivity[
  sensitivity$period == "competition" & sensitivity$estimand == "dose_difference",
]
model_labels <- c(
  primary_association = "Primary nonlinear share",
  sensitivity_linear = "Linear share",
  sensitivity_wide_prior = "Wider priors"
)
robust$Model <- factor(unname(model_labels[robust$model]), rev(unname(model_labels)))
robust$probability <- sprintf("P(positive) %.1f%%", 100 * robust$probability_above_null)

p_sensitivity <- ggplot(robust, aes(mean, Model, color = Model)) +
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.2) +
  geom_vline(xintercept = 0, color = "grey35", linetype = 2) +
  geom_segment(aes(x = lower_95, xend = upper_95, yend = Model), linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(x = 0.325, label = probability), color = "grey25", hjust = 0, size = 3.8) +
  scale_color_manual(values = c(
    "Primary nonlinear share" = "#0072B2", "Linear share" = "#E69F00",
    "Wider priors" = "#6A3D9A"
  ), guide = "none") +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1), limits = c(-0.04, 0.43)) +
  labs(
    title = "Competition association is directionally consistent, not size-stable",
    subtitle = "Adjusted injury-prevalence difference: Q75 vs Q25 carbon share among users",
    x = "Adjusted prevalence difference", y = NULL,
    caption = "Points are posterior means; bars are 95% credible intervals."
  ) +
  theme_results

ggsave(file.path(plot_dir, "competition_dose_sensitivity.png"), p_sensitivity,
       width = 9.5, height = 5.2, dpi = 300, bg = "white")

# 3. Absolute-scale figure: adjusted prevalence under each exposure scenario.
prevalence <- primary[primary$measure == "prevalence", ]
scenario_labels <- c(
  prevalence_nonuser = "No carbon",
  prevalence_low_user = "User - Q25 share",
  prevalence_median_user = "User - median share",
  prevalence_high_user = "User - Q75 share"
)
prevalence$Scenario <- factor(
  unname(scenario_labels[prevalence$estimand]), unname(scenario_labels)
)
prevalence$Period <- factor(
  tools::toTitleCase(prevalence$period), c("Preparation", "Competition")
)

p_prevalence <- ggplot(prevalence, aes(Scenario, mean, color = Period)) +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0.12, linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~Period) +
  scale_color_manual(values = period_colors, guide = "none") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0.18, 0.60)) +
  labs(
    title = "Adjusted injury prevalence by carbon-use scenario",
    subtitle = "Standardized over the observed athlete and training covariates",
    x = NULL, y = "Adjusted injury prevalence",
    caption = "Intervals are 95% credible intervals; user quartiles are period-specific."
  ) +
  theme_results +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(plot_dir, "adjusted_prevalence.png"), p_prevalence,
       width = 10, height = 5.8, dpi = 300, bg = "white")

stopifnot(all(file.info(list.files(plot_dir, full.names = TRUE))$size > 1000))
cat("Saved 3 result figures to", plot_dir, "\n")
