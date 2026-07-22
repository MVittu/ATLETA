# Face-validity analysis.
# Rates survey-item clarity with FVI; valid items have at least 80% acceptable
# responses ("Chiara" or "Abbastanza chiara").

data_path <- "data/raw/face_validity/Face_validity_SURVEY_ATLETICA.csv"
output_dir <- "data/processed/face_validity"
plot_dir <- "plots/face_validity"

if (!file.exists(data_path)) {
  stop("Face validity dataset not found in the data folder.")
}

face_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = TRUE)

if (ncol(face_data) < 3) {
  stop("The dataset must include timestamp, item responses, and comments columns.")
}

item_cols <- names(face_data)[2:(ncol(face_data) - 1)]
items_data <- face_data[, item_cols, drop = FALSE]

item_ids <- sprintf("I%d", seq_along(item_cols))
item_legend <- data.frame(
  Item_ID = item_ids,
  Original_Column_Name = item_cols,
  stringsAsFactors = FALSE
)

colnames(items_data) <- item_ids

convert_to_numeric <- function(x) {
  x_clean <- trimws(tolower(as.character(x)))

  ifelse(is.na(x_clean) | x_clean == "" | x_clean == "na", NA,
    ifelse(x_clean == "chiara", 3,
      ifelse(x_clean == "abbastanza chiara", 2,
        ifelse(x_clean == "poco chiara", 1,
          ifelse(x_clean %in% c("non chiara", "per niente chiara"), 0, NA)
        )
      )
    )
  )
}

items_data_numeric <- as.data.frame(lapply(items_data, convert_to_numeric))

acceptable_values <- c(2, 3)
clear_value <- 3
validity_cutoff <- 0.80

results <- data.frame(
  Item_ID = character(),
  Total_Respondents_N = integer(),
  Clear_Count = integer(),
  Acceptable_Count = integer(),
  Mean_Clarity_Score = numeric(),
  Item_FVI_Clear = numeric(),
  Item_FVI_Acceptable = numeric(),
  Is_Valid = logical(),
  stringsAsFactors = FALSE
)

for (col_name in colnames(items_data_numeric)) {
  responses <- na.omit(items_data_numeric[[col_name]])
  N <- length(responses)

  if (N == 0) next

  clear_count <- sum(responses == clear_value)
  acceptable_count <- sum(responses %in% acceptable_values)
  item_fvi_clear <- clear_count / N
  item_fvi_acceptable <- acceptable_count / N
  is_valid <- item_fvi_acceptable >= validity_cutoff

  results <- rbind(results, data.frame(
    Item_ID = col_name,
    Total_Respondents_N = N,
    Clear_Count = clear_count,
    Acceptable_Count = acceptable_count,
    Mean_Clarity_Score = round(mean(responses), 3),
    Item_FVI_Clear = round(item_fvi_clear, 3),
    Item_FVI_Acceptable = round(item_fvi_acceptable, 3),
    Is_Valid = is_valid,
    stringsAsFactors = FALSE
  ))
}

print(results)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_path <- file.path(output_dir, "Face_Validity_Results.csv")
legend_path <- file.path(output_dir, "Face_Validity_Item_Legend.csv")
items_to_revise_path <- file.path(output_dir, "Face_Validity_Items_To_Revise.csv")

write.csv(results, output_path, row.names = FALSE)
write.csv(item_legend, legend_path, row.names = FALSE)

items_to_revise <- results[!results$Is_Valid, ]
write.csv(items_to_revise, items_to_revise_path, row.names = FALSE)

cat("\nAnalysis complete, results saved to:", output_path, "\n")
cat("Item legend saved to:", legend_path, "\n")
cat("Items to revise saved to:", items_to_revise_path, "\n")

library(ggplot2)
library(dplyr)

plot_data <- results %>%
  mutate(
    Validity_Label = ifelse(Is_Valid, "Valid", "Revise")
  ) %>%
  arrange(desc(Item_FVI_Acceptable))

p_face_main <- ggplot(plot_data, aes(x = reorder(Item_ID, Item_FVI_Acceptable),
                                     y = Item_FVI_Acceptable,
                                     fill = Validity_Label)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = validity_cutoff, linetype = "dashed",
             color = "#2c3e50", size = 0.5) +
  geom_text(aes(label = round(Item_FVI_Acceptable, 2)),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Valid" = "#2ecc71", "Revise" = "#e74c3c"),
                    name = "Item Status") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Face Validity Index (FVI) - Clarity Analysis",
       subtitle = "Valid items have at least 80% Chiara or Abbastanza chiara responses",
       x = "Survey Items",
       y = "Item-level FVI") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "top")

print(p_face_main)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%m%d_%H%M")

ggsave(filename = file.path(plot_dir, sprintf("Face_validity_analysis_%s.png", timestamp)),
       plot = p_face_main, width = 12, height = 6, dpi = 300, bg = "white")

cat("Plot saved successfully with timestamp:", timestamp, "\n")

if (nrow(items_to_revise) > 0) {
  cat("\n")
  cat("========================================\n")
  cat("ITEMS TO REVIEW (Failed Face Validity):\n")
  cat("========================================\n")
  for (i in seq_len(nrow(items_to_revise))) {
    cat(sprintf("%s: FVI = %.3f, Mean clarity = %.3f\n",
                items_to_revise$Item_ID[i],
                items_to_revise$Item_FVI_Acceptable[i],
                items_to_revise$Mean_Clarity_Score[i]))
  }
  cat(sprintf("\nTotal items to revise: %d out of %d\n",
              nrow(items_to_revise), nrow(results)))
} else {
  cat("All items passed face-validity validation!\n")
}

cat("\nScale: Chiara = 3, Abbastanza chiara = 2, Poco chiara = 1, Non chiara = 0.\n")
cat("Face-validity cutoff:", validity_cutoff, "\n")
