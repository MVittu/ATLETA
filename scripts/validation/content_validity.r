# Content Validity Ratio (CVR) analysis - Bayesian approach.
# Based on Baghestani et al. (2019), Bayesian Critical Values for Lawshe's CVR.

data_path <- "data/raw/content_validity/Content validity SURVEY ATLETICA.csv"
output_dir <- "data/processed/content_validity"
plot_dir <- "plots/content_validity"

cvr_data <- read.csv(data_path, stringsAsFactors = FALSE)

cols_to_exclude <- c("X", # CSV row index, not a survey item
                     "Informazioni.cronologiche",
                     "Quale.professione.svolge.",
                     "Hai.altre.osservazioni.da.fare.")

item_cols <- setdiff(names(cvr_data), cols_to_exclude)
items_data <- cvr_data[, item_cols]

convert_to_numeric <- function(x) {
  x_clean <- trimws(tolower(as.character(x)))
  mapping <- c("non necessario" = 1, "non utile" = 1,
               "utile ma non indispensabile" = 2, "indispensabile" = 3)
  unknown <- setdiff(unique(x_clean[!is.na(x_clean) & nzchar(x_clean)]), names(mapping))
  if (length(unknown)) stop("Unknown CVR response: ", paste(unknown, collapse = ", "))
  unname(mapping[x_clean])
}

items_data_numeric <- as.data.frame(lapply(items_data, convert_to_numeric))

essential_value <- 3

get_bayesian_critical_ne <- function(N, alpha = 0.05) {
  for (ne in N:0) {
    prob_null <- pbeta(0.5, ne + 1, N - ne + 1)

    if (prob_null > alpha) {
      return(ne + 1)
    }
  }
  return(N)
}

results <- data.frame(
  Item = character(),
  Total_Experts_N = integer(),
  Essential_Count_ne = integer(),
  Observed_CVR = numeric(),
  Critical_ne_Bayesian = integer(),
  Critical_CVR_Bayesian = numeric(),
  Is_Valid = logical(),
  stringsAsFactors = FALSE
)

for (col_name in colnames(items_data_numeric)) {
  responses <- na.omit(items_data_numeric[[col_name]])
  N <- length(responses)

  if (N == 0) next

  n_e <- sum(responses == essential_value)
  observed_cvr <- (n_e - (N/2)) / (N/2)
  critical_ne <- get_bayesian_critical_ne(N)
  critical_cvr <- (critical_ne - (N/2)) / (N/2)
  is_valid <- n_e >= critical_ne

  results <- rbind(results, data.frame(
    Item = col_name,
    Total_Experts_N = N,
    Essential_Count_ne = n_e,
    Observed_CVR = round(observed_cvr, 3),
    Critical_ne_Bayesian = critical_ne,
    Critical_CVR_Bayesian = round(critical_cvr, 3),
    Is_Valid = is_valid
  ))
}

print(head(results))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, "CVR_Bayesian_Results.csv")
write.csv(results, output_path, row.names = FALSE)
cat("\nAnalysis complete, results saved to:", output_path, "\n")

library(ggplot2)
library(dplyr)

results <- results %>%
  mutate(
    Item_ID = paste0("I", row_number()),
    Item_Short = Item_ID,
    Validity_Label = ifelse(Is_Valid, "Valid", "Invalid")
  ) %>%
  arrange(desc(Observed_CVR))

p_cvr_main <- ggplot(results, aes(x = reorder(Item_Short, Observed_CVR),
                                   y = Observed_CVR,
                                   fill = Validity_Label)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  geom_text(aes(label = round(Observed_CVR, 2)),
            vjust = ifelse(results$Observed_CVR > 0, -0.5, 1.5),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Valid" = "#2ecc71", "Invalid" = "#e74c3c"),
                    name = "Item Status") +
  scale_y_continuous(limits = c(min(0, min(results$Observed_CVR) - 0.15),
                                max(results$Observed_CVR) + 0.15)) +
  labs(title = "Content Validity Ratio (CVR) - Bayesian Analysis",
       subtitle = "Green = Valid items (passed Bayesian critical value); Red = Items to revise",
       x = "Survey Items",
       y = "Observed CVR") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "top")

print(p_cvr_main)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%m%d_%H%M")

ggsave(filename = file.path(plot_dir, sprintf("CVR_analysis_%s.png", timestamp)),
       plot = p_cvr_main, width = 12, height = 6, dpi = 300, bg = "white")

cat("Plot saved successfully with timestamp:", timestamp, "\n")
cat("\nItem naming: I1, I2, I3, ... correspond to survey items for easy reference.\n")

invalid_items <- results %>% filter(!Is_Valid)

if (nrow(invalid_items) > 0) {
  cat("\n")
  cat("========================================\n")
  cat("ITEMS TO REFUTE (Failed CVR):\n")
  cat("========================================\n")
  for (i in seq_len(nrow(invalid_items))) {
    cat(sprintf("%s: %s (CVR = %.3f, Critical = %.3f)\n",
                invalid_items$Item_ID[i],
                invalid_items$Item[i],
                invalid_items$Observed_CVR[i],
                invalid_items$Critical_CVR_Bayesian[i]))
  }
  cat(sprintf("\nTotal items to revise: %d out of %d\n",
              nrow(invalid_items), nrow(results)))
} else {
  cat("All items passed CVR validation!\n")
}
