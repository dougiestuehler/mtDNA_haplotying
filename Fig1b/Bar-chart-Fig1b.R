# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read the data from the tab-separated file

setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots")

data <- read.table("Haplotype-counts-sample-proportions-impact-all-datasets-for-R.txt", header = TRUE, sep = "\t")


# Ensure the Dataset column is treated as a factor
data$Dataset <- as.factor(data$Dataset)
data$kmer <- as.factor(data$kmer)
data$Proportion <- as.factor(data$Proportion)

horizontal_lines <- data.frame(
  Dataset = levels(data$Dataset),
  y_intercept = c(17, 124, 17, 0)  # Replace with your specific y-intercepts
)

# Define the colors for the bars
bar_colors <- c("black", "grey70", "white")

outline_color <- "black"

line_colors <- c("red", "red", "red", "#FFFFFF00")

y_axis_limits <- data.frame(
  Dataset = levels(data$Dataset),
  y_min = c(0, 0, 0, 0),  # Replace with your specific minimum y-values
  y_max = c(17, 124, 23, 58)  # Replace with your specific maximum y-values
)

ggplot(data, aes(x = Proportion, y = Haplotypes, fill = kmer)) +
  geom_bar(stat = "identity", position = "dodge", color = outline_color) +
  facet_wrap(~ Dataset, scales = "free", labeller = as_labeller(function(x) {
    # Get the corresponding y-axis limits for each plot
    y_limit <- y_axis_limits[y_axis_limits$Dataset == x, ]
    c(paste0(x, "\n(y: ", y_limit$y_min, " - ", y_limit$y_max, ")"))
  })) +
  labs(title = "Bar Plots of Haplotypes by Dataset and Proportion",
       x = "Proportion",
       y = "Haplotypes") +
  geom_hline(data = horizontal_lines, aes(yintercept = y_intercept), color = line_colors, linetype = "dashed") +
  geom_hline(data = horizontal_lines[2, ], aes(yintercept = y_intercept), alpha = 0) +
  scale_fill_manual(values = bar_colors) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.5, "cm"),
    axis.text.x = element_text(size = 16, color = "black", margin = margin(t = -7)), # Move x-axis text closer
    axis.text.y = element_text(size = 16, color = "black", margin = margin(r = -12))  # Move y-axis text closer
  )