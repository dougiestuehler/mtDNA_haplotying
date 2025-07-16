# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots")

# Read the data from the provided table
data <- read.table("Haplotype-counts-samle-proportions-impact-all-datasets-for-R.txt", header = TRUE, sep = "\t")

# Melt the data frame to long format for ggplot2
data_long <- melt(data, id.vars = c("Dataset", "kmer"), 
                  measure.vars = c("Proportion", "Haplotypes"))

# Define custom color for kmer = 2
data_long$color <- ifelse(data_long$kmer == 2, "#F95335",
                          ifelse(data_long$kmer == 4, "#00BFC4",
                          ifelse(data_long$kmer == 6, "#C77CFF",
                          ifelse(data_long$kmer == 15, "limegreen", "black"))))


black_points <- data_long[data_long$color == "black", ]
colored_points <- data_long[data_long$color != "black", ]

# Define custom colors for each dataset
dataset_colors <- c("Combined" = "#FCAF38", "Carlson" = "#F95335", "Ovis" = "#C77CFF", "Thrips" = "#00BFC4")
dataset_colors <- c("Combined" = "grey", "Carlson" = "grey", "Ovis" = "grey", "Thrips" = "grey")


# Plot the violin plot
ggplot() +
  geom_violin(data = data_long, aes(x = Proportion, y = value, fill = Proportion), trim = FALSE) +
  geom_jitter(data = black_points, aes(x = Proportion, y = value), width = 0.2, size = 2, color = "black") +
  geom_jitter(data = colored_points, aes(x = Proportion, y = value, color = color), width = 0.2, size = 2) +
  scale_fill_manual(values = dataset_colors) +
  scale_color_identity() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Violin Plot for Different Statistics",
       x = "Proportion",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot() +
  geom_violin(data = data_long, aes(x = Dataset, y = value, fill = Dataset), trim = FALSE) +
  geom_jitter(data = black_points, aes(x = Dataset, y = value), width = 0.2, size = 1.5, color = "black") +
  geom_jitter(data = colored_points, aes(x = Dataset, y = value, color = color), width = 0.2, size = 1.5) +
  scale_fill_manual(values = dataset_colors) +
  scale_color_identity() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Violin Plot for Different Statistics",
       x = "Dataset",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Plot the violin plot
ggplot(data_long, aes(x = Dataset, y = value, fill = Dataset)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(color = color), width = 0.2, size = 1.5) +
  scale_color_identity() +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Violin Plot for Different Statistics",
       x = "Dataset",
       y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the data from the provided file
data <- read.table("Statistics-against-current-protocol-table-for-R.txt", header = TRUE, sep = "\t")

# Melt the data frame to long format for ggplot2
data_long <- melt(data, id.vars = c("Dataset", "kmer"), 
                  measure.vars = c("Haplotype_diversity", "Nucleotide_diversity", "Number_of_haplotypes"))

# Define custom colors for specific kmer values
data_long$color <- ifelse(data_long$kmer == 2 | data_long$kmer == 15, "black", NA)

# Define custom colors for each dataset
dataset_colors <- c("Combined" = "#FCAF38", "Carlson" = "#F95335", "Ovis" = "#C77CFF", "Thrips" = "#00BFC4")

# Plot the violin plot
ggplot() +
  geom_violin(data = data_long, aes(x = Dataset, y = value, fill = "grey95"), trim = FALSE) +
  geom_jitter(data = data_long[is.na(data_long$color), ], aes(x = Dataset, y = value, color = Dataset), width = 0.2, size = 1.7) +
  geom_jitter(data = data_long[!is.na(data_long$color), ], aes(x = Dataset, y = value), width = 0.2, size = 1.7, color = "black") +
  scale_fill_manual(values = dataset_colors) +
  scale_color_manual(values = dataset_colors) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Violin Plot for Different Statistics",
       x = "Dataset",
       y = "Value") +
  theme_minimal(base_size = 13) + # Increase base font size
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Increase x-axis text size
        axis.text.y = element_text(size = 10),  # Increase y-axis text size
        axis.title.x = element_text(size = 12),  # Increase x-axis title size
        axis.title.y = element_text(size = 12),  # Increase y-axis title size
        plot.title = element_text(size = 1),  # Increase plot title size
        strip.text = element_text(size = 10))  # Increase facet title size


