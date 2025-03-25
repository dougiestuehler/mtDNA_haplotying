setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Combined_ACP_titer_analysis")

# Load required packages
library(ggplot2)

# Example data frame (replace with your actual data)
data <- read.csv('Sample-hap-titer-info.csv', header = TRUE, sep = ",")

# Load required packages
library(ggplot2)
library(dplyr)

# Add a column to categorize samples by starting letter
data <- data %>%
  mutate(Group = case_when(
    grepl("^A", Sample) ~ "A",
    grepl("^B", Sample) ~ "B",
    grepl("^C", Sample) ~ "C",
    grepl("^D", Sample) ~ "D",
    TRUE ~ "Other"
  ))

# Filter haplotypes with 3 or more samples
haplo_counts <- table(data$Haplotype)
valid_haplos <- names(haplo_counts[haplo_counts >= 3])
filtered_data <- data %>% filter(Haplotype %in% valid_haplos)

# Calculate sample counts per haplotype for labeling
sample_counts <- filtered_data %>%
  group_by(Haplotype) %>%
  summarise(n = n(), max_titer = max(Titer)) %>%
  mutate(label_pos = max_titer + 1)  # Adjust offset as needed

# Create the jittered boxplot with colored points and sample counts
plot <- ggplot(filtered_data, aes(x = Haplotype, y = Titer)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.2, alpha = 0.6, size = 2) +  # Color by Group
  geom_text(data = sample_counts, aes(x = Haplotype, y = label_pos, label = paste0("n=", n)), 
            vjust = -0.5, size = 3) +
  scale_color_manual(values = c("A" = "darkgoldenrod2", "B" = "#7CAE00", "C" = "#00BFC4", "D" = "#F8766D", "Other" = "grey")) +  # Define colors
  labs(title = "CLas Titer by Diaphorina citri mtDNA Haplotype (≥ 3 Samples)",
       x = "mtDNA Haplotype",
       y = "CLas Titer (e.g., ln(CLas_abs +1)",
       color = "Field") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export the plot as a high-resolution image
ggsave("jitter_boxplot_with_counts_colors_highres.png", 
       plot = plot, 
       width = 10, 
       height = 6, 
       dpi = 300, 
       units = "in", 
       bg = "white")

# Optional: View the plot in R
print(plot)

kruskal.test(Titer ~ Haplotype, data = filtered_data)
# Create violin plot with jittered points
ggplot(filtered_data, aes(x = Haplotype, y = Titer)) +
  geom_violin(trim = FALSE, fill = "lightgray") +
  geom_jitter(width = 0.2, alpha = 0.6, color = "red") +
  labs(title = "CLas Titer Distribution by mtDNA Haplotype (≥ 3 Samples)",
       x = "mtDNA Haplotype",
       y = "CLas Titer") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
