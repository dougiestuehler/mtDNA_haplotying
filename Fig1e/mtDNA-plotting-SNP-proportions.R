# Stacked Bar Plot

setwd("C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots")

snps <- read.table("SNP_region_percentages-table.txt", header = TRUE, sep = "\t", quote = "")

counts <- table(snps$CDS, snps$Intergenic, snps$tRNA, snps$rRNA)
barplot(counts, main="SNP Distribution by mtDNA Region",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = rownames(counts))

install.packages("wesanderson")
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(wesanderson)

# Read the data from the file
snps <- read.table("SNP_region_percentages-table.txt", header = TRUE, sep = "\t", quote = "")

# Melt the data frame to long format for ggplot2
snps_long <- melt(snps, id.vars = "Sample_names")

# Define custom colors for each region
custom_colors <- c("CDS" = "#8DD3C7", "Intergenic" = "#BEBADA", "tRNA" = "#FFFFB3", "rRNA" = "#FB8072")

# Plot the stacked bar chart
ggplot(snps_long, aes(x = Sample_names, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP Region Percentages", x = "Sample Names", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = "Set3", name = "Regions")

# Display the plot
print(ggplot(snps_long, aes(x = Sample_names, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        labs(title = "SNP Region Percentages", x = "Sample Names", y = "Percentage") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values = custom_colors, name = "Regions"))

# Load necessary libraries
library(ggplot2)
library(reshape2)

# Read the data from the file
snps <- read.table("SNP_region_percentages-table.txt", header = TRUE, sep = "\t", quote = "")

# Melt the data frame to long format for ggplot2
snps_long <- melt(snps, id.vars = "Sample_names")

# Plot the stacked bar chart
ggplot(snps_long, aes(x = Sample_names, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP Region Percentages", x = "Sample Names", y = "Percentage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3", name = "Regions")

# Display the plot
print(ggplot(snps_long, aes(x = Sample_names, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        labs(title = "SNP Region Percentages", x = "Sample Names", y = "Percentage") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_brewer(palette = "Set3", name = "Regions"))

