library(ggplot2)

setwd("C:/Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/DnaSP_statistics")
All_Hd <- read.csv("All-Hd-tot.txt", header = TRUE, row.names = NULL)
All_number_of_haplotypes <- read.csv("All-number-of-haplotypes-tot.txt", header = TRUE, row.names = NULL)
All_number_of_sites <- read.csv("All-number-of-sites-tot.txt", header = TRUE, row.names = NULL)
All_pi <- read.csv("All-pi-tot.txt", header = TRUE, row.names = NULL)
All_segregating <- read.csv("All-segregating-tot.txt", header = TRUE, row.names = NULL)
All_Strobecks <- read.csv("All-StrobecksS-tot.txt", header = TRUE, row.names = NULL)
All_K <- read.csv("All-nucl-differences-tot.txt", header = TRUE, row.names = NULL)

dataset_colors <- c("Combined" = "#FCAF38", "Carlson" = "#F95335", "Ovis spp." = "#C77CFF", "Frankliniella intonsa" = "#00BFC4")
dataset_colors_num_sampl <- c("Combined (681)" = "#FCAF38", "Carlson (31)" = "#F95335", "Ovis spp. (17)" = "#C77CFF", "Frankliniella intonsa (149)" = "#00BFC4")


ggplot(All_Hd, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Haplotype Diversity", title = "Haplotype Diversity")

ggplot(All_number_of_haplotypes, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Number of Haplotypes", title = "Number of Haplotypes")

ggplot(All_number_of_sites, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Total Nucleotide Sites", title = "Total Nucleotide Sites")

ggplot(All_pi, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Nucleotide Diversity", title = "Nucleotide Diversity")

ggplot(All_segregating, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Total Segregating Sites", title = "Total Segregating Sites")

ggplot(All_K, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Pairwise Nucleotide Differences (k)", title = "Pairwise Nucleotide Differences (k)")

ggplot(All_Strobecks, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Strobeck's S", title = "Strobeck's S")


# Plot with specific colors for each dataset
ggplot(All_Hd, aes(x = X, y = Y, color = Dataset)) +
  geom_point() +
  labs(x = "Split k-mer length", y = "Haplotype Diversity", title = "Haplotype Diversity")
