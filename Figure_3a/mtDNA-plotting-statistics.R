library(ggplot2)
library(patchwork)

setwd("C:/Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/DnaSP_statistics")
All_Hd <- read.csv("All-Hd-tot.txt", header = TRUE, row.names = NULL)
All_number_of_haplotypes <- read.csv("All-number-of-haplotypes-tot.txt", header = TRUE, row.names = NULL)
All_number_of_sites <- read.csv("All-number-of-sites-tot.txt", header = TRUE, row.names = NULL)
All_pi <- read.csv("All-pi-tot.txt", header = TRUE, row.names = NULL)
All_segregating <- read.csv("All-segregating-tot.txt", header = TRUE, row.names = NULL)
All_Strobecks <- read.csv("All-StrobecksS-tot.txt", header = TRUE, row.names = NULL)
All_K <- read.csv("All-nucl-differences-tot.txt", header = TRUE, row.names = NULL)

dataset_colors <- c("ACP Global Expanded" = "#FCAF38", "ACP Global" = "#F95335", "Ovis spp." = "#C77CFF", "Frankliniella intonsa" = "#00BFC4")
dataset_colors_num_sampl <- c("Combined (681)" = "#FCAF38", "Carlson (31)" = "#F95335", "Ovis spp. (17)" = "#C77CFF", "Frankliniella intonsa (149)" = "#00BFC4")


# Define the individual plots with bold titles and lines
p1 <- ggplot(All_Hd, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +  # Add lines connecting points
  scale_color_manual(values = dataset_colors) +
  labs(x = NULL, y = "Haplotype Diversity", title = "Haplotype Diversity") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))  # Bold title

p2 <- ggplot(All_number_of_haplotypes, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +
  scale_color_manual(values = dataset_colors) +
  labs(x = NULL, y = "Number of Haplotypes", title = "Number of Haplotypes") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))

p3 <- ggplot(All_number_of_sites, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +
  scale_color_manual(values = dataset_colors) +
  labs(x = NULL, y = "Total Nucleotide Sites", title = "Total Nucleotide Sites") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))

p4 <- ggplot(All_pi, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +
  scale_color_manual(values = dataset_colors) +
  labs(x = NULL, y = "Nucleotide Diversity", title = "Nucleotide Diversity") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))

p5 <- ggplot(All_segregating, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +
  scale_color_manual(values = dataset_colors) +
  labs(x = "Split k-mer length", y = "Total Segregating Sites", title = "Total Segregating Sites") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))

p6 <- ggplot(All_K, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  geom_line() +
  scale_color_manual(values = dataset_colors) +
  labs(x = "Split k-mer length", y = "Pairwise Nucleotide Differences (k)", title = "Pairwise Nucleotide Differences (k)") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 13),              # Larger y-axis label
        axis.text.x = element_text(face = "bold", colour = "black", size = 10),           # Bold x-axis tick labels
        axis.text.y = element_text(face = "bold", colour = "black", size = 10),
        plot.margin = unit(c(0.3, 0.8, 0.3, 0.8), "cm"))

# Combine the plots into a 3x2 grid using patchwork
combined_plot <- (p1 | p4) /
                 (p2 | p3) /
                 (p5 | p6)

combined_plot <- (p1 | p4 | p2) /
                 (p3 | p5 | p6)

# Display the combined plot
combined_plot

ggsave("combined_plot.png", 
       plot = combined_plot, 
       width = 290,        # Width in mm (max for Genome Biology)
       height = 290,       # Height in mm (approx. full page)
       units = "mm",       # Specify units as millimeters
       dpi = 300)          # Resolution at 300 DPI

ggsave("combined_plot_horiz.png", 
       plot = combined_plot, 
       width = 525,        # Width in mm (max for Genome Biology)
       height = 210,       # Height in mm (approx. full page)
       units = "mm",       # Specify units as millimeters
       dpi = 300)          # Resolution at 300 DPI

p1 <- ggplot(All_Strobecks, aes(x = X, y = Y, color = Dataset)) +
  geom_point(size = 1.5, shape = 19) +
  scale_color_manual(values = dataset_colors) +  # Specify colors
  labs(x = "Split k-mer length", y = "Strobeck's S", title = "Strobeck's S")


# Plot with specific colors for each dataset
ggplot(All_Hd, aes(x = X, y = Y, color = Dataset)) +
  geom_point() +
  labs(x = "Split k-mer length", y = "Haplotype Diversity", title = "Haplotype Diversity")
