# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(splines)

# Input the data directly
element_data <- data.frame(
  Sample_names = c("ACP_Global", "ACP_Combined", "Thrips", "Ovis"),
  rRNA = c(12.59, 12.6, 11.17, 12.94),
  tRNA = c(9.08, 9.2, 9.14, 9.20),
  CDS = c(71.85, 72.0, 72.52, 69.26),
  Intergenic = c(6.48, 6.2, 7.17, 8.6)
)

snp_data <- data.frame(
  Sample_names = c("ACP_Global", "ACP_Combined", "Thrips", "Ovis"),
  rRNA = c(12.86, 14.29, 8.08, 7.76),
  tRNA = c(8.57, 5.71, 4.85, 4.44),
  CDS = c(75.71, 77.14, 83.37, 80.93),
  Intergenic = c(2.86, 2.86, 3.70, 6.87)
)

library(tidyverse)

#Start widened rectangles to fit more text
# Load required libraries
library(tidyverse)

# Reshape the data into long format
element_long <- pivot_longer(element_data, cols = c("rRNA", "tRNA", "CDS", "Intergenic"),
                             names_to = "Element", values_to = "Percentage")
element_long$Dimension <- "Sequence"

snp_long <- pivot_longer(snp_data, cols = c("rRNA", "tRNA", "CDS", "Intergenic"),
                         names_to = "Element", values_to = "Percentage")
snp_long$Dimension <- "SNPs"

# Combine the two data frames
combined_data <- rbind(element_long, snp_long)

# Ensure the Element column is a factor with consistent levels
combined_data$Element <- factor(combined_data$Element, 
                                levels = c("CDS", "Intergenic", "rRNA", "tRNA"))

# Ensure Dimension is a factor with correct order
combined_data$Dimension <- factor(combined_data$Dimension, levels = c("Sequence", "SNPs"))

# Define colors for the elements
element_colors <- c("tRNA" = "#FFFD70", "rRNA" = "#FE8C84", 
                    "CDS" = "#76E5C8", "Intergenic" = "#A095C8")

# Function to compute cumulative positions, spline curves, and y-axis ticks
create_proportion_data <- function(data) {
  gap <- 0.5  # Gap between bands (vertical)
  min_percentage <- 5  # Minimum percentage for narrow ribbons
  
  data <- data %>%
    group_by(Dimension) %>%
    mutate(
      Scaled_Percentage = pmax(Percentage, min_percentage),
      Scaled_Percentage = Scaled_Percentage / sum(Scaled_Percentage) * 100
    ) %>%
    ungroup()
  
  data <- data %>%
    group_by(Dimension) %>%
    arrange(Element) %>%
    mutate(
      Base_Cumulative = cumsum(Scaled_Percentage),
      Lower = lag(Base_Cumulative, default = 0),
      Cumulative = Base_Cumulative + (row_number() - 1) * gap,
      Lower = Lower + (row_number() - 1) * gap
    ) %>%
    ungroup()
  
  seq_ticks <- data %>%
    filter(Dimension == "Sequence") %>%
    arrange(Element) %>%
    mutate(Tick = Base_Cumulative + (row_number() - 1) * gap) %>%
    pull(Tick)
  seq_ticks <- c(0, seq_ticks)
  
  snp_ticks <- data %>%
    filter(Dimension == "SNPs") %>%
    arrange(Element) %>%
    mutate(Tick = Base_Cumulative + (row_number() - 1) * gap) %>%
    pull(Tick)
  snp_ticks <- c(0, snp_ticks)
  
  spline_data <- data.frame()
  x_positions <- c(0.7, 2.3)  # Ribbons span 1.6 units
  
  base_splines <- list()
  for (elem in unique(data$Element)) {
    elem_data <- data[data$Element == elem, ]
    seq_lower <- elem_data$Lower[elem_data$Dimension == "Sequence"]
    seq_upper <- elem_data$Cumulative[elem_data$Dimension == "Sequence"]
    snp_lower <- elem_data$Lower[elem_data$Dimension == "SNPs"]
    snp_upper <- elem_data$Cumulative[elem_data$Dimension == "SNPs"]
    
    x_spline <- seq(0.7, 2.3, length.out = 100)
    lower_spline <- spline(x_positions, c(seq_lower, snp_lower), xout = x_spline, method = "natural")$y
    upper_spline <- spline(x_positions, c(seq_upper, snp_upper), xout = x_spline, method = "natural")$y
    
    base_splines[[elem]] <- list(
      x = x_spline,
      Lower = lower_spline,
      Upper = upper_spline,
      Percentage_Sequence = elem_data$Percentage[elem_data$Dimension == "Sequence"],
      Percentage_SNPs = elem_data$Percentage[elem_data$Dimension == "SNPs"]
    )
  }
  
  spline_data <- data.frame()
  for (elem in unique(data$Element)) {
    x_spline <- base_splines[[elem]]$x
    lower_spline <- base_splines[[elem]]$Lower
    upper_spline <- base_splines[[elem]]$Upper
    
    elem_spline <- data.frame(
      x = x_spline,
      Lower = lower_spline,
      Upper = upper_spline,
      Element = elem,
      Percentage_Sequence = base_splines[[elem]]$Percentage_Sequence,
      Percentage_SNPs = base_splines[[elem]]$Percentage_SNPs
    )
    spline_data <- rbind(spline_data, elem_spline)
  }
  
  spline_data <- spline_data %>%
    arrange(Element, x) %>%
    group_by(x) %>%
    mutate(
      Thickness = Upper - Lower,
      Base_Lower = lag(cumsum(Thickness + gap), default = 0),
      Base_Upper = Base_Lower + Thickness
    ) %>%
    ungroup()
  
  wave_amplitude <- 2
  wave_frequency <- 1
  spline_data <- spline_data %>%
    group_by(x) %>%
    mutate(
      Wave = wave_amplitude * sin(wave_frequency * pi * (x - 0.7) / 1.6),  # Adjusted for new range
      Lower = if_else(Element == "CDS", Base_Lower, Base_Lower + Wave),
      Upper = if_else(Element == "tRNA", Base_Upper, Base_Upper + Wave),
      Lower = pmax(Lower, lag(Upper, default = -Inf) + gap),
      Upper = pmin(Upper, lead(Lower, default = Inf) - gap),
      Upper = pmax(Upper, Lower + 0.01)
    ) %>%
    ungroup() %>%
    select(-Thickness, -Base_Lower, -Base_Upper, -Wave)
  
  return(list(spline_data = spline_data, seq_ticks = seq_ticks, snp_ticks = snp_ticks))
}

# Create a list to store the plots
plots <- list()

# Create a plot for each sample
samples <- unique(combined_data$Sample_names)

for (sample in samples) {
  # Subset data for the current sample
  sample_data <- combined_data[combined_data$Sample_names == sample, ]
  
  # Compute spline data and y-axis ticks
  plot_info <- create_proportion_data(sample_data)
  plot_data <- plot_info$spline_data
  seq_ticks <- plot_info$seq_ticks
  snp_ticks <- plot_info$snp_ticks
  
  # Calculate y-axis limits
  wave_amplitude <- 2
  total_height <- 100 + (length(unique(plot_data$Element)) - 1) * 0.5  # Match gap
  y_min <- 0 - wave_amplitude
  y_max <- total_height + wave_amplitude
  
  # Define horizontal gaps between rectangles and ribbons
  left_gap <- 0.015
  right_gap <- 0.015
  
  # Calculate rectangle boundaries (widened)
  left_xmax <- 0.7 - left_gap  # 0.7 - 0.05 = 0.65
  left_xmin <- left_xmax - 0.69  # 0.65 - 0.69 = -0.04
  left_text_x <- (left_xmin + left_xmax) / 2  # (-0.04 + 0.65) / 2 = 0.305
  
  right_xmin <- 2.3 + right_gap  # 2.3 + 0.05 = 2.35
  right_xmax <- right_xmin + 0.69  # 2.35 + 0.69 = 3.04
  right_text_x <- (right_xmin + right_xmax) / 2  # (2.35 + 3.04) / 2 = 2.695
  
  # Create the proportion plot
  p <- ggplot(plot_data, aes(x = x)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Element), alpha = 0.8) +
    # Left (Sequence) side rectangles
    geom_rect(data = plot_data[plot_data$x == 0.7, ],  # Match ribbon start
              aes(xmin = left_xmin, xmax = left_xmax, ymin = Lower, ymax = Upper, fill = Element),
              alpha = 0.8) +
    geom_text(data = plot_data[plot_data$x == 0.7, ],
              aes(x = left_text_x, y = (Lower + Upper) / 2, 
                  label = paste0(Element, "\n", round(Percentage_Sequence, 1), "%")),
              hjust = 0.5, size = 3.7, color = "black", lineheight = 0.9, alpha = 1) +
    # Right (SNPs) side rectangles
    geom_rect(data = plot_data[plot_data$x == 2.3, ],  # Match ribbon end
              aes(xmin = right_xmin, xmax = right_xmax, ymin = Lower, ymax = Upper, fill = Element),
              alpha = 0.8) +
    geom_text(data = plot_data[plot_data$x == 2.3, ],
              aes(x = right_text_x, y = (Lower + Upper) / 2, 
                  label = paste0(Element, "\n", round(Percentage_SNPs, 1), "%")),
              hjust = 0.5, size = 3.7, color = "black", lineheight = 0.9, alpha = 1) +
    scale_x_continuous(breaks = c(0.7, 2.3), labels = c("Sequence", "SNPs"), 
                       limits = c(-0.05, 3.04), expand = c(0, 0)) +  # Adjusted left limit
    scale_y_continuous(breaks = seq_ticks, limits = c(y_min, y_max),
                       sec.axis = sec_axis(~., breaks = snp_ticks)) +
    scale_fill_manual(values = element_colors) +
    labs(x = NULL, y = "Percentage", 
         title = paste(sample, "Sequence vs. SNP Proportions")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.caption = element_text(size = 10, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  plots[[sample]] <- p
}
p
# Define the directory path
save_path <- "C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots"  # Windows example

# Create the directory if it doesn't exist
dir.create(save_path, showWarnings = FALSE)

for (sample in names(plots)) {
  ggsave(filename = paste0(save_path, "proportion_", sample, ".png"), 
         plot = plots[[sample]], 
         width = 8, height = 6, dpi = 300)
}

# Return the plots list
plots


#Start Vertical flow proportions
# Load required libraries
library(tidyverse)

# Reshape the data into long format
element_long <- pivot_longer(element_data, cols = c("rRNA", "tRNA", "CDS", "Intergenic"),
                             names_to = "Element", values_to = "Percentage")
element_long$Dimension <- "Sequence"

snp_long <- pivot_longer(snp_data, cols = c("rRNA", "tRNA", "CDS", "Intergenic"),
                         names_to = "Element", values_to = "Percentage")
snp_long$Dimension <- "SNPs"

# Combine the two data frames
combined_data <- rbind(element_long, snp_long)

# Ensure the Element column is a factor with consistent levels
combined_data$Element <- factor(combined_data$Element, 
                                levels = c("CDS", "Intergenic", "rRNA", "tRNA"))

# Ensure Dimension is a factor with correct order
combined_data$Dimension <- factor(combined_data$Dimension, levels = c("Sequence", "SNPs"))

# Define colors for the elements
element_colors <- c("tRNA" = "#FFFD70", "rRNA" = "#FE8C84", 
                    "CDS" = "#76E5C8", "Intergenic" = "#A095C8")

# Function to compute cumulative positions, spline curves, and x-axis ticks
create_proportion_data <- function(data) {
  gap <- 0.5  # Gap between bands (now horizontal)
  min_percentage <- 5  # Minimum percentage for narrow ribbons
  
  # Scale percentages to ensure minimum width, then normalize to sum to 100
  data <- data %>%
    group_by(Dimension) %>%
    mutate(
      Scaled_Percentage = pmax(Percentage, min_percentage),
      Scaled_Percentage = Scaled_Percentage / sum(Scaled_Percentage) * 100
    ) %>%
    ungroup()
  
  # Use Scaled_Percentage for stacking (now along x-axis)
  data <- data %>%
    group_by(Dimension) %>%
    arrange(Element) %>%
    mutate(
      Base_Cumulative = cumsum(Scaled_Percentage),
      Lower = lag(Base_Cumulative, default = 0),
      Cumulative = Base_Cumulative + (row_number() - 1) * gap,
      Lower = Lower + (row_number() - 1) * gap
    ) %>%
    ungroup()
  
  # Compute x-axis tick positions (cumulative boundaries)
  seq_ticks <- data %>%
    filter(Dimension == "Sequence") %>%
    arrange(Element) %>%
    mutate(Tick = Base_Cumulative + (row_number() - 1) * gap) %>%
    pull(Tick)
  seq_ticks <- c(0, seq_ticks)
  
  snp_ticks <- data %>%
    filter(Dimension == "SNPs") %>%
    arrange(Element) %>%
    mutate(Tick = Base_Cumulative + (row_number() - 1) * gap) %>%
    pull(Tick)
  snp_ticks <- c(0, snp_ticks)
  
  # Create a data frame for spline interpolation (y-axis now controls flow)
  spline_data <- data.frame()
  y_positions <- c(0.5, 2.5)  # Top (Sequence) to bottom (SNPs)
  
  # Compute base splines (along y-axis)
  base_splines <- list()
  for (elem in unique(data$Element)) {
    elem_data <- data[data$Element == elem, ]
    seq_lower <- elem_data$Lower[elem_data$Dimension == "Sequence"]
    seq_upper <- elem_data$Cumulative[elem_data$Dimension == "Sequence"]
    snp_lower <- elem_data$Lower[elem_data$Dimension == "SNPs"]
    snp_upper <- elem_data$Cumulative[elem_data$Dimension == "SNPs"]
    
    y_spline <- seq(0.5, 2.5, length.out = 100)  # Vertical flow
    lower_spline <- spline(y_positions, c(seq_lower, snp_lower), xout = y_spline, method = "natural")$y
    upper_spline <- spline(y_positions, c(seq_upper, snp_upper), xout = y_spline, method = "natural")$y
    
    base_splines[[elem]] <- list(
      y = y_spline,  # Swapped from x to y
      Lower = lower_spline,
      Upper = upper_spline,
      Percentage_Sequence = elem_data$Percentage[elem_data$Dimension == "Sequence"],
      Percentage_SNPs = elem_data$Percentage[elem_data$Dimension == "SNPs"]
    )
  }
  
  # Stack the ribbons (now horizontally)
  spline_data <- data.frame()
  for (elem in unique(data$Element)) {
    y_spline <- base_splines[[elem]]$y
    lower_spline <- base_splines[[elem]]$Lower
    upper_spline <- base_splines[[elem]]$Upper
    
    elem_spline <- data.frame(
      y = y_spline,
      Lower = lower_spline,
      Upper = upper_spline,
      Element = elem,
      Percentage_Sequence = base_splines[[elem]]$Percentage_Sequence,
      Percentage_SNPs = base_splines[[elem]]$Percentage_SNPs
    )
    spline_data <- rbind(spline_data, elem_spline)
  }
  
  # Stack with gaps (along x-axis)
  spline_data <- spline_data %>%
    arrange(Element, y) %>%
    group_by(y) %>%
    mutate(
      Thickness = Upper - Lower,
      Base_Lower = lag(cumsum(Thickness + gap), default = 0),
      Base_Upper = Base_Lower + Thickness
    ) %>%
    ungroup()
  
  # Apply wavy effect (now horizontal waves)
  wave_amplitude <- 2
  wave_frequency <- 1
  spline_data <- spline_data %>%
    group_by(y) %>%
    mutate(
      Wave = wave_amplitude * sin(wave_frequency * pi * (y - 0.5) / 2),  # Wave along y
      Lower = if_else(Element == "CDS", Base_Lower, Base_Lower + Wave),
      Upper = if_else(Element == "tRNA", Base_Upper, Base_Upper + Wave),
      Lower = pmax(Lower, lag(Upper, default = -Inf) + gap),
      Upper = pmin(Upper, lead(Lower, default = Inf) - gap),
      Upper = pmax(Upper, Lower + 0.01)
    ) %>%
    ungroup() %>%
    select(-Thickness, -Base_Lower, -Base_Upper, -Wave)
  
  return(list(spline_data = spline_data, seq_ticks = seq_ticks, snp_ticks = snp_ticks))
}

# Create a list to store the plots
plots <- list()

# Create a plot for each sample
samples <- unique(combined_data$Sample_names)

for (sample in samples) {
  # Subset data for the current sample
  sample_data <- combined_data[combined_data$Sample_names == sample, ]
  
  # Compute spline data and x-axis ticks
  plot_info <- create_proportion_data(sample_data)
  plot_data <- plot_info$spline_data
  seq_ticks <- plot_info$seq_ticks
  snp_ticks <- plot_info$snp_ticks
  
  # Calculate x-axis limits (formerly y-axis)
  wave_amplitude <- 2
  total_width <- 100 + (length(unique(plot_data$Element)) - 1) * 0.5  # Match gap = 0.5
  x_min <- 0 - wave_amplitude
  x_max <- total_width + wave_amplitude
  
  # Define vertical gaps between rectangles and ribbons
  top_gap <- 0.05
  bottom_gap <- 0.05
  
  # Calculate rectangle boundaries (now on y-axis)
  top_ymin <- 0.12
  top_ymax <- 0.5 - top_gap  # 0.45
  top_text_y <- (top_ymin + top_ymax) / 2
  
  bottom_ymin <- 2.5 + bottom_gap  # 2.55
  bottom_ymax <- bottom_ymin + (2.88 - 2.51)
  bottom_text_y <- (bottom_ymin + bottom_ymax) / 2
  
  # Create the proportion plot (flipped)
  p <- ggplot(plot_data, aes(y = y)) +  # Swap x to y
    geom_ribbon(aes(xmin = Lower, xmax = Upper, fill = Element), alpha = 0.8) +  # Swap ymin/ymax to xmin/xmax
    # Top (Sequence) rectangles
    geom_rect(data = plot_data[plot_data$y == 0.5, ],
              aes(ymin = top_ymin, ymax = top_ymax, xmin = Lower, xmax = Upper, fill = Element),
              alpha = 2) +
    geom_text(data = plot_data[plot_data$y == 0.5, ],
              aes(y = top_text_y, x = (Lower + Upper) / 2, 
                  label = paste0(Element, "\n", round(Percentage_Sequence, 1), "%")),
              hjust = 0.5, vjust = 0.5, size = 3.7, color = "black", lineheight = 0.9, alpha = 0, angle = 90) +  # Rotate text
    # Bottom (SNPs) rectangles
    geom_rect(data = plot_data[plot_data$y == 2.5, ],
              aes(ymin = bottom_ymin, ymax = bottom_ymax, xmin = Lower, xmax = Upper, fill = Element),
              alpha = 2) +
    geom_text(data = plot_data[plot_data$y == 2.5, ],
              aes(y = bottom_text_y, x = (Lower + Upper) / 2, 
                  label = paste0(Element, "\n", round(Percentage_SNPs, 1), "%")),
              hjust = 0.5, vjust = 0.5, size = 3.7, color = "black", lineheight = 0.9, alpha = 0, angle = 90) +  # Rotate text
    scale_y_continuous(breaks = c(0.5, 2.5), labels = c("Sequence", "SNPs"), 
                       limits = c(0, 3), expand = c(0, 0)) +  # Now vertical axis
    scale_x_continuous(breaks = seq_ticks, limits = c(x_min, x_max),  # Now horizontal axis
                       sec.axis = sec_axis(~., breaks = snp_ticks)) +
    scale_fill_manual(values = element_colors) +
    labs(y = NULL, x = "Percentage",  # Swap x and y labels
         title = paste(sample, "Sequence vs. SNP Proportions")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 12, face = "bold"),  # Now vertical (Sequence/SNPs)
          axis.text.x = element_text(size = 12),                # Now horizontal (percentages)
          axis.title.x = element_text(size = 14),               # Now horizontal
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.caption = element_text(size = 10, hjust = 0),
          panel.grid.major.y = element_blank(),  # Now vertical grid
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),  # Now horizontal grid
          panel.grid.minor.x = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
  plots[[sample]] <- p
}

# Define the directory path (fixed double slashes for Windows)
save_path <- "C://Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Rplots"

# Create the directory if it doesn't exist
dir.create(save_path, showWarnings = FALSE)

for (sample in names(plots)) {
  ggsave(filename = paste0(save_path, "/proportion_vert_", sample, ".png"), 
         plot = plots[[sample]], 
         width = 6, height = 8, dpi = 300)  # Swapped width and height
}

# Return the plots list
plots
