# Load necessary libraries
library(ggplot2)
library(patchwork)

# List of CSV files and corresponding titles
files <- c(
  "stats_no_maf_regular_coverage.csv",
  "stats_no_maf_low_coverage.csv", 
  "stats_maf_regular_coverage.csv", 
  "stats_maf_low_coverage.csv"
)

# Corresponding plot titles for each dataset
titles <- c(
  "No MAF threshold + regular coverage",
  "No MAF threshold + simulated low coverage",
  "MAF > 1% + regular coverage",
  "MAF > 1% + simulated low coverage"
)

# Base path where your CSV files are stored
base_path <- "~/Desktop/IBDGem Followup/unique joint genotypes/"

# Initialize a list to hold the plots
plot_list <- list()

# Loop over the files to read and plot each one
for (i in seq_along(files)) {
  
  # Read each CSV file into a data frame
  file_path <- paste0(base_path, files[i])
  stats_df <- read.csv(file_path, row.names = 1)
  
  # Clean up the column names by removing the 'X' prefix
  colnames(stats_df) <- gsub("^X", "", colnames(stats_df))
  
  # Generate the plot for each file
  p <- ggplot(stats_df, aes(x = as.numeric(rownames(stats_df)))) +
    
    # Line plot for Panel size = 500 individuals
    geom_line(aes(y = `500.Median`, color = "Panel size = 500 individuals"), size = 0.5) +
    
    # Fill area between 500 Min and 500 Max
    geom_ribbon(aes(ymin = `500.Min`, ymax = `500.Max`), fill = "blue", alpha = 0.1) +
    
    # Line plot for All unrelated 1kGP individuals (Panel 2500)
    geom_line(aes(y = `2500.Median`, color = "All unrelated 1kGP individuals"), size = 0.5) +
    
    # Fill area between 2500 Min and 2500 Max
    geom_ribbon(aes(ymin = `2500.Min`, ymax = `2500.Max`), fill = "#D55E00", alpha = 0.11) +
    
    # Axis labels
    labs(x = "window size", y = "fraction of unique genotypes", title = titles[i]) +
    
    # Add a legend and adjust its position
    scale_color_manual(values = c("blue", "#D55E00")) +
    
    # Disable grid lines and enable axis lines
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.key.height = unit(0.4, "cm"),
      legend.spacing.y = unit(0.1, "cm"),    
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black",linewidth =0.3),
      axis.ticks = element_line(color="black",linewidth=0.3),
      axis.title = element_text(size =8),
      plot.title = element_text(face = "bold", size=8)
    )
  
  # Remove the legend for all but the first plot
  if (i != 1) {
    p <- p + theme(legend.position = "none")
  } else {
    # For the first plot, place the legend at the upper right and make the text bigger
    p <- p + theme(
      legend.position = c(0.95, 0.3),  # Upper right corner
      legend.justification = c("right", "top"),  # Align to the top right corner
      legend.text = element_text(size = 7)  # Increase the legend text size
    )
  }
  
  # Add the plot to the list
  plot_list[[i]] <- p
}

# Combine the four plots into one
combined_plot <- (plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]]) +
  plot_layout(ncol = 2, nrow = 2) + 
  plot_annotation(tag_levels = 'A') 

# Save the combined plot as a PDF file
ggsave("~/Desktop/Final Final Figures/supp/FigS5.pdf", combined_plot, width = 8, height = 6)



######### Not using ggplot
# List of CSV files and corresponding titles
files <- c(
  "stats_no_maf_regular_coverage.csv",
  "stats_no_maf_low_coverage.csv", 
  "stats_maf_regular_coverage.csv", 
  "stats_maf_low_coverage.csv"
)

# Corresponding plot titles for each dataset
titles <- c(
  "No MAF threshold + regular coverage",
  "No MAF threshold + simulated low coverage",
  "MAF > 1% + regular coverage",
  "MAF > 1% + simulated low coverage"
)

# Base path where your CSV files are stored
base_path <- "~/Desktop/IBDGem Followup/unique joint genotypes/"

labels <- c("A", "B", "C", "D")

# Set up a 2x2 plot layout
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

#par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(0, 0, 2, 0))

# Loop over the files to read and plot each one
for (i in seq_along(files)) {
  
  # Read each CSV file into a data frame
  file_path <- paste0(base_path, files[i])
  stats_df <- read.csv(file_path, row.names = 1)
  
  # Clean up the column names by removing the 'X' prefix
  colnames(stats_df) <- gsub("^X", "", colnames(stats_df))
  
  # Extract x and y values for plotting
  x_values <- as.numeric(rownames(stats_df))
  y_median_500 <- stats_df$`500.Median`
  y_min_500 <- stats_df$`500.Min`
  y_max_500 <- stats_df$`500.Max`
  y_median_2500 <- stats_df$`2500.Median`
  y_min_2500 <- stats_df$`2500.Min`
  y_max_2500 <- stats_df$`2500.Max`
  
  # Plot for 500 individuals
  plot(x_values, y_median_500, type = "l", col = "blue", lwd = 1,
       xlab = "Window size", ylab = "Fraction of unique genotypes",
       main = titles[i], cex.main = 0.8,ylim = range(c(y_min_500, y_max_500, y_min_2500, y_max_2500)),bty="n")
  
  # Add subplot label (A, B, C, or D)
  mtext(labels[i], side = 3, line = 1, adj = -0.2, cex = 1.2)
  
  # Fill area between min and max for 500 individuals
  polygon(c(x_values, rev(x_values)), c(y_min_500, rev(y_max_500)), col = rgb(0, 0, 1, 0.1), border = NA)
  
  # Plot for 2500 individuals
  lines(x_values, y_median_2500, col = "#D55E00", lwd = 1)
  
  # Fill area between min and max for 2500 individuals
  polygon(c(x_values, rev(x_values)), c(y_min_2500, rev(y_max_2500)), col = rgb(1, 0.5, 0, 0.1), border = NA)
  
  # Add a legend to the first plot
  if (i == 1) {
    legend("bottomright", legend = c("Panel size = 500 individuals", "All unrelated 1kGP individuals"),
           col = c("blue", "#D55E00"), lty = 1, lwd = 1, bty = "n", cex = 0.8,inset = c(0.1, 0.1))
  }
}

# Save the combined plot as a PDF file
dev.copy(pdf, "~/Desktop/Final Final Figures/supp/FigS5.pdf", width = 8, height = 6)
dev.off()




