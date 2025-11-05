# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope_sensitivity_test.R
# Last updated: 2025-10-15
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
# -----------------------------------------------------------------------
# Define percentiles to analyze
percentiles <- c("q50", "q60", "q75", "q90", "q95")

# Compute slopes separately for qD_normalized and qD for Northern and Southern Hemispheres
compute_slopes <- function(df, richness_column) {
  results <- list()
  
  for (perc in percentiles) {
    slope_data <- df %>%
      group_by(bin_midpoint, hemisphere) %>%
      filter(n() > 2, length(unique(abs_lat_bin_mid)) > 1) %>%  # Ensure enough data points
      summarise(
        model = list(lm(!!sym(perc) ~ abs_lat_bin_mid, data = cur_data())), 
        .groups = "drop"
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract Theil-Sen slope
        intercept = coef(model)[1],  # Extract Theil-Sen intercept
        quantile = perc,  # Store percentile for reference
        richness_type = richness_column  # Store richness type (qD or qD_normalized)
      ) %>%
      select(-model)
    
    results[[perc]] <- slope_data
  }
  
  return(bind_rows(results))
}
# -----------------------------------------------------------------------
# Compute slopes for both qD_normalized and qD
slope_qD_normalized <- compute_slopes(rich_df, "qD_normalized")
slope_qD <- compute_slopes(rich_df, "qD")

# Combine the percentile (quantile) and richness type (qD or qD_normalized)
LDG_slope_all <- bind_rows(slope_qD_normalized, slope_qD) %>%
  mutate(slope_type = paste0(quantile, ifelse(richness_type == "qD_normalized", "_Norm", "")))

# Convert `slope_type` into wide format for correlation matrix
slope_data_wide <- LDG_slope_all %>%
  select(bin_midpoint, hemisphere, slope_type, slope) %>%
  pivot_wider(names_from = slope_type, values_from = slope)

# Drop bin_midpoint column for correlation calculation
slope_matrix <- slope_data_wide %>% 
  select(-bin_midpoint) %>%
  mutate(hemisphere = ifelse(hemisphere == "Northern", "N", "S"))

# Create the correlation scatterplot matrix
correlation_plot <- ggpairs(
  slope_matrix, 
  columns = 2:ncol(slope_matrix),  # Exclude "hemisphere" from numeric calculations
  aes(color = hemisphere, alpha = 0.7),  # Color points by hemisphere
  upper = list(continuous = wrap("cor", size = 2)),  # Display correlation coefficient in upper panel
  lower = list(continuous = wrap("points", alpha = 0.5, size = 1)),  # Scatter plots in lower panel
  diag = list(continuous = wrap("densityDiag"))  # Density plots on diagonal
) +
  theme_minimal() +  # Remove default ggplot2 background
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Keep black borders for subplots
    strip.background = element_blank(),  # Remove background (border) for facet labels
    strip.text = element_text(size = 8, face = "bold"),  # Reduce facet label text size
    axis.text = element_text(size = 6),  # Reduce axis text size
    axis.title = element_text(size = 8),  # Reduce axis title size
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Ensure no title border
    plot.background = element_blank()  # Ensure no background box for the whole plot
  )

cor_path <- sprintf("./figures/test/%skm %squota %slatitude band correlation plot.jpg", 
                    params$spacing, params$level, rich_params$n_lat_bins)
# Save the figure
ggsave(cor_path, correlation_plot, width = 8, height = 8, dpi = 300)
