library(tidyr)
# Convert q50, q60, q75, q90, q95 into long format for faceted plotting
df_long <- rich_df %>%
  select(bin_midpoint, lat_band_mid_15, q50, q60, q75, q90, q95, color) %>%
  pivot_longer(cols = c(q50, q60, q75, q90, q95), 
               names_to = "quantile", values_to = "richness_value")

# Ensure consistent color mapping
rich_df <- rich_df %>%
  mutate(color = factor(color, levels = c("Northern", "Southern", "Bad hemisphere")))

# Get unique colors present in the dataset
color_levels <- c("Northern", "Southern", "Bad hemisphere")
color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Bad hemisphere" = "#D3D3D3")
color_palette <- color_palette[color_levels]  # Keep only the necessary colors

# Define linetypes for different percentiles
quantile_palette <- c("q50" = "solid", "q60" = "dashed", "q75" = "dotdash", "q90" = "twodash", "q95" = "longdash")

# Create the faceted plot
combined_rich_fig <- ggplot(rich_df, aes(x = abs_lat, y = qD_normalized, color = color)) +
  geom_point(alpha = 0.7, size = 1) +
  
  # Add raw richness percentile points
  geom_point(data = df_long, 
             aes(x = lat_band_mid_15, y = richness_value, color = color, shape = quantile),
             size = 1, inherit.aes = FALSE) +
  
  # Add Theil-Sen fitted lines for each percentile
  geom_line(data = ols_lines, 
            aes(x = lat_band_mid_15, y = fitted_values, 
                linetype = quantile, color = color), 
            linewidth = 1, inherit.aes = FALSE) +
  
  # Use viridis color palette (same for points and lines)
  scale_color_manual(name = "LDG slope",
                     values = c("Bad hemipshere" = "#D3D3D3",
                                "Northern" = "#0072B2",
                                "Southern" = "#E69F00")) +
  scale_linetype_manual(name = "Percentile Fit", values = quantile_palette) +
  scale_shape_manual(name = "Percentiles", values = c("q50" = 16, "q60" = 17, "q75" = 15, "q90" = 18, "q95" = 19)) +
  
  guides(
    color = guide_legend(title = "LDG slope"),
    linetype = guide_legend(title = "Percentile Fit"),
    shape = guide_legend(title = "Percentiles")
  ) +
  facet_wrap(~ bin_midpoint,
             labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)])),
             scales = "free_y", ncol = 6) +
  # Facet by bin_midpoint with 8 columns
  scale_y_continuous(
    breaks = function(y) {
      max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10
      mid_val <- ceiling(max_val / 2 / 10) * 10
      return(c(0, mid_val, max_val))
    }
  ) +
  # Labels
  labs(
    x = "Absolute Latitude (°)",
    y = "Generic Richness (qD)"
  ) +
  theme_minimal() +
  # Reduce spacing between facet plots
  theme(
    strip.text = element_text(size = 8, face = "bold", margin = margin(1,1,1,1)),
    strip.placement = "inside",
    panel.spacing = unit(0.01, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
    axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
    axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
    axis.title.y = element_text(size = 12, color = "black", angle = 90),  # Show y-axis title
    legend.position = "right"  # Remove redundant legend
  )

print(combined_rich_fig)
# Save the faceted figure
com_path  <- sprintf("./figures/LDG_slope_facet_%s_km_%s_quota_%s_latitude_band_combined_figure.jpg", 
                       params$spacing, params$level,rich_params$lat_band_width)
ggsave(com_path, combined_rich_fig ,  width = 8, height = 9, dpi = 300)

# Display the plot
print(p)

for (bin in unique(rich_df$bin_midpoint)) {
  
  df_bin <- rich_df %>% 
    filter(bin_midpoint == bin, !is.na(color)) %>%
    mutate(color = factor(color, levels = c("Northern", "Southern", "Bad hemisphere")))
  
  olsl_data_bin <- ols_lines %>% 
    filter(bin_midpoint == bin, !is.na(color)) %>%
    mutate(color = factor(color, levels = c("Northern", "Southern", "Bad hemisphere")))
  
  # Convert q50, q60, q75, q90, q95 into long format
  df_long <- df_bin %>%
    select(lat_band_mid_15, q50, q60, q75, q90, q95, color) %>%
    pivot_longer(cols = c(q50, q60, q75, q90, q95), 
                 names_to = "quantile", values_to = "richness_value")
  
  # Get unique colors present in the current bin
  color_levels <- unique(c(df_bin$color, olsl_data_bin$color))
  
  # Ensure we only keep colors that exist in the data
  color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Bad hemisphere" = "#D3D3D3")
  color_palette <- color_palette[color_levels]  # Only use colors that exist in the data
  
  # Define linetypes for different percentiles
  quantile_palette <- c("q50" = "solid", "q60" = "dashed", "q75" = "dotdash", "q90" = "twodash", "q95" = "longdash")
  
  p <- ggplot(df_bin, aes(x = abs_lat, y = qD_normalized, color = color)) +
    geom_point(alpha = 0.7, size = 2) +
    
    # Add raw richness percentile points
    geom_point(data = df_long, 
               aes(x = lat_band_mid_15, y = richness_value, color = color, shape = quantile),
               size = 2, inherit.aes = FALSE) +
    
    # Add Theil-Sen fitted lines for each percentile
    geom_line(data = olsl_data_bin, 
              aes(x = lat_band_mid_15, y = fitted_values, 
                  linetype = quantile, color = color), 
              linewidth = 1, inherit.aes = FALSE) +
    
    # Add OLS fit for qD_normalized
    geom_smooth(data = df_bin, aes(x = abs_lat, y = qD_normalized, color = color), 
                method = "lm", se = FALSE, linetype = "dotted", linewidth = 1.2) +
    
    scale_color_manual(name = "Hemisphere type", values = color_palette) +  
    scale_linetype_manual(name = "Percentile fit", values = quantile_palette) +
    scale_shape_manual(name = "Percentiles", values = c("q50" = 16, "q60" = 17, "q75" = 15, "q90" = 18, "q95" = 19)) +
    
    guides(
      color = guide_legend(override.aes = list(color = unname(color_palette))),
      linetype = guide_legend(title = "Percentile fit"),
      shape = guide_legend(title = "Percentiles"),
    ) +
    
    labs(
      title = sprintf("LDG Slope for Time Bin %s", bin),
      x = "Absolute Latitude (°)",
      y = "Generic richness"
    ) +
    
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  # Save the figure in a new output directory
  output_dir <- sprintf("./figures/LDG_slope_combined/%s km %squota %slatitude band", 
                        params$spacing, params$level, rich_params$lat_band_width)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s_combined.jpg", output_dir, bin)
  ggsave(file_name, p, width = 7, height = 5, dpi = 200)
}
