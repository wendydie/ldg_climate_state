library(ggplot2)
library(dplyr)

# Filter out "bad" data and ensure bin_midpoint is numeric
slope_data <- LDG_slope %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint)),
         slope = ifelse(label == "bad", NA, slope))

# Define colors for different quantiles
quantile_colors <- c(
  "q50" = "#1F77B4",   # Blue
  "q60" = "#2CA02C",   # Green
  "q75" = "#FF7F0E",   # Orange
  "q90" = "#D62728",   # Red
  "q95" = "#9467BD"    # Purple
)

# Create a line plot to show slope changes over time
slope_vTimesensitivity_plot <- ggplot(slope_data, aes(x = bin_midpoint, y = slope, color = quantile)) +
  geom_line(linewidth = 1) +  # Draw trend lines
  geom_point(size = 2) +  # Add points for each bin
  # Add "Cooler Climate" & "Warmer Climate" labels (ONLY in the first facet)
  geom_rect(data = climate_legend, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y * 0.8, ymax = y * 0.85, fill = climate_color),
            inherit.aes = FALSE) +
  geom_text(data = climate_labels, 
            aes(x = x, y = y, label = label, hjust=hjust),
            size = 4,  inherit.aes = FALSE) +
  # Add climate state color bars
  geom_rect(data = climate_states, 
            aes(xmin = top, xmax = bottom, ymin = y_max_val*0.9, ymax = y_max_val, fill = climate_color),
            linewidth = 0.3, inherit.aes = FALSE) +
  # Add geological system color bands
  geom_rect(data = system_labels, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y_min_val*0.9, ymax = y_min_val, fill = fill_color),
            color = "black", linewidth = 0.3, inherit.aes = FALSE) +  
  # Add geological stage color bands
  geom_rect(data = time_bins, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y_min_val*0.9, ymax = y_min_val*0.85, fill = stageCol),
            color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  # Force ggplot2 to use the exact colors in the dataset
  scale_fill_identity() +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8) + 
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Geological time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val),
    expand = c(0, 0)
  ) +
  
  # Add geological system names above the color bands
  geom_text(
    data = system_labels, 
    aes(x = mid_ma, y = y_text_pos, label = system),
    angle = 0, 
    vjust = 0.5,      
    hjust = 0.5,      
    size = 4, inherit.aes = FALSE
  ) +
  # # Set custom colors for slope types
  # scale_color_manual(
  #   values = c("Southern" = "#E69F00", 
  #              "Northern" = "#0072B2"),
  #   labels = c("Southern", "Northern")
  # ) +
  # **Facet the plot by slope type (separate subplots for Northern and Southern slopes)**
  facet_wrap(~ hemisphere, scales = "free_y", ncol=1) +
  
  # Add custom facet labels inside the plot (Top-right position)
  geom_text(data = slope_data %>%
              group_by(hemisphere) %>%
              summarise(x_label = 20,  # Position at the rightmost point
                        y_label = y_max_val*0.85),  # Position above max slope
            aes(x = x_label, y = y_label, label = hemisphere),
            hjust = 1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE) +
  # Other aesthetic improvements
  labs(
    x = "Geological Time (Ma)",
    y = "Slope Value"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = "right"  # Remove redundant legend
  )
# Display the plot
print(slope_vTimesensitivity_plot)

# Save the boxplot figure
svTimeall_path <- sprintf("./figures/test/%s km %slevel %s latitude band LDG time series figure.jpg", 
                   params$spacing, params$level, rich_params$lat_band_width)

ggsave(svTimeall_path, slope_vTimesensitivity_plot, width = 8, height = 8, dpi = 300)

# Filter out "bad" data and ensure bin_midpoint is numeric
slope_data <- LDG_slope %>%
  filter(label != "bad") %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint)))

# Merge slope data with climate states and ensure correct order of climate states
slope_boxplot_data <- slope_data %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid")) %>%
  mutate(
    climate_state = factor(climate_state, levels = c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"))
  )

# Define quantile colors
quantile_colors <- c(
  "q50" = "#1F77B4",   # Blue
  "q60" = "#2CA02C",   # Green
  "q75" = "#FF7F0E",   # Orange
  "q90" = "#D62728",   # Red
  "q95" = "#9467BD"    # Purple
)

# Create boxplot showing LDG slope distribution across climate states
slope_sensitivity_boxplot <- ggplot(slope_boxplot_data, aes(x = climate_state, y = slope, fill = quantile)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Draw boxplot without outliers
  geom_jitter(aes(color = quantile), width = 0.2, size = 1, alpha = 0.5) +  # Add jitter points for data visualization
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) + 
  scale_fill_manual(name = "Quantile", values = quantile_colors) +  # Assign fill colors
  scale_color_manual(name = "Quantile", values = quantile_colors) +  # Assign point colors
  
  # Facet by hemisphere
  facet_wrap(~ hemisphere, ncol=1) +  
  
  # Add custom facet labels inside the plot (Top-right position)
  geom_text(data = slope_boxplot_data %>%
              group_by(hemisphere) %>%
              summarise(x_label = 'Hothouse',  # Position at the rightmost point
                        y_label = y_max_val*0.85),  # Position above max slope
            aes(x = x_label, y = y_label, label = hemisphere),
            hjust = 1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE) +
  # Titles and formatting
  labs(
    x = "Climate State",
    y = "Slope Value",
    title = "LDG Slope Distributions Across Climate States"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "bottom"  # Remove redundant legend
    
  ) +
  guides(fill = guide_legend(direction = "horizontal"))

# Display the plot
print(slope_sensitivity_boxplot)

# Save the boxplot figure
sb_path <- sprintf("./figures/test/%s km %slevel %s latitude band LDG boxplot figure.jpg", 
                   params$spacing, params$level, rich_params$lat_band_width)

ggsave(sb_path, slope_sensitivity_boxplot, width = 6, height = 8, dpi = 300)
