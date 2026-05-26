# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 03_LDG_sensitivity_test.R
# Last updated: 2025-10-15
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(ggplot2)
library(dplyr)
#-------------------------------------------------------------------------
# Drawing the slope figures of LDG of time-bins in one figure-------------
# -- 1. Prepare slope data for plotting ---------------------
slope_data <- slope_cli_df %>%
  mutate(
    slope_type = case_when(
      hemisphere == "Northern" ~ "Northern",
      hemisphere == "Southern" ~ "Southern",
      TRUE ~ NA_character_  # Handle unexpected cases
    ),
    bin_midpoint = as.numeric(as.character(bin_midpoint)),
    slope_value = ifelse(label == "bad", NA, slope)
  )
# slope_data <- LDG_slope %>%
#   mutate(bin_midpoint = as.numeric(as.character(bin_midpoint)),
#          slope = ifelse(label == "bad", NA, slope))

# Define colors for different quantiles
# quantile_colors <- c(
#   "q50" = "#1F77B4",   # Blue
#   "q60" = "#2CA02C",   # Green
#   "q75" = "#FF7F0E",   # Orange
#   "q90" = "#D62728",   # Red
#   "q95" = "#9467BD"    # Purple
# )
# quantile_colors <- setNames(gray.colors(5, start = 0.8, end=0.1),
#                             c("q50", "q60", "q75", "q90", "q95"))
# -- 2. Define axis ranges and theme -----------------------
x_max_val <- max(time_bins$min_ma) 
x_min_val <- min(time_bins$max_ma) 

y_min_val <- min(slope_data$slope_value, na.rm = TRUE) * 1.3
y_max_val <- max(slope_data$slope_value, na.rm = TRUE) *1.3
x_text_pos <- 15
y_text_pos <- y_min_val *0.94  # system text 

# -- 3. Draw the figure ----------------------------------
P_north_sensitivity <- ggplot(filter(slope_data, slope_type == "Northern"), 
                  aes(x = bin_midpoint, y = slope_value, color = quantile)) +
  geom_vline(
    xintercept = major_boundaries,
    color = "black", linewidth = 0.4, alpha = 0.8
  ) +
  # Add slope curves for different types
  geom_line(alpha = 0.6, linewidth = 0.8) +
  geom_point(
    aes(shape = abs(slope_value) < 0.1),
    size = 2, stroke = 0.55, alpha = 0.6
  ) +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
  # scale_color_manual(
  #   name = "Percentile",
  #   values = quantile_colors
  # ) +
  # Add climate state color bars
  geom_rect(
    data = climate_states, 
    aes(xmin = bottom, xmax = top, ymin = y_max_val, ymax = y_max_val * 1.12),
    fill = I(climate_states$climate_color), 
    color = "black", linewidth = 0.3, inherit.aes = FALSE
  ) +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val*1.12), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8, linetype = "dashed") + 
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val*1.12),
    breaks = seq(-3, 2, 1),
    expand = c(0, 0)
  )  +
  labs(
    x = NULL,
    y = 'Slope value',
    tag = "A",
    subtitle = "Northern Hemisphere"
  ) +
  guides(
    color = guide_legend(title = "Percentile"),
    linetype = "none",
    shape = "none"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = 'none',  # Remove redundant legend
    legend.title = element_text(size = 12, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 10)  # Customize legend text
  )

P_south_sensitivity <- ggplot(filter(slope_data, slope_type == "Southern"), 
                  aes(x = bin_midpoint, y = slope_value, color = quantile)) +
  geom_vline(
    xintercept = major_boundaries,
    color = "black", linewidth = 0.4, alpha = 0.8
  )  +
  # Add slope curves for different types
  geom_line(alpha = 0.6,linewidth = 1) +
  geom_point(
    aes(shape = abs(slope_value) < 0.1),
    size = 2, alpha = 0.6,stroke = 0.55
  ) +
  scale_shape_manual(values = c(`TRUE` = 1, `FALSE` = 16)) +
  # scale_color_manual(
  #   name = "Percentile",
  #   values = quantile_colors
  # ) +
  scale_fill_identity() +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val*1.12), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8, linetype = "dashed") + 
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val*1.12),
    breaks = seq(-3, 2, 1),
    expand = c(0, 0)
  )  +
  labs(
    x = "Time (Ma)",
    y = "Slope value",
    tag = "B",
    subtitle = "Southern Hemisphere"
  ) +
  guides(
    color = guide_legend(title = "Percentile"),
    linetype = "none",
    shape = "none"
  ) +
  coord_geo(
    xlim = c(x_max_val, 0),
    pos = "bottom",
    dat = list("periods", "epochs"),
    height = unit(1.5, "lines")
  ) + 
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = "bottom"  # Remove redundant legend
  )

climate_bar <- ggplot(climate_legend, aes(x = 0, y = 0, fill = climate_state)) +
  # Invisible points to trigger the legend (not plotted on canvas)
  geom_point(shape = 22, size = 5, alpha = 0) +
  scale_fill_manual(
    values = climate_colors,
    name   = "Climate state",
    drop   = FALSE
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1, shape = 22, size = 5, colour = NA) # show color blocks in legend
    )
  ) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 12),
    plot.margin     = margin(0, 10, 0, 10)
  )
# Extract only the legend as a grob
climate_legend_grob <- cowplot::get_legend(climate_bar)

# -- 4. Combine plots and print ---------------------------
slope_vTimesensitivity_plot <- ((P_north_sensitivity / P_south_sensitivity) / climate_legend_grob) +
  plot_layout(heights = c(10, 10, 1)) &
  theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing = unit(0, "cm")
  )

# Display the plot
print(slope_vTimesensitivity_plot)
# Save the boxplot figure
svTimeall_path <- sprintf("./figures/test/%s km %slevel %s equal-area latitude bins LDG time series figure.jpg", 
                          params$spacing, params$level, rich_params$n_lat_bins)
ggsave(svTimeall_path, slope_vTimesensitivity_plot, width = 8, height = 8, dpi = 300)
# Save the boxplot figure
svTimeall_path <- sprintf("./figures/test/%s km %slevel %s equal-area latitude bins LDG time series figure.pdf", 
                          params$spacing, params$level, rich_params$n_lat_bins)
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


# Create boxplot showing LDG slope distribution across climate states
slope_sensitivity_boxplot <- ggplot(slope_boxplot_data, aes(x = climate_state, y = slope, fill = quantile)) +
  # Fill the "near zero" range from -0.1 to 0.1 with a light blue rectangle
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.1, ymax = 0.1,
           fill = "lightblue", alpha = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8,linetype = "dashed") + 
  geom_boxplot(aes(fill = quantile), outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = quantile), size = 1, alpha = 0.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75))+
  # scale_fill_manual(name = "Percentile", values = quantile_colors) +  # Assign fill colors
  # scale_color_manual(name = "Percentile", values = quantile_colors) +  # Assign point colors
  
  # Facet by hemisphere
  facet_wrap(~ hemisphere, ncol=1,labeller = labeller(
    hemisphere = c("Northern" = "Northern Hemisphere", "Southern" = "Southern Hemisphere")
  )) +  
  
  # Add custom facet labels inside the plot (Top-right position)
  geom_text(data = slope_boxplot_data %>%
              group_by(hemisphere) %>%
              summarise(x_label = 'Hothouse',  # Position at the rightmost point
                        y_label = y_max_val*0.9),  # Position above max slope
            aes(x = x_label, y = y_label, label =paste(hemisphere, "Hemisphere")),
            hjust = 0.7, vjust = 0.5, size = 4, inherit.aes = FALSE) +
  # Titles and formatting
  labs(
    x = "Climate state",
    y = "Slope value"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(color = "black", linewidth=0.6),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.margin = margin(0,0,0,0),
    legend.box.spacing = unit(0.1, "cm")
    
  ) +
  guides(fill = guide_legend(direction = "horizontal"))

# Display the plot
print(slope_sensitivity_boxplot)

# Save the boxplot figure
sb_path <- sprintf("./figures/test/%s km %slevel %s equal-area latitude bins LDG boxplot figure.jpg", 
                   params$spacing, params$level, rich_params$n_lat_bins)

ggsave(sb_path, slope_sensitivity_boxplot, width = 6.5, height = 8, dpi = 300)



