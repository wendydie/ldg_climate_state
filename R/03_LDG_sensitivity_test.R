# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 03_LDG_sensitivity_test.R
# Last updated: 2025-01-21
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
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
quantile_colors <- c(
  "q50" = "#1F77B4",   # Blue
  "q60" = "#2CA02C",   # Green
  "q75" = "#FF7F0E",   # Orange
  "q90" = "#D62728",   # Red
  "q95" = "#9467BD"    # Purple
)
# -- 2. Define axis ranges and theme -----------------------
x_max_val <- max(time_bins$min_ma) 
x_min_val <- min(time_bins$max_ma) 

y_min_val <- min(slope_data$slope_value, na.rm = TRUE) * 1.3
y_max_val <- max(slope_data$slope_value, na.rm = TRUE) *1.3
y_text_pos <- y_min_val *0.94  # system text 

# Create a bar chart dataframe for climate state colors
climate_legend <- data.frame(
  climate_state = factor(c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"),
                         levels = c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")),
  max_ma = seq(x_max_val / 2 + 30,x_max_val / 2 - 30, length.out = 5),  # X-axis values
  min_ma = seq(x_max_val / 2 + 15,x_max_val / 2 - 45, length.out = 5),
  y = rep(y_max_val, 5),  # Ensure placement above the main plot
  climate_color = c("#005344", "#007d65", "#c8c7c7", "#b57a51", "#95484b")  # Color mapping
)

# Create a separate dataset for "Cooler Climate" and "Warmer Climate" labels
climate_labels <- data.frame(
  x = c(x_max_val / 2 + 35, x_max_val / 2 - 50),
  y = c(y_max_val * 0.8, y_max_val * 0.8),
  label = c("Cooler climate", "Warmer climate"),
  hjust = c(1,0)
)
# -- 3. Draw the figure ----------------------------------
P_north_sensitivity <- ggplot(filter(slope_data, slope_type == "Northern"), 
                  aes(x = bin_midpoint, y = slope_value, color = quantile)) +
  # Fill the "near zero" range from -0.1 to 0.1 with a light blue rectangle
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.1, ymax = 0.1,
           fill = "lightblue", alpha = 0.3) +
  # Add slope curves for different types
  geom_line(linewidth = 1) +
  geom_point(aes(shape = abs(slope_value) < 0.1), size = 2) +
  # Add climate state color bars
  geom_rect(data = climate_states, 
            aes(xmin = top, xmax = bottom, 
                ymin = y_max_val*0.9, 
                ymax = y_max_val *1.02, 
                fill = climate_color),
            linewidth = 0.3, inherit.aes = FALSE) +
  # Force ggplot2 to use the exact colors in the dataset
  scale_fill_identity(name = "Climate State", guide = "none") +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val*1.02), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8, linetype = "dashed") + 
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Geological time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val*1.02),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = 'Slope value',
    tag = "A"
  )+
  # annotate ("text", x = 540, y = y_max_val*0.85,label="A",
  #           size = 4, fontface = "bold")+
  annotate ("text", x = 80, y = y_max_val*0.8,label="Northern Hemisphere",
            size = 4, fontface = "bold")+
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = 'right',  # Remove redundant legend
    legend.title = element_text(size = 12, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 10)  # Customize legend text
  )
P_south_sensitivity <- ggplot(filter(slope_data, slope_type == "Southern"), 
                  aes(x = bin_midpoint, y = slope_value, color = quantile)) +
  # Fill the "near zero" range from -0.1 to 0.1 with a light blue rectangle
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.1, ymax = 0.1,
           fill = "lightblue", alpha = 0.3) +
  # Add slope curves for different types
  geom_line(linewidth = 1) +
  # **Use geom_smooth() for smoother lines**
  # geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, span = 0.3) +  # Loess smoothing
  geom_point(size = 2) +
  # geom_rect(data = climate_legend, 
  #           aes(xmin = min_ma, xmax = max_ma, ymin = y * 0.8, ymax = y * 0.85, fill = climate_color),
  #           inherit.aes = FALSE) +
  # # Add "Cooler Climate" & "Warmer Climate" labels (ONLY in the first facet)
  # geom_text(data = climate_labels, 
  #           aes(x = x, y = y, label = label, hjust=hjust),
  #           size = 4,  inherit.aes = FALSE) +
  # Force ggplot2 to use the exact colors in the dataset
  scale_fill_identity() +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8, linetype = "dashed") + 
  # Set x and y axis ranges
  scale_x_reverse(
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val),
    expand = c(0, 0)
  ) +
  # # Set custom colors for slope types
  # scale_color_manual(
  #   values = c("Southern" = "#E69F00", 
  #              "Northern" = "#0072B2"),
  #   labels = c("Southern", "Northern")
  # ) +
  labs(
    x = "Time (Ma)",
    y = "Slope Value",
    tag = "B"
  ) +
  # annotate ("text", x = 540, y = y_max_val*0.85,label="B",
  #           size = 4, fontface = "bold")+
  annotate ("text", x = 80, y = y_max_val*0.8,label="Southern Hemisphere",
            size = 4, fontface = "bold")+
  coord_geo(
    xlim = c(x_max_val, 0),
    pos = "bottom",
    dat = list("periods", "epochs"),
    height = unit(1.5, "lines")
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = "none"  # Remove redundant legend
  )

climate_bar <- ggplot() +
  geom_rect(data = climate_legend, 
            aes(xmin = min_ma, xmax = max_ma, ymin = 0.8, ymax = 0.85, fill = climate_color)) +
  geom_text(data = climate_labels, 
            aes(x = x, y = 0.83, label = label, hjust = hjust),
            size = 4) +
  scale_x_reverse(limits = c(x_max_val, 0)) +
  scale_fill_identity() +
  theme_void() +
  theme(
    plot.margin = margin(0, 10, 0, 10)
  ) 

# -- 4. Combine plots and print ---------------------------
slope_vTimesensitivity_plot <- ((P_north_sensitivity / P_south_sensitivity) / climate_bar) +
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
  # Fill the "near zero" range from -0.1 to 0.1 with a light blue rectangle
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -0.1, ymax = 0.1,
           fill = "lightblue", alpha = 0.3) +
  geom_boxplot(aes(fill = quantile), outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = quantile), size = 1, alpha = 0.5,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75))+
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
            aes(x = x_label, y = y_label, label =paste(hemisphere, "Hemisphere")),
            hjust = 1, vjust = 1, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
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
sb_path <- sprintf("./figures/test/%s km %slevel %s equal-area latitude bins LDG boxplot figure.jpg", 
                   params$spacing, params$level, rich_params$n_lat_bins)

ggsave(sb_path, slope_sensitivity_boxplot, width = 6, height = 8, dpi = 300)



