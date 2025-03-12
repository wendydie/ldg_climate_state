# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(viridis)
# 
# source("./R/options.R")
# # Generate OLS fitted values for visualization, grouped by percentile
# ols_lines <- rich_df %>%
#   select(bin_midpoint, lat_band_mid_15, hemisphere, q50, q60, q75, q90, q95) %>%  # Include percentiles
#   pivot_longer(cols = c(q50, q60, q75, q90, q95), names_to = "quantile", values_to = "qD_value") %>%  # Reshape to long format
#   left_join(LDG_slope, by = c("bin_midpoint", "hemisphere", "quantile")) %>%  # Join with slopes based on bin, hemisphere, and quantile
#   mutate(fitted_values = intercept + slope * lat_band_mid_15) %>%  # Compute Theil-Sen fitted values
#   filter(!is.na(fitted_values))  # Remove rows where fitted values could not be computed
# 
# 
# # Step 4: Create the scatter plot with LDG slopes------------------------
# LDG_s_plot <- ggplot(rich_df, aes(x = abs_lat, y = qD_normalized,
#                                   color = ifelse(label == "bad", "Bad hemipshere", hemisphere),
#                                   shape = hemisphere
# )) +
#   # Scatter plot points (Northern & Southern hemispheres get automatic colors)
#   geom_point(alpha = 0.7, size = 1) +
#   geom_point(data = filter(rich_df, label == "good"), 
#              aes(x = lat_band_mid_15, y = q75,shape = 'q75',
#                  color = ifelse(label == "bad", "Bad hemipshere", hemisphere)), 
#              size = 2) +
#   # Theil-Sen fitted line for Northern Hemisphere
#   geom_line(data = filter(ols_lines, quantile=='q75'),
#             aes(x = lat_band_mid_15, y = fitted_values, linetype = hemisphere, color = color),
#             linewidth = 1, inherit.aes = FALSE) +
#   # Use viridis color palette (same for points and lines)
#   scale_color_manual(name = "LDG slope",
#                      values = c("Bad hemipshere" = "#D3D3D3",
#                                 "Northern" = "#0072B2",
#                                 "Southern" = "#E69F00")) +
#   scale_shape_manual(name = "Point type",
#                      values = c("q75" = 4, "Northern" = 16, "Southern" = 17)) +
#   scale_linetype_manual(name = "Legend", values = c("Northern" = "solid", 
#                                                     "Southern" = "solid")) +
#   guides(
#     shape = guide_legend(override.aes = list(size = c(2, 3, 3))),  
#     color = guide_legend(override.aes = list(color = c("#D3D3D3", "#0072B2", "#E69F00"), shape = c(15, 16, 17))),
#     linetype = "none"
#   ) +
#   # Facet by bin_midpoint with 8 columns
#   scale_y_continuous(
#     breaks = function(y) {
#       max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10
#       mid_val <- ceiling(max_val / 2 / 10) * 10
#       return(c(0, mid_val, max_val))
#     }
#   )  +
#   facet_wrap(~ bin_midpoint,
#              labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)])),
#              scales = "free_y", ncol = 6) +
#   # Labels
#   labs(
#     x = "Absolute Latitude (°)",
#     y = "Normalized generic richness"
#   ) +
#   theme_minimal() +
#   # Reduce spacing between facet plots
#   theme(
#     strip.text = element_text(size = 8, face = "bold", margin = margin(1,1,1,1)),
#     strip.placement = "inside",
#     panel.spacing = unit(0.01, "lines"),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
#     panel.grid = element_blank(),
#     axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
#     axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
#     axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
#     axis.title.y = element_text(size = 12, color = "black", angle = 90),  # Show y-axis title
#     legend.position = "bottom"  # Remove redundant legend
#   )
# 
# print(LDG_s_plot)
# gg_path <- sprintf("./figures/%s km %squota %slatitude band LDG slopes figure.jpg", 
#                    params$spacing, params$level, rich_params$lat_band_width)
# 
# ggsave(gg_path, LDG_s_plot, width = 8, height = 9, dpi = 300)
# 
# for (bin in unique(rich_df$bin_midpoint)) {
#   
#   df_bin <- rich_df %>% 
#     filter(bin_midpoint == bin)
#   
#   for (perc in rich_params$percentiles) {
#     
#     slope_data <- LDG_slope %>% filter(bin_midpoint == bin, quantile == perc)
#     olsl_data <- ols_lines %>% filter(bin_midpoint == bin, quantile == perc)
#     
#     # Get unique color levels
#     color_levels <- unique(df_bin$color)
#     
#     # Define a fixed color palette, but only use the colors needed
#     color_palette <- c("Northern" = "#0072B2", "Southern" = "#E69F00", "Bad hemisphere" = "#D3D3D3")
#     color_palette <- color_palette[color_levels]  # Keep only the required colors
#     
#     p <- ggplot(df_bin, aes(x = abs_lat, y = qD_normalized,
#                             color = color)) +  # Use new color column
#       geom_point(alpha = 0.7, size = 2) +
#       geom_point(data = filter(df_bin, label == "good"), 
#                  aes(x = lat_band_mid_15, y = get(perc),
#                      color = color), 
#                  shape = 4, 
#                  size = 2) +
#       # Theil-Sen fitted lines (Merged for Northern & Southern Hemispheres)
#       geom_line(data = olsl_data,
#                 aes(x = lat_band_mid_15, y = fitted_values, linetype = hemisphere, color = color),
#                 linewidth = 1, inherit.aes = FALSE) +
#       scale_color_manual(name = "LDG slope", values = color_palette) +  # Dynamically adjust colors
#       scale_linetype_manual(name = "Legend", values = c("Northern" = "solid", 
#                                                         "Southern" = "solid")) +
#       guides(
#         color = guide_legend(override.aes = list(color = unname(color_palette))),
#         linetype = "none"
#       ) +
#       labs(
#         title = sprintf("LDG Slope for Bin %s - Percentile %s", bin, perc),
#         x = "Absolute Latitude (°)",
#         y = sprintf("Richness (%s)", perc)
#       ) +
#       theme_minimal() +
#       theme(
#         legend.position = "bottom",
#         legend.box = "horizontal",
#         plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
#       )
#     # Save the figure
#     output_dir <- sprintf("./figures/LDG slope per stage/%s km %squota %slatitude band", 
#                           params$spacing, params$level, rich_params$lat_band_width)
#     if (!dir.exists(output_dir)) {
#       dir.create(output_dir, recursive = TRUE)
#     }
#     file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s_percentile_%s.jpg", output_dir, bin,perc)
#     ggsave(file_name, p, width = 6, height = 5, dpi = 150)
#   }
# }

# Bin latitude into intervals
rich_df <- rich_df %>%
  mutate(lat_band_mid = floor(cell_lat / rich_params$lat_band_width) * rich_params$lat_band_width)

# Compute richness percentiles across all bins
richness_percentiles <- rich_df %>%
  group_by(bin_midpoint, lat_band_mid) %>%
  summarise(across(qD_normalized, list(q0 = ~quantile(., 0, na.rm = TRUE),
                                       q25 = ~quantile(., 0.25, na.rm = TRUE),
                                       q50 = ~quantile(., 0.50, na.rm = TRUE),
                                       q60 = ~quantile(., 0.60, na.rm = TRUE),
                                       q75 = ~quantile(., 0.75, na.rm = TRUE),
                                       q90 = ~quantile(., 0.90, na.rm = TRUE),
                                       q100 = ~quantile(., 1, na.rm = TRUE)),
                   .names = "{.fn}")) %>%
  pivot_longer(cols = starts_with("q"), names_to = "percentile", values_to = "richness") %>%
  ungroup() %>%
  mutate(percentile = factor(percentile, levels = c("q0", "q25", "q50", "q60", "q75", "q90", "q100")))  # Ensure correct legend order

# Use viridis color scheme (colorblind-friendly)
percentile_colors <- setNames(viridis(7, option = "D"), c("q0", "q25", "q50", "q60", "q75", "q90", "q100"))

# Plot species richness with faceting
LDG_fig <- ggplot() +
  geom_point(data = rich_df, aes(x = cell_lat, y = qD_normalized), color = "black", alpha = 0.5) +  # Raw richness data
  geom_line(data = richness_percentiles, aes(x = lat_band_mid, y = richness, color = percentile, group = percentile), size = 1) +
  scale_color_manual(name = "Percentile", values = percentile_colors) +
  labs(x = "Latitude",
       y = "Normalized generic richness") +
  
  # Facet by bin_midpoint with stage names
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
  )  +
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
    legend.position = "right"
  )
print(LDG_fig)
# Save the faceted figure
LDG_f_path <- sprintf("./figures/LDG_per_stage_facet_%s_km_%s_quota %slatitude_band.jpg",
                       params$spacing, params$level,rich_params$lat_band_width)
ggsave(LDG_f_path, LDG_fig, width = 8, height = 9, dpi = 300)


# for (bin in unique(rich_df$bin_midpoint)) {
#   
#   df_bin <- rich_df %>% 
#     filter(bin_midpoint == bin) %>%
#     mutate(lat_band_mid = floor(cell_lat / rich_params$lat_band_width) * rich_params$lat_band_width)  # Bin latitude into intervals
#   
#   # Compute richness percentiles
#   richness_percentiles <- df_bin %>%
#     group_by(lat_band_mid) %>%
#     summarise(across(qD_normalized, list(q0 = ~quantile(., 0, na.rm = TRUE),
#                                          q25 = ~quantile(., 0.25, na.rm = TRUE),
#                                          q50 = ~quantile(., 0.50, na.rm = TRUE),
#                                          q60 = ~quantile(., 0.60, na.rm = TRUE),
#                                          q75 = ~quantile(., 0.75, na.rm = TRUE),
#                                          q90 = ~quantile(., 0.90, na.rm = TRUE),
#                                          q100 = ~quantile(., 1, na.rm = TRUE)), 
#                      .names = "{.fn}")) %>%
#     pivot_longer(cols = starts_with("q"), names_to = "percentile", values_to = "richness") %>%
#     ungroup() %>%
#     mutate(percentile = factor(percentile, levels = c("q0", "q25", "q50", "q60", "q75", "q90", "q100")))  # Ensure correct legend order
#   
#   # Use viridis color scheme (colorblind-friendly)
#   percentile_colors <- setNames(viridis(7, option = "D"), c("q0", "q25", "q50", "q60", "q75", "q90", "q100"))
#   
#   # Plot species richness
#   p <- ggplot() +
#     geom_point(data = df_bin, aes(x = cell_lat, y = qD_normalized), color = "black", alpha = 0.5) +  # Raw richness data
#     geom_line(data = richness_percentiles, aes(x = lat_band_mid, y = richness, color = percentile, group = percentile), size = 1) +
#     scale_color_manual(name = "Percentile", values = percentile_colors) +
#     labs(title = sprintf("Genus Richness vs. Latitude (Bin: %s Ma)", bin),
#          x = "Latitude",
#          y = "Normalized richness") +
#     theme_minimal() +
#     theme(
#       legend.position = "right",
#       plot.title = element_text(size = 14, face = "bold"),
#       axis.text = element_text(size = 10),
#       axis.title = element_text(size = 12),
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
#     )
#   
#   print(p)
#   
#   # Save the figure
#   output_dir <- sprintf("./figures/LDG per stage/%s km %squota %slatitude band", 
#                         params$spacing, params$level, rich_params$lat_band_width)
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
#   file_name <- sprintf("%s/Richness_vs_Latitude_Bin_%s.jpg", output_dir, bin)
#   ggsave(file_name, p, width = 6, height = 5, dpi = 150)
# }

