#geowendywen@outlook.com

library(dplyr)
library(tidyr)
library(ggplot2)
source("./R/options.R")
# 02 calculating the slope of LDG in time bins (N,S,combined)---------------

# Read dataset--------------------------------------------------------------
rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv", 
                            params$spacing, params$level))
time_bins <- readRDS("./data/time_bins.RDS")

# Step 1 : Data filtering --------------------------------------------------
# Our time period starts from 486.8500 because the climate state data begins at 486.85.
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>% 
  filter(bin_midpoint < 486.8500)%>%
  mutate(bin_midpoint = factor(bin_midpoint, 
                               levels = sort(unique(bin_midpoint), decreasing = TRUE)))

# If we use the Multiton rate, defined as (Mrate = (Sobs - p1)/p1), as the filtering rule, then (Mrate_filter_no) is greater than 0. 
# However, in the main text, we do not apply the Multiton rate for filtering, so Mrate_filter_no is 0.
rich_df <- rich_df %>%
  filter(Mrate >= rich_params$Mrate_filter_no) # Mrate >= 0.3

# In the main text, we apply the collection for filtering, so the col_filter_no is greater than 0.
# Here, we test different levels of filtering in our work, including 0, 1, 2, 3, and 4.
rich_df <- rich_df %>%
  filter(collection_no >= rich_params$col_filter_no)

# The richness values are normalized to 100 in each stage by dividing the richness in each cell by the maximum richness within that stage.
rich_df <- rich_df %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized = qD*100 / max(qD)) %>%
  ungroup()

# Convert latitude to absolute values and classify hemispheres.
rich_df <- rich_df %>%
  mutate(
    abs_lat = abs(cell_lat),  # Convert latitude to absolute value
    hemisphere = case_when(
      cell_lat >= 0 ~ "Northern",  # good cells & cell_lat >= 0 → Northern
      cell_lat < 0 ~ "Southern"  # good cells & cell_lat < 0 → Southern
    ))

# Divide the cells into 30-degree latitude bands.
expected_lat_bands <- seq(rich_params$lat_band_width / 2, 90, rich_params$lat_band_width)  # Midpoints
rich_df <- rich_df %>%
  mutate(lat_band_mid = (floor(abs_lat / rich_params$lat_band_width) * 
                           rich_params$lat_band_width) + (rich_params$lat_band_width / 2))

# Check if each bin_midpoint & hemisphere contains all required lat_band_mid values.
rich_df <- rich_df %>%
  group_by(bin_midpoint, hemisphere) %>%
  mutate(
    label = case_when(
      all(expected_lat_bands %in% lat_band_mid) ~ "good",
      TRUE ~ "bad"
    )
  ) %>%
  ungroup()

# Step 2: Calculate LDG slopes-------------------------------------------
# Calculate the richness at the 25th, 50th, 75th, and 95th percentiles for testing. However, in the paper, I will only use the median results.
rich_df <- rich_df %>%
  group_by(bin_midpoint, hemisphere, lat_band_mid) %>%
  mutate(
    q25 = quantile(qD_normalized, 0.25, na.rm = TRUE),
    q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
    q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
    q95 = quantile(qD_normalized, 0.95, na.rm = TRUE)
  ) %>%
  ungroup()

# The function used to calculate the LDG slope.
calculate_slope <- function(df, percentile) {
  df %>%
    group_by(bin_midpoint) %>%
    summarise(
      slope = coef(lm(!!sym(percentile) ~ lat_band_mid, data = cur_data()))[2],
      quantile = percentile, 
      .groups = "drop"
    )
}

northern_slope <- bind_rows(
  calculate_slope(filter(rich_df, hemisphere == "Northern"), "q25"),
  calculate_slope(filter(rich_df, hemisphere == "Northern"), "q50"),
  calculate_slope(filter(rich_df, hemisphere == "Northern"), "q75"),
  calculate_slope(filter(rich_df, hemisphere == "Northern"), "q95")
) %>%
  mutate(hemisphere = "Northern")

southern_slope <- bind_rows(
  calculate_slope(filter(rich_df, hemisphere == "Southern"), "q25"),
  calculate_slope(filter(rich_df, hemisphere == "Southern"), "q50"),
  calculate_slope(filter(rich_df, hemisphere == "Southern"), "q75"),
  calculate_slope(filter(rich_df, hemisphere == "Southern"), "q95")
) %>%
  mutate(hemisphere = "Southern")

# Combine the results
LDG_slope <- bind_rows(northern_slope, southern_slope)
LDG_slope$label <- rich_df$label[match(
  paste(LDG_slope$bin_midpoint, LDG_slope$hemisphere),
  paste(rich_df$bin_midpoint, rich_df$hemisphere)
)]

rm(northern_slope, southern_slope)
# View the LDG slope
View(LDG_slope)

# Step 4: Create the scatter plot with LDG slopes------------------------
LDG_s_plot <- ggplot(rich_df, aes(x = abs_lat, y = qD_normalized,
                    color = ifelse(label == "bad", "Bad hemipshere", hemisphere),
                    shape = hemisphere
                    )) +
  # Scatter plot points (Northern & Southern hemispheres get automatic colors)
  geom_point(alpha = 0.7, size = 1) +
  # Overlay precomputed LDG slopes with increased line thickness
  geom_smooth(method = "lm", se = FALSE, aes(x = abs_lat, y = q50, linetype = "Northern Slope",
                                             color = ifelse(label == "bad", "Bad hemipshere", hemisphere)),
              linewidth = 1, data = filter(rich_df, hemisphere == "Northern"), inherit.aes = FALSE) +  # Northern slope
  geom_smooth(method = "lm", se = FALSE, aes(x = abs_lat, y = q50, linetype = "Southern Slope",
                                             color = ifelse(label == "bad", "Bad hemipshere", hemisphere)),
              linewidth = 1, data = filter(rich_df, hemisphere == "Southern"), inherit.aes = FALSE) +  # Southern slope
  # Define legend for slope lines
  # Use viridis color palette (same for points and lines)
  scale_color_manual(name = "LDG slope",
                     values = c("Bad hemipshere" = "#D3D3D3",
                                "Northern" = "#0072B2",
                                "Southern" = "#E69F00")) +
  scale_shape_manual(name = "Hemishpere",
                     values = c("Northern" = 16, "Southern" = 17)) +
  scale_linetype_manual(name = "Legend", values = c("Northern Slope" = "solid",
                                                    "Southern Slope" = "solid"))+
  guides(
    color = guide_legend(override.aes = list(color = c("#D3D3D3", "#0072B2", "#E69F00"), shape = c(15, 16, 17))),
    linetype = "none"
  ) +
  # Facet by bin_midpoint with 8 columns
  scale_y_continuous(
    breaks = function(y) {
      max_val <- ceiling(max(y, na.rm = TRUE) / 10) * 10
      mid_val <- ceiling(max_val / 2 / 10) * 10
      return(c(0, mid_val, max_val))
    }
  )  +
  facet_wrap(~ bin_midpoint,
             labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)])),
             scales = "free_y", ncol = 6) +
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
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  # Restore x-axis labels ONLY on the bottom row & y-axis labels ONLY on the left column
  theme(
    axis.text.x = element_text(size = 8, color = "black"),  # Show x-axis labels at the bottom
    axis.title.x = element_text(size = 12, color = "black"),  # Show x-axis title
    axis.text.y = element_text(size = 8, color = "black"),  # Show y-axis labels on the left
    axis.title.y = element_text(size = 12, color = "black", angle = 90)  # Show y-axis title
  )

print(LDG_s_plot)
gg_path <- sprintf("./figures/%s km LDG slopes figure Mrate is %s and col is %s.jpg", 
                   params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)

ggsave(gg_path, LDG_s_plot, width = 8, height = 9, dpi = 300)
