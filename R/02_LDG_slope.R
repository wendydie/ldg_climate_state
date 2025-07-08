# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 02_LDG_slope.R
# Last updated: 2025-01-21
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(palaeoverse)
source("./R/options.R")
source("./R/functions/calculate_LDG_slope.R")
# Read dataset-----------------------------------------------------------
rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv", 
                            params$spacing, params$level))
time_bins <- readRDS("./data/time_bins.RDS")
lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins)
lat_bins <- lat_bins %>%
  arrange(min)
# Step 1 : Data filtering -----------------------------------------------
# Filter the incomplete cells.
rich_df <- rich_df %>%
  filter (nT>=5 & t <= 2*nT)
# Our time period starts from 486.8500 because the climate state data begins at 486.85.
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>% 
  filter(bin_midpoint < 486.8500)%>%
  mutate(bin_midpoint = factor(bin_midpoint, 
                               levels = sort(unique(bin_midpoint), decreasing = TRUE)))
# Create 12 equal-area latitude bands for LDG slope calculation
rich_df <- rich_df %>%
  mutate(bin = findInterval(cell_lat, vec = lat_bins$min))
rich_df <- rich_df %>%
  left_join(lat_bins %>% select(bin, lat_bin_mid = mid), by = "bin")
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
# Divide the cells into tropical-temperate-pole latitude bands for discarding reliable hemispheres.
expected_lat_bands <- seq(30 / 2, 90, 30)  # Midpoints
rich_df <- rich_df %>%
  mutate(lat_band_mid = (floor(abs_lat / 30) *
                           30) + (30 / 2))
# Classify data as "good" or "bad" using 30-degree bands
rich_df <- rich_df %>%
  group_by(bin_midpoint, hemisphere) %>%
  mutate(
    label = case_when(
      all(c(15, 45) %in% lat_band_mid) ~ "good",  # 15 (0-30 midpoint), 45 (30-60 midpoint)
      TRUE ~ "bad"
    )
  ) %>%
  ungroup()
rich_df <- rich_df %>% 
  mutate(color = ifelse(label == "bad", "Bad hemisphere", hemisphere))  # Add color column
# Create 15-degree latitude bands for LDG slope calculation
# rich_df <- rich_df %>%
#   mutate(lat_band_mid_15 = (floor(abs_lat / rich_params$lat_band_width) * 
#                               rich_params$lat_band_width) + rich_params$lat_band_width / 2)  # Midpoints of 15-degree bands
rich_df <- rich_df %>%
  mutate(abs_lat_bin_mid = abs(lat_bin_mid))  # Midpoints of 15-degree bands

# Step 2 : Calculate richness -------------------------------------------
# Calculate the richness at different percentiles, including q50 (median)
rich_df_summary <- rich_df %>%
  group_by(bin_midpoint, hemisphere, abs_lat_bin_mid) %>%
  summarise(
    q50 = quantile(qD_normalized, 0.50, na.rm = TRUE),
    q60 = quantile(qD_normalized, 0.60, na.rm = TRUE),
    q75 = quantile(qD_normalized, 0.75, na.rm = TRUE),
    q90 = quantile(qD_normalized, 0.90, na.rm = TRUE),
    q95 = quantile(qD_normalized, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Initialize an empty list to store results
slope_results <- list()

# Loop through each percentile and compute Theil-Sen slopes separately
for (perc in rich_params$percentiles) {
  # Compute slopes for the Northern Hemisphere
  northern_slope <- calculate_slope_stats(filter(rich_df_summary, hemisphere == "Northern"), perc) %>%
    mutate(hemisphere = "Northern", quantile = perc)
  # Compute slopes for the Southern Hemisphere
  southern_slope <- calculate_slope_stats(filter(rich_df_summary, hemisphere == "Southern"), perc) %>%
    mutate(hemisphere = "Southern", quantile = perc)
  # Store the results
  slope_results[[perc]] <- bind_rows(northern_slope, southern_slope)
}

# Combine all percentile results into one `LDG_slope` dataframe
LDG_slope <- bind_rows(slope_results)

# Assign "good" or "bad" labels from `rich_df`
LDG_slope$label <- rich_df$label[match(
  paste(LDG_slope$bin_midpoint, LDG_slope$hemisphere),
  paste(rich_df$bin_midpoint, rich_df$hemisphere)
)]
LDG_slope$color <- rich_df$color[match(
  paste(LDG_slope$bin_midpoint, LDG_slope$hemisphere),
  paste(rich_df$bin_midpoint, rich_df$hemisphere)
)]
# View the updated LDG_slope with slope, intercept, and classification
View(LDG_slope)
# Save results to CSV
Ls_path <- sprintf("./results/%skm %squota %s equal-area latitude bins LDG slope.csv", 
                   params$spacing, params$level, rich_params$n_lat_bins)
write.csv(LDG_slope, Ls_path, row.names = FALSE)

# ---------------------------------------------------------------
# Filter for the 50th percentile (q75), remove NA slopes, and keep only "good" data
LDG_slope_q75 <- LDG_slope %>%
  filter(quantile=='q75' & !is.na(slope))  # Remove NA values

# Summarize separately for Northern and Southern hemispheres
slope_summary <- LDG_slope_q75 %>%
  group_by(hemisphere) %>%
  summarise(
    good_count = sum(label == "good"),  # Count "good" labels
    bad_count = sum(label == "bad"),  # Count "bad" labels
    total = n()
  ) %>%
  # Only consider "good" data for slope calculations
  left_join(
    LDG_slope_q75 %>%
      filter(label == "good") %>%
      group_by(hemisphere) %>%
      summarise(
        greater_than_0 = sum(slope > 0),
        less_than_0 = sum(slope < 0),
        equal_to_0 = sum(slope == 0),
        na_count = sum(is.na(slope)),
        .groups = "drop"
      ),
    by = "hemisphere"
  ) %>%
  mutate(
    reverse_LDG = (greater_than_0 / good_count) * 100,
    normal_LDG = (less_than_0 / good_count) * 100,
    flat_LDG = (equal_to_0 / good_count) * 100,
    na_LDG = (na_count / good_count) * 100,
    good_stage = (good_count / total) * 100,
    bad_stage = (bad_count / total) * 100
  ) %>%
  ungroup()

ss_path <- sprintf("./results/%skm %squota %s equal-area latitude bins LDG slope summary.csv", 
                   params$spacing, params$level, rich_params$n_lat_bins)
write.csv(slope_summary, ss_path, row.names = FALSE)
