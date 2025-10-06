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
source("./R/functions/check_hemisphere_good.R")
# Read dataset-----------------------------------------------------------
rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv", 
                            params$spacing, params$level))
time_bins <- readRDS("./data/time_bins.RDS")
lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins) %>% 
  arrange(min)
# Step 1 : Data filtering -----------------------------------------------
# Filter the incomplete cells.
rich_df <- rich_df %>%
  filter (nT>=5 & t <= 2*nT)
# Our time period starts from 486.8500 because the climate state data begins at 486.85.
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>% 
  filter(bin_midpoint <= 486.8500) %>%
  mutate(bin_index = findInterval(cell_lat, vec = c(lat_bins$min, Inf)),
         bin = lat_bins$bin[bin_index],
         abs_lat = abs(cell_lat),  # Convert latitude to absolute value
         hemisphere = case_when(
           cell_lat >= 0 ~ "Northern",  # good cells & cell_lat >= 0 → Northern
           cell_lat < 0 ~ "Southern"  # good cells & cell_lat < 0 → Southern
           ),
         lat_band_mid = (floor(abs_lat / 30) *
                           30) + (30 / 2) # Divide the cells into tropical(0,30)-temperate(30,60)-pole(60,90) latitude bands for discarding reliable hemispheres.
         ) %>%
  left_join(lat_bins %>% select(bin, lat_bin_mid = mid), by = "bin") %>%
  mutate(abs_lat_bin_mid = abs(lat_bin_mid))  # Midpoints of 12 equal-area latitude bands

# The richness values are normalized to 100 in each stage by dividing the richness in each cell by the maximum richness within that stage.
rich_df <- rich_df %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized = qD*100 / max(qD)) %>%
  ungroup()
# Classify data as "good" or "bad" -----------------------------------
if (rich_params$n_lat_bins == 6){
  rich_df <- rich_df %>%
    group_by(bin_midpoint, hemisphere) %>%
    mutate(
      label = if (
        # n_distinct(abs_lat_bin_mid) >= 3 &
        sum(lat_band_mid == 15) >= 1 & sum(lat_band_mid == 45) >= 1 # & sum(lat_band_mid == 75) >= 1
      ) "good" else "bad"
    ) %>%
    ungroup()
  south_bin_no <- lat_bins$bin[lat_bins$mid >= -40 & lat_bins$mid <= 0]
  north_bin_no <- lat_bins$bin[lat_bins$mid <= 40 & lat_bins$mid >= 0]
  rich_df <- rich_df %>%
    group_by(bin_midpoint, hemisphere) %>%
    mutate(
      label = ifelse(
        (hemisphere == "Southern" & all(south_bin_no %in% bin)) |
          (hemisphere == "Northern" & all(north_bin_no %in% bin)),
        "good", "bad"
      )
    ) %>%
    ungroup()
} else{
  rich_df <- has_adjacent_bins (rich_df, lat_bins)
}

rich_df <- rich_df %>% 
  mutate(color = ifelse(label == "bad", "Bad hemisphere", hemisphere))  # Add color column

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
all_bin__hemi_perc <- expand.grid(
  bin_midpoint = unique(rich_df$bin_midpoint),
  hemisphere = c("Northern", "Southern"),
  quantile = rich_params$percentiles
)
LDG_slope <- LDG_slope %>%
  right_join(all_bin__hemi_perc, by = c("bin_midpoint", "hemisphere", "quantile"))

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
