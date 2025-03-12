
# Function to calculateOLS slope for LDG using 15-degree bands
calculate_slope_stats <- function(df, percentile) {
  df %>%
    group_by(bin_midpoint) %>%
    filter(n() > 2, length(unique(lat_band_mid_15)) > 1) %>%  # Use 15-degree bands
    summarise(
      model = list(lm(!!sym(percentile) ~ lat_band_mid_15, data = cur_data())), 
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      slope = coef(model)[2],  # Extract Theil-Sen slope
      intercept = coef(model)[1],  # Extract Theil-Sen intercept
      quantile = percentile  # Store percentile for reference
    ) %>%
    select(-model)
}