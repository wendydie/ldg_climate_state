
# Function to calculate OLS slope for LDG 
calculate_slope_stats <- function(df, percentile) {
  df %>%
    group_by(bin_midpoint) %>%
    filter(n() > 2, length(unique(abs_lat_bin_mid)) > 1) %>%  
    summarise(
      model = list(lm(!!sym(percentile) ~ abs_lat_bin_mid, data = cur_data())), 
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