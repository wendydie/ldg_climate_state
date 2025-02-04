###TEST

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
library(spdep)
library(spatialreg)
library(ncf)
library(viridis)
library(sf)
library(patchwork)
source("./R/02_LDG_slope.R")

# rich_df <- read.csv("./results/LDG/250_cell_0.7_richness.csv") #_nonunique_occurrence
# time_bins <- readRDS("./data/time_bins.RDS")
TAS <- readRDS("./data/processed/TAS_Scenario8.RDS")

# Step 1 : data filter ---------------------------------------------------
good_stage <- rich_df %>%
  filter(label == "good") %>%
  distinct(bin_midpoint)

# Step 2 : match simulation temperature-------------------------------------
# Function to calculate mean TAS for each cell_lat, cell_lng, and Bin
calculate_mean_tas <- function(cell_lat, cell_lng, bin_midpoint, radius_deg, tas_data) {
  # Filter TAS data within the 2.5° range and matching Bin
  nearby_points <- tas_data %>%
    filter(
      Latitude >= (cell_lat - radius_deg) & Latitude <= (cell_lat + radius_deg),
      Longitude >= (cell_lng - radius_deg) & Longitude <= (cell_lng + radius_deg),
      bin_midpoint == bin_midpoint
    )
  
  # Compute mean TAS if there are nearby points, otherwise return NA
  if (nrow(nearby_points) > 0) {
    return(mean(nearby_points$TAS, na.rm = TRUE))
  } else {
    return(NA)
  }
}

# Apply function to each row in rich_df
rich_df$matched_TAS <- mapply(calculate_mean_tas, 
                              rich_df$cell_lat, 
                              rich_df$cell_lng, 
                              rich_df$bin_midpoint, 
                              radius_deg=2.5,
                              MoreArgs = list(tas_data = TAS))

# calculate a fully rescaled version
rich_df_frs <- rich_df %>%
  select(bin_midpoint, cell_lat, cell_lng, qD, matched_TAS)
# %>%
#   group_by(bin_midpoint) %>%
#   mutate(across(c(qD, matched_TAS), ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

# rich_df_frs <- rich_df_frs %>%
#   filter(cell_lat < 60)

# Convert to spatial object
rich_sf <- st_as_sf(rich_df_frs, coords = c("cell_lng", "cell_lat"), crs = 4326) %>%
  st_transform(3857)  # Convert to planar coordinates for easier calculations

# Step 3: SAR correlation-------------------------------------------------

analyze_spatial_model <- function(df_subset) {
  # Input: Subset of data for a single time period
  # Output: List of model results with complete statistics
  
  # Handle cases where data points are insufficient -----------------------
  if(nrow(df_subset) < 10) return(NULL)
  
  # Construct spatial weight matrix ---------------------------------------
  coords <- st_coordinates(df_subset)
  k <- min(5, nrow(df_subset)-1)  # Dynamically adjust the number of neighboring points
  knn_nb <- knearneigh(coords, k = k)
  nb <- knn2nb(knn_nb)
  w <- nb2listw(nb, style = "W")  # Row-standardized weight matrix
  
  # Spatial autocorrelation test (Moran's I) ---------------------------
  moran_test <- tryCatch(
    moran.test(df_subset$qD, w),
    error = function(e) {
      message("Moran's I test failed:", e$message)
      list(estimate = c(NA, NA, NA), p.value = NA)
    }
  )
  
  # OLS model construction ---------------------------------------------
  ols <- try(lm(qD ~ matched_TAS, data = df_subset), silent = TRUE)
  if(inherits(ols, "try-error")) {
    message("OLS model fitting failed:", geterrmessage())
    return(NULL)
  }
  
  # Extract OLS model statistics ---------------------------------------
  ols_summary <- summary(ols)
  ols_coef <- coef(ols)["matched_TAS"]
  ols_se <- ols_summary$coefficients["matched_TAS", "Std. Error"]
  ols_p <- ols_summary$coefficients["matched_TAS", "Pr(>|t|)"]
  
  # Test for spatial dependence (Rao's score tests) --------------------
  RS_test <- tryCatch(
    lm.RStests(ols, listw = w, test = c("RSlag", "RSerr")),
    error = function(e) {
      message("Spatial dependence test failed:", e$message)
      NULL
    }
  )
  
  # Initialize model parameters ----------------------------------------
  model_type <- "OLS"
  spatial_coef <- NA
  spatial_se <- NA
  spatial_p <- NA
  
  # Model selection logic ----------------------------------------------
  if(!is.null(RS_test)) {
    # Safely extract p-values (handle potential NULL values)
    p_lag <- tryCatch(RS_test$RSlag$p.value, error = function(e) NA)
    p_err <- tryCatch(RS_test$RSerr$p.value, error = function(e) NA)
    
    # Handle NA values (treat as insignificant)
    p_lag <- ifelse(is.na(p_lag), 1, p_lag)
    p_err <- ifelse(is.na(p_err), 1, p_err)
    
    # Core decision logic
    if(p_lag < 0.05 | p_err < 0.05) {
      # Spatial model selection branch
      if(p_lag < p_err) {
        model <- tryCatch(
          lagsarlm(qD ~ matched_TAS, data = df_subset, listw = w),
          error = function(e) {
            message("SAR model fitting failed:", e$message)
            NULL
          }
        )
        if(!is.null(model)) {
          model_summary <- summary(model)
          spatial_coef <- model_summary$Coef["matched_TAS", "Estimate"]
          spatial_se <- model_summary$Coef["matched_TAS", "Std. Error"]
          spatial_p <- model_summary$Coef["matched_TAS", "Pr(>|z|)"]
          model_type <- "SAR"
        }
      } else {
        model <- tryCatch(
          errorsarlm(qD ~ matched_TAS, data = df_subset, listw = w),
          error = function(e) {
            message("SEM model fitting failed:", e$message)
            NULL
          }
        )
        if(!is.null(model)) {
          model_summary <- summary(model)
          spatial_coef <- model_summary$Coef["matched_TAS", "Estimate"]
          spatial_se <- model_summary$Coef["matched_TAS", "Std. Error"]
          spatial_p <- model_summary$Coef["matched_TAS", "Pr(>|z|)"]
          model_type <- "SEM"
        }
      }
    }
  }
  
  # Confidence interval calculation ------------------------------------
  calculate_ci <- function(estimate, se) {
    if(is.na(estimate) | is.na(se)) return(c(NA, NA))
    c(
      estimate - 1.96 * se,
      estimate + 1.96 * se
    )
  }
  
  # OLS confidence interval
  ols_ci <- calculate_ci(ols_coef, ols_se)
  
  # Spatial model confidence interval
  spatial_ci <- calculate_ci(spatial_coef, spatial_se)
  
  # Significance marker generation -------------------------------------
  get_significance <- function(p) {
    case_when(
      p <= 0.001 ~ "***",
      p <= 0.01 ~ "**",
      p <= 0.05 ~ "*",
      TRUE ~ ""
    )
  }
  
  # Return structured results -----------------------------------------
  list(
    # Basic information
    bin_midpoint = unique(df_subset$bin_midpoint),
    n_obs = nrow(df_subset),
    
    # Spatial autocorrelation test
    moran_I = moran_test$estimate[1],
    moran_p = moran_test$p.value,
    
    # OLS results
    ols = list(
      estimate = ols_coef,
      se = ols_se,
      ci_low = ols_ci[1],
      ci_high = ols_ci[2],
      p_value = ols_p,
      significance = get_significance(ols_p)
    ),
    
    # Spatial model results
    spatial = list(
      type = model_type,
      estimate = spatial_coef,
      se = spatial_se,
      ci_low = spatial_ci[1],
      ci_high = spatial_ci[2],
      p_value = spatial_p,
      significance = get_significance(spatial_p)
    )
  )
}

# Results Processing and Visualization -----------------------------------------------------------
# Process grouped results
results <- rich_sf %>%
  group_by(bin_midpoint) %>%
  group_modify(~ {
    res <- analyze_spatial_model(.x)
    if(is.null(res)) return(data.frame())
    
    data.frame(
      n_obs = res$n_obs,
      
      # Spatial Autocorrelation
      moran_I = res$moran_I,
      moran_p = res$moran_p,
      
      # OLS Results
      ols_estimate = res$ols$estimate,
      ols_se = res$ols$se,
      ols_ci_low = res$ols$ci_low,
      ols_ci_high = res$ols$ci_high,
      ols_p = res$ols$p_value,
      ols_sig = res$ols$significance,
      
      # Spatial Model Results
      spatial_type = res$spatial$type,
      spatial_estimate = res$spatial$estimate,
      spatial_se = res$spatial$se,
      spatial_ci_low = res$spatial$ci_low,
      spatial_ci_high = res$spatial$ci_high,
      spatial_p = res$spatial$p_value,
      spatial_sig = res$spatial$significance
    )
  }) %>%
  ungroup()

results <- results %>%  # Convert to numeric
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint))) %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid"))
results <- results %>%
  mutate(climate_state = as.character(climate_state))
results_filter <- results %>%
  filter(bin_midpoint %in% good_stage$bin_midpoint)

plot_spatial_results <- function(results) {
  # Generate climate period data (new addition)
  climate_periods <- results %>%
    arrange(bin_midpoint) %>%
    select(bin_midpoint, climate_state) %>%
    distinct() %>%
    mutate(
      change = climate_state != lag(climate_state, default = first(climate_state)),
      group = cumsum(change)
    ) %>%
    group_by(group, climate_state) %>%
    summarise(
      start = min(bin_midpoint) - 2.5,  # Assume each period width is 5 million years
      end = max(bin_midpoint) + 2.5,
      .groups = 'drop'
    )
  
  # Data preparation (original code)
  plot_data <- results %>%
    pivot_longer(
      cols = c(ols_estimate, spatial_estimate),
      names_to = "model_type",
      values_to = "estimate",
      names_pattern = "(.*)_estimate"
    ) %>%
    mutate(
      ci_low = case_when(
        model_type == "ols" ~ ols_ci_low,
        model_type == "spatial" ~ spatial_ci_low
      ),
      ci_high = case_when(
        model_type == "ols" ~ ols_ci_high,
        model_type == "spatial" ~ spatial_ci_high
      ),
      significance = case_when(
        model_type == "ols" ~ ols_sig,
        model_type == "spatial" ~ spatial_sig
      ),
      model_type = factor(
        case_when(
          model_type == "ols" ~ "OLS",
          spatial_type == "SAR" ~ "SAR",
          spatial_type == "SEM" ~ "SEM"
        ),
        levels = c("OLS", "SAR", "SEM")
      )
    )
  
  # Create visualization (add background layer)
  ggplot(plot_data, aes(x = bin_midpoint, y = estimate, color = model_type)) +
    # Climate state background layer (new addition)
    geom_rect(
      data = climate_periods,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = climate_state),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    # Original visualization elements
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.3,
      linewidth = 0.8,
      alpha = 0.7
    ) +
    geom_point(
      aes(shape = model_type),
      size = 3,
      fill = "white"
    ) +
    geom_text(
      aes(label = significance, y = ci_high + 0.1*sd(estimate, na.rm = TRUE)),
      size = 4,
      vjust = 0,
      color = "black",
      show.legend = FALSE
    ) +
    # Color and fill scales (new fill scale)
    scale_color_manual(
      values = c("#1f77b4", "#ff7f0e", "#2ca02c"),
      labels = c("Ordinary Least Squares", "Spatial Lag Model", "Spatial Error Model")
    ) +
    scale_fill_manual(
      values = c(Warm = "#FFB6C1", Cold = "#87CEFA", Stable = "#98FB98"),
      name = "Climate State"
    ) +
    scale_shape_manual(
      values = c(16, 17, 15),
      labels = c("Ordinary Least Squares", "Spatial Lag Model", "Spatial Error Model")
    ) +
    # Adjust x and y axis ranges
    scale_x_reverse(
      name = "Geological time (Ma)",
      limits = c(x_max_val, 0),  # Explicitly set limits to ensure alignment
      breaks = seq(500, 0, -50),
      expand = c(0, 0)  # Avoid extra whitespace
    )  +
    # Labels and theme settings
    labs(
      x = "Geological time (Ma)",
      y = "Estimated temperature effect coefficient (95% CI)",
      color = "Model type",
      shape = "Model type",
      title = "Spatiotemporal effects of temperature on marine invertebrate diversity",
      subtitle = "Error bars indicate 95% confidence intervals; * p<0.05, ** p<0.01, *** p<0.001"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      axis.title = element_text(face = "bold"),
      legend.box = "horizontal",
      legend.spacing = unit(0.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    guides(
      color = guide_legend(nrow = 1),
      shape = guide_legend(nrow = 1),
      fill = guide_legend(nrow = 1)
    )
}

# Generate final plot
final_plot <- plot_spatial_results(results_filter)
print(final_plot)
ggsave("./figures/temperature_effect_analysis.png", final_plot, 
       width = 12, height = 7, dpi = 300)

# Ensure climate_state is a character type--------------------------------
results_filter <- results_filter %>%
  mutate(climate_state = as.character(climate_state))

# Remove rows where climate_state is missing
results_filter <- results_filter %>%
  filter(!is.na(climate_state))  

# Check unique climate states
unique(results_filter$climate_state)

# Define comparison groups to match unique(results_filter$climate_state)
comparisons_list <- list(
  c("Hothouse", "Warmhouse"),
  c("Hothouse", "Coolhouse"),
  c("Hothouse", "Transitional"),
  c("Hothouse", "Coldhouse"),
  c("Warmhouse", "Coolhouse"),
  c("Warmhouse", "Transitional"),
  c("Warmhouse", "Coldhouse"),
  c("Coolhouse", "Transitional"),
  c("Coolhouse", "Coldhouse"),
  c("Transitional", "Coldhouse")
)

library(ggpubr)

# Create boxplot and perform statistical tests
boxplot_plot <- ggplot(results_filter, aes(x = climate_state, y = ols_estimate, fill = climate_state)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparency 0.7, remove outliers
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add scatter points to avoid overlapping data
  stat_compare_means(method = "kruskal.test", label.y = max(results_filter$ols_estimate, na.rm = TRUE) * 1.1) +  # Kruskal-Wallis test for group comparison
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_list) +  # Pairwise Wilcoxon test
  scale_fill_manual(values = c("Warm" = "#FFB6C1", "Cold" = "#87CEFA", "Stable" = "#98FB98")) +  # Custom colors
  labs(
    x = "Climate State",
    y = "OLS Estimate",
    fill = "Climate State",
    title = "Distribution of Temperature Effects Across Climate States",
    subtitle = "Boxplots represent data distribution, Wilcoxon test indicates significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.title = element_text(face = "bold")
  )

# Save the plot
ggsave("./figures/SAR_climate_state_boxplot.jpg", boxplot_plot, width = 8, height = 6, dpi = 300)

# Display the plot
print(boxplot_plot)

# Count significance levels within each climate_state
sig_counts <- results_filter %>%
  mutate(ols_sig = ifelse(ols_sig %in% c("*", "**", "***"), "Significant", "Not_Significant")) %>%
  group_by(climate_state, ols_sig) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = ols_sig, values_from = count, values_fill = list(count = 0))  # Convert to wide format, filling missing values

# Display the results
print(sig_counts)

# Calculate percentage of significant results
sig_counts <- sig_counts %>%
  mutate(percentage_significant = Significant / (Significant + Not_Significant))
print(sig_counts)

# Create a bar plot of significant OLS results
ggplot(sig_counts, aes(x = climate_state, y = Significant, fill = climate_state)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Significant OLS Results Across Climate States", y = "Significant OLS Results") +
  theme_minimal()

