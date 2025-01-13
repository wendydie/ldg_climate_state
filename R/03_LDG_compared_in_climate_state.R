### Compare LDG with different climate states

# 3.1 calculating the slope of LDG in time bins (N,S,combined)

# 3.2 t-test analysis in different climate states (N,S,combined)
climate_states <- read.csv("./data/climate_states.csv")

# Step 1: Merge final_results and climate_states by stbin
slope_cli_df <- final_results %>%
  left_join(climate_states, by = "stbin")

# Step 2: Group by climate_state and perform t-tests
t_test_results <- slope_cli_df %>%
  group_by(climate_state) %>%
  summarise(
    # Perform t-test between southern_slope and northern_slope
    t_test_result = list(t.test(
      southern_slope, northern_slope,
      alternative = "two.sided",   # Two-tailed test
      paired = FALSE              # Non-paired data
    )),
    # Extract p-value from the t-test
    p_value = t_test_result[[1]]$p.value,
    # Extract t-statistic from the t-test
    t_statistic = t_test_result[[1]]$statistic,
    # Calculate the mean of southern_slope
    mean_southern_slope = mean(southern_slope, na.rm = TRUE),
    # Calculate the mean of northern_slope
    mean_northern_slope = mean(northern_slope, na.rm = TRUE),
    .groups = 'drop' # Remove grouping after summarizing
  )

# 3.3 drawing the slope figures of LDG of time-bins in one figure
# Step 1: Prepare the data
# Ensure final_results contains time-bin (stbin) and slopes for combined, southern, and northern groups
slope_data <- slope_cli_df %>%
  pivot_longer(
    cols = c(combined_slope, southern_slope, northern_slope),
    names_to = "slope_type",
    values_to = "slope_value"
  )

# Step 2: Create the slope plot
slope_plot <- ggplot(slope_data, aes(x = stbin, y = slope_value, color = slope_type)) +
  geom_line(size = 1) +           # Add lines for slopes
  geom_point(size = 2) +          # Add points for better visibility
  scale_color_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  labs(
    title = "Slopes of Latitudinal Diversity Gradients Over Time",
    x = "Time Bins (stbin)",
    y = "Slope Value",
    color = "Slope Type"
  ) +
  theme_minimal() +               # Use a minimal theme for clean visualization
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title
    axis.title = element_text(size = 14),                             # Adjust axis titles
    axis.text = element_text(size = 12),                              # Adjust axis text
    legend.title = element_text(size = 14),                           # Adjust legend title
    legend.text = element_text(size = 12)                             # Adjust legend text
  )

# 3.4 drawing the slope figures of LDG in each time bin
facet_plot <- ggplot(slope_data, aes(x = slope_type, y = slope_value, fill = slope_type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Use bar plots for clarity
  facet_wrap(~ stbin, scales = "free", ncol = 3) +                # Create one facet per time bin
  scale_fill_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  labs(
    title = "Slope Figures of Latitudinal Diversity Gradients by Time Bin",
    x = "Slope Type",
    y = "Slope Value",
    fill = "Slope Type"
  ) +
  theme_minimal() +                                                # Use a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center the title
    axis.title = element_text(size = 14),                            # Adjust axis title size
    axis.text = element_text(size = 12),                             # Adjust axis text size
    legend.title = element_text(size = 14),                          # Adjust legend title
    legend.text = element_text(size = 12),                           # Adjust legend text
    strip.text = element_text(size = 14, face = "bold")              # Adjust facet label size
  )

# 3.5 drawing the 箱型图 of slope of LDG in different climate states
boxplot <- ggplot(slope_data, aes(x = climate_state, y = slope_value, fill = slope_type)) +
  geom_boxplot(outlier.size = 2, outlier.shape = 21) +   # Add boxplots with outliers highlighted
  scale_fill_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  labs(
    title = "Boxplot of LDG Slopes in Different Climate States",
    x = "Climate State",
    y = "Slope Value",
    fill = "Slope Type"
  ) +
  theme_minimal() +                                         # Apply a minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center title
    axis.title = element_text(size = 14),                              # Adjust axis titles
    axis.text = element_text(size = 12),                               # Adjust axis text size
    legend.title = element_text(size = 14),                            # Adjust legend title
    legend.text = element_text(size = 12)                              # Adjust legend text
  )

# 3.6 the changes in climate states need to be shown on the figure, which can be placed at the top. For example, changes can be represented using bars of different colors.
# Prepare climate states for overlay
climate_state_data <- slope_cli_df %>%
  distinct(stbin, climate_state) %>%
  mutate(climate_state_numeric = as.numeric(factor(climate_state)))  # Assign numeric codes for states

# Step 2: Create the plot
slope_climate_plot <- ggplot() +
  # Add the climate state bars (as rectangles)
  geom_rect(data = climate_state_data, aes(
    xmin = stbin - 0.5, xmax = stbin + 0.5,
    ymin = max(slope_data$slope_value) + 0.1, ymax = max(slope_data$slope_value) + 0.3,
    fill = climate_state
  )) +
  # Add the slope lines
  geom_line(data = slope_data, aes(
    x = stbin, y = slope_value, color = slope_type, group = slope_type
  ), size = 1) +
  # Add points for better visualization
  geom_point(data = slope_data, aes(
    x = stbin, y = slope_value, color = slope_type
  ), size = 2) +
  # Customize colors for slope lines
  scale_color_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  # Customize colors for climate states
  scale_fill_manual(
    values = c("Warm" = "orange", "Cold" = "lightblue", "Transition" = "purple"),
    name = "Climate State"
  ) +
  # Adjust axis labels and title
  labs(
    title = "LDG Slopes and Climate State Changes Through Time",
    x = "Time Bins (stbin)",
    y = "Slope Value",
    color = "Slope Type"
  ) +
  # Adjust theme for clarity
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

