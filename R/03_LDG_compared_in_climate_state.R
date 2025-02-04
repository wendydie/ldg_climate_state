### Compare LDG with different climate states

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)

# 3.2 t-test analysis in different climate states-------------------------
climate_states <- read.csv("./data/climate_states.csv")
time_bins <- readRDS("./data/time_bins.RDS")
# Order climate states and climate colors
climate_levels <- c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")
climate_colors <- c(
  "Coldhouse" = "#005344",
  "Coolhouse" = "#007d65",
  "Transitional" = "#c8c7c7",
  "Warmhouse" = "#b57a51",
  "Hothouse" = "#95484b"
)

# Step 1: Merge final_results and climate_states by bin_midpoint----------
climate_states <- climate_states %>% filter(bottom <= 486.8500 & top >= 0)
climate_states <- climate_states %>%
  mutate(
    climate_color = as.character(recode(climate_state, !!!climate_colors))  # 添加颜色列
  )
time_bins <- time_bins %>% filter(min_ma <= 486.8500 & max_ma >= 0)
slope_cli_df <- LDG_slope %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint))) %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid"))

slope_cli_df_filter <- slope_cli_df[(slope_cli_df$climate_state != ''&
                                slope_cli_df$label != 'bad'), ]

# Calculate the mean of LDG slope in climate states
slope_mean_df <- slope_cli_df_filter %>%
  group_by(climate_state) %>%
  summarise(
    mean_combined_slope = mean(combined_slope, na.rm = TRUE),
    mean_northern_slope = mean(northern_slope, na.rm = TRUE),
    mean_southern_slope = mean(southern_slope, na.rm = TRUE)
  )
print(slope_mean_df)
write.csv(slope_mean_df, "./results/LDG slopes in climate states.csv", row.names = FALSE)

# Step 2: Group by climate_state and perform t-tests----------------------
# group by climate_state
slope_groups <- slope_cli_df_filter %>%
  select(climate_state, combined_slope) %>%
  drop_na() %>%  # remove NA
  group_by(climate_state) %>%
  group_split()

# Get all unique pairs of climate states
climate_pairs <- combn(unique(slope_cli_df$climate_state), 2, simplify = FALSE)

# t tests
t_results <- map_df(climate_pairs, function(pair) {
  group1 <- slope_cli_df %>% filter(climate_state == pair[1]) %>% pull(combined_slope)
  group2 <- slope_cli_df %>% filter(climate_state == pair[2]) %>% pull(combined_slope)
  
  t_test <- t.test(group1, group2) 
  
  # Compute mean difference safely
  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)
  mean_diff <- ifelse(!is.na(mean1) & !is.na(mean2), mean1 - mean2, NA)
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    p_value = t_test$p.value,
    mean_diff = mean_diff,  
    t_statistic = t_test$statistic
  )
})
print(t_results)

write.csv(t_results, "./results/LDG slope t test pairs of climate states.csv", row.names = FALSE)

# ANOVA test
anova_result <- aov(combined_slope ~ climate_state, data = slope_cli_df)
summary(anova_result)
TukeyHSD(anova_result)

# 3.3 drawing the slope figures of LDG of time-bins in one figure-----------
# Step 1: Prepare the data
# Ensure final_results contains time-bin (bin) and slopes for combined, southern, and northern groups
slope_data <- slope_cli_df %>%
  pivot_longer(
    cols = c(combined_slope, southern_slope, northern_slope),
    names_to = "slope_type",
    values_to = "slope_value"
  )%>%
  mutate(slope_value = ifelse(label == "bad", NA, slope_value))

# Step 2: Create the slope plot
x_max_val <- max(time_bins$min_ma) 
x_min_val <- min(time_bins$max_ma) 

y_min_val <- min(slope_data$slope_value, na.rm = TRUE) - 0.2
y_max_val <- max(slope_data$slope_value, na.rm = TRUE) + 0.2
y_text_pos <- y_min_val - 0.3  # system text 

system_labels <- time_bins %>%
  group_by(system = sys) %>%
  summarise(
    min_ma = min(min_ma),
    max_ma = max(max_ma),
    mid_ma = (min(min_ma) + max(max_ma)) / 2,
    fill_color = first(systemCol)
  ) %>%
  distinct(system, .keep_all = TRUE)

# Create a bar chart dataframe for climate state colors
climate_legend <- data.frame(
  climate_state = factor(c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"),
                         levels = c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse")),
  max_ma = seq(350,290, length.out = 5),  # X-axis values
  min_ma = seq(335,275, length.out = 5),
  y = rep(y_max_val - 0.3, 5),  # Ensure placement above the main plot
  climate_color = c("#005344", "#007d65", "#c8c7c7", "#b57a51", "#95484b")  # Color mapping
)

# Create the slope variation curve plot
slope_vTime_plot <- ggplot() +
  # climate state color bars (Legend bar)
  geom_rect(data = climate_legend, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y-0.1, ymax=y+0.1),
            fill = climate_legend$climate_color) +  # Control bar size
  # Add "Cooler Climate" & "Warmer Climate" labels
  geom_text(aes(x = 355, y = y_max_val - 0.3, label = "Cooler climate"), 
            size = 3, hjust = 1) +
  geom_text(aes(x = 270, y = y_max_val - 0.3, label = "Warmer climate"), 
            size = 3, hjust = 0) +
  geom_rect(data = climate_states,
            aes(xmin = top, xmax = bottom, 
                ymin = y_max_val, ymax = y_max_val + 0.2),
            fill = climate_states$climate_color,
            linewidth = 0.3) +
  # Add geological system color bands
  geom_rect(data = system_labels,
            aes(xmin = min_ma, xmax = max_ma, 
                ymin = y_min_val - 0.4, ymax = y_min_val - 0.2),
            fill = system_labels$fill_color,
            color = "black",
            linewidth = 0.3) +  
  # Add geological stage color bands
  geom_rect(data = time_bins,
            aes(xmin = min_ma, xmax = max_ma, 
                ymin = y_min_val - 0.2, ymax = y_min_val),
            fill = time_bins$stageCol,
            color = "black",
            linewidth = 0.2) +
  # Remove fill color legend (to avoid redundancy)
  guides(fill = "none") +  
  # Add slope curves (different types)
  geom_line(data = slope_data, aes(x = bin_midpoint, y = slope_value, color = slope_type), linewidth = 1) +
  geom_point(data = slope_data, aes(x = bin_midpoint, y = slope_value, color = slope_type), size = 2) +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val+0.2), 
            color = "black", fill = NA, linewidth = 1) +
  # Customize colors
  scale_color_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined slope", "Southern slope", "Northern slope")
  ) +
  # Adjust x and y axis ranges
  scale_x_reverse(
    name = "Geological time (Ma)",
    limits = c(x_max_val, 0),  # Explicitly set limits to ensure alignment
    breaks = seq(500, 0, -50),
    expand = c(0, 0)  # Avoid extra whitespace
  ) +
  scale_y_continuous(limits = c(y_min_val-0.4, y_max_val+0.2),
                     breaks = c(-2, -1, 0, 1),        # Specify tick positions
                     expand = c(0, 0)                 # Disable axis expansion
  ) +
  # Add geological system names above the color bands
  geom_text(
    data = system_labels, 
    aes(x = mid_ma, y = y_text_pos, label = system),  # Use system names and computed midpoints
    angle = 0, 
    vjust = 0.5,      # Vertically centered
    hjust = 0.5,      
    size = 3
  ) +
  # Legends and labels
  labs(
    title = "Slopes of LDGs over geological time",
    x = "Geological time (Ma)",
    y = "Slope value",
    color = "Slope Type"
  ) +
  # Improve aesthetics using a clean theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "top",
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5)
  ) 

# Save the plot as an image file
ggsave("./figures/LDG slopes through geological time.jpg", slope_vTime_plot, width = 8, height = 5, dpi = 300)
# Display the plot
print(slope_vTime_plot)


# 3.5 drawing the boxplot of slope of LDG in different climate states
# Filter LDG slope
slope_data_filtered <- slope_data %>%
  filter(climate_state %in% climate_levels) %>%
  mutate(climate_state = factor(climate_state, levels = climate_levels))  # 设定因子顺序

# Draw boxplot
boxplot <- ggplot(slope_data_filtered, aes(x = climate_state, y = slope_value, fill = slope_type)) +
  # Boxplot layer
  geom_boxplot(outlier.size = 2, outlier.shape = 21, position = position_dodge(width = 0.75)) + 
  # Jittered points layer - ensuring they are placed within the correct box
  geom_jitter(color="grey",  
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1, alpha = 0.6) +  
  geom_hline(yintercept = 0, color="black", linewidth=0.8) + 
  # Customizing colors for both boxplot and points
  scale_fill_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  scale_color_manual(
    values = c("combined_slope" = "blue", 
               "southern_slope" = "red", 
               "northern_slope" = "green"),
    labels = c("Combined Slope", "Southern Slope", "Northern Slope")
  ) +
  # Labels and titles
  labs(
    title = "Boxplot of LDG Slopes in Different Climate States",
    x = "Climate State",
    y = "Slope Value",
    fill = "Slope Type",
    color = "Slope Type"
  ) +
  # Minimal theme with some customizations
  theme_minimal() +                                       
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14),                             
    axis.text = element_text(size = 12),                              
    legend.title = element_text(size = 14),                         
    legend.text = element_text(size = 12),                         
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
ggsave("./figures/LDG slope in climate states boxplot.jpg", boxplot, width = 8, height = 5, dpi = 300)
print(boxplot)

### 3.4 vilion figure--------------------------------------------------------
# Convert data to long format and remove missing values (NA)
slope_melt <- melt(slope_cli_df_filter, 
                   id.vars = "climate_state", 
                   measure.vars = c("combined_slope", "southern_slope", "northern_slope"),
                   variable.name = "Slope_Type", 
                   value.name = "Slope_Value") %>%
  na.omit()  # Remove missing values

# Reorder the factor levels of climate_state
slope_melt$climate_state <- factor(slope_melt$climate_state, levels = climate_levels)

# Compute the total sample count for each climate state (combining all slope types)
sample_counts <- slope_cli_df_filter %>%
  select(climate_state, bin_midpoint) %>%
  group_by(climate_state) %>%  # Group only by climate state
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", count))  # Format label for display

# Adjust Y-axis position for text labels (using a fixed offset)
label_y_pos <- max(slope_melt$Slope_Value, na.rm = TRUE) + 0.3

# Create the violin plot with adjustments
violin_plot <- ggplot(slope_melt, aes(x = climate_state, y = Slope_Value, fill = Slope_Type)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot without trimming
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # Overlay boxplot with no outliers
  # Add total sample count labels (centered above the violins)
  geom_text(
    data = sample_counts,
    aes(x = climate_state, y = label_y_pos, label = label),  # Use a fixed Y position
    inherit.aes = FALSE,  # Do not inherit fill mapping from Slope_Type
    size = 4,
    color = "black",
    vjust = -0.5  # Slight vertical adjustment
  ) +
  # Apply minimal theme
  theme_minimal() +
  labs(
    title = "Distribution of LDG slopes by climate state",
    x = "Climate state",
    y = "Slope calue",
    fill = "Slope type"
  ) +
  # Customize text and border styles
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14),                             
    axis.text = element_text(size = 12),                              
    legend.title = element_text(size = 14),                         
    legend.text = element_text(size = 12),                         
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

# Save the violin plot as an image
ggsave("./figures/LDG slope in climate states violin plot.jpg", violin_plot, width = 8, height = 5, dpi = 300)

# Display the plot
print(violin_plot)
