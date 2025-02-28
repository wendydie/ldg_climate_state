### Compare LDG with different climate states

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
source("./R/options.R")

# Wilcoxon signed-rank test analysis in different climate states-------------------------------

# Read dataset--------------------------------------------------------------
climate_states <- read.csv("./data/climate_states.csv")
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
# Filter and assign colors
climate_states <- climate_states %>%
  filter(bottom <= 486.8500 & top >= 0) %>%
  mutate(climate_color = climate_colors[climate_state])  # Direct mapping
time_bins <- time_bins %>% filter(min_ma <= 486.8500 & max_ma >= 0)

slope_cli_df <- LDG_slope %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint))) %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid"))

slope_cli_df_filter <- slope_cli_df[(slope_cli_df$climate_state != ''&
                                     slope_cli_df$label != 'bad' &
                                     slope_cli_df$quantile == 'q50' # used the median to do the Wilcoxon signed-rank test
                                     ), ]

# Step 2: Group by climate_state and perform Wilcoxon signed-rank test----------------------

# Get all unique pairs of climate states
climate_pairs <- combn(unique(slope_cli_df_filter$climate_state), 2, simplify = FALSE)

# Whether the LDG slope is normally distributed
test_results <- map_df(climate_pairs, function(pair) {
  group1 <- slope_cli_df_filter %>% filter(climate_state == pair[1]) %>% pull(slope)
  group2 <- slope_cli_df_filter %>% filter(climate_state == pair[2]) %>% pull(slope)
  
  # Perform Shapiro-Wilk normality test
  shapiro1_p <- shapiro.test(group1)$p.value
  shapiro2_p <- shapiro.test(group2)$p.value
  
  # Check if both groups are normally distributed
  normal <- shapiro1_p > 0.05 & shapiro2_p > 0.05
  test_type <- ifelse(normal, "t-test", "Wilcoxon")
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    normal = normal,
    test_type = test_type,
    shapiro_p1 = shapiro1_p,
    shapiro_p2 = shapiro2_p
  )
})

print(test_results)

test_path <- sprintf("./results/%skm test selection of climate states Mrate %s and col is %s.csv", 
                     params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)
write.csv(test_results, test_path, row.names = FALSE)

# Wilcoxon test results
wil_results <- map_df(climate_pairs, function(pair) {
  group1 <- slope_cli_df_filter %>% filter(climate_state == pair[1]) %>% pull(slope)
  group2 <- slope_cli_df_filter %>% filter(climate_state == pair[2]) %>% pull(slope)
  
  wilcox_test <- wilcox.test(group1, group2)
  
  # Compute median difference safely
  median1 <- median(group1, na.rm = TRUE)
  median2 <- median(group2, na.rm = TRUE)
  median_diff <- ifelse(!is.na(median1) & !is.na(median2), median1 - median2, NA)
  
  # Compute additional statistics
  n1 <- sum(!is.na(group1))
  n2 <- sum(!is.na(group2))
  iqr1 <- IQR(group1, na.rm = TRUE)
  iqr2 <- IQR(group2, na.rm = TRUE)
  
  tibble(
    group1 = pair[1],
    group2 = pair[2],
    n1 = n1,
    n2 = n2,
    p_value = wilcox_test$p.value,
    median_diff = median_diff,
    w_statistic = wilcox_test$statistic,
    iqr1 = iqr1,
    iqr2 = iqr2
  )
}) %>% 
  arrange(p_value)

print(wil_results)

wilr_path <- sprintf("./results/%skm wilcoxon test pairs of climate states Mrate %s and col is %s.csv", 
                     params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)
write.csv(wil_results, wilr_path, row.names = FALSE)

# Drawing the slope figures of LDG of time-bins in one figure------------
# Step 1: Prepare the data
# Ensure final_results contains time-bin (stage) and slopes for southern, and northern groups
slope_data <- slope_cli_df %>%
  filter(quantile == 'q50') %>%
  mutate(
    slope_type = case_when(
      hemisphere == "Northern" ~ "Northern",
      hemisphere == "Southern" ~ "Southern",
      TRUE ~ NA_character_  # Handle unexpected cases
    )
  ) %>%
  rename(slope_value = slope)%>%
  mutate(slope_value = ifelse(label == "bad", NA, slope_value))

# Step 2: Create the slope plot
x_max_val <- max(time_bins$min_ma) 
x_min_val <- min(time_bins$max_ma) 

y_min_val <- min(slope_data$slope_value, na.rm = TRUE) * 1.3
y_max_val <- max(slope_data$slope_value, na.rm = TRUE) *1.3
y_text_pos <- y_min_val *0.94  # system text 

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
  y = rep(y_max_val, 5),  # Ensure placement above the main plot
  climate_color = c("#005344", "#007d65", "#c8c7c7", "#b57a51", "#95484b"),  # Color mapping
  facet = "Northern"
)

# Create a separate dataset for "Cooler Climate" and "Warmer Climate" labels
climate_labels <- data.frame(
  facet = "Northern",  # Ensure it only appears in the first facet
  x = c(355, 270),
  y = c(y_max_val * 0.8, y_max_val * 0.8),
  label = c("Cooler climate", "Warmer climate"),
  hjust = c(1,0)
)

slope_vTime_plot <- ggplot(slope_data, aes(x = bin_midpoint, y = slope_value, color = slope_type)) +
  # Add slope curves for different types
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Add climate state color bars (ONLY in the first facet)
  geom_rect(data = climate_legend, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y * 0.8, ymax = y * 0.85, fill = climate_color),
            inherit.aes = FALSE) +
  
  # Add "Cooler Climate" & "Warmer Climate" labels (ONLY in the first facet)
  geom_text(data = climate_labels, 
            aes(x = x, y = y, label = label, hjust=hjust),
            size = 4,  inherit.aes = FALSE) +
  # Add climate state color bars
  geom_rect(data = climate_states, 
            aes(xmin = top, xmax = bottom, ymin = y_max_val*0.9, ymax = y_max_val, fill = climate_color),
            linewidth = 0.3, inherit.aes = FALSE) +
  # Add geological system color bands
  geom_rect(data = system_labels, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y_min_val*0.9, ymax = y_min_val, fill = fill_color),
            color = "black", linewidth = 0.3, inherit.aes = FALSE) +  
  # Add geological stage color bands
  geom_rect(data = time_bins, 
            aes(xmin = min_ma, xmax = max_ma, ymin = y_min_val*0.9, ymax = y_min_val*0.85, fill = stageCol),
            color = "black", linewidth = 0.2, inherit.aes = FALSE) +
  # Force ggplot2 to use the exact colors in the dataset
  scale_fill_identity() +
  # Add a black border, enclosing only the data range
  geom_rect(aes(xmin = x_min_val, xmax = x_max_val, ymin = y_min_val, ymax = y_max_val), 
            color = "black", fill = NA, linewidth = 1) +
  geom_hline(yintercept = 0, color="black", linewidth=0.8) + 
  # Set x and y axis ranges
  scale_x_reverse(
    name = "Geological time (Ma)",
    limits = c(x_max_val, 0),
    breaks = seq(500, 0, -50),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(y_min_val, y_max_val),
    expand = c(0, 0)
  ) +

  # Add geological system names above the color bands
  geom_text(
    data = system_labels, 
    aes(x = mid_ma, y = y_text_pos, label = system),
    angle = 0, 
    vjust = 0.5,      
    hjust = 0.5,      
    size = 4, inherit.aes = FALSE
  ) +
  # Set custom colors for slope types
  scale_color_manual(
    values = c("Southern" = "#E69F00", 
               "Northern" = "#0072B2"),
    labels = c("Southern", "Northern")
  ) +
  # **Facet the plot by slope type (separate subplots for Northern and Southern slopes)**
  facet_wrap(~ slope_type, scales = "free_y", ncol=1) +
  
  # Add custom facet labels inside the plot (Top-right position)
  geom_text(data = slope_data %>%
              group_by(slope_type) %>%
              summarise(x_label = 20,  # Position at the rightmost point
                        y_label = y_max_val*0.85),  # Position above max slope
            aes(x = x_label, y = y_label, label = slope_type),
            hjust = 1, vjust = 1, size = 4, fontface = "bold", inherit.aes = FALSE) +
  # Other aesthetic improvements
  labs(
    x = "Geological Time (Ma)",
    y = "Slope Value"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_blank(),  # Remove facet labels
    strip.background = element_blank(),  # Remove facet label background
    panel.spacing = unit(0.3, "lines"),  # Reduce spacing between facets
    legend.position = "none"  # Remove redundant legend
  )

# Display the plot
print(slope_vTime_plot)

sT_path <- sprintf("./figures/%skm time series Mrate is %s col is %s.jpg", 
                  params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)
ggsave(sT_path, slope_vTime_plot, width = 8, height = 7, dpi = 300)

# Drawing the boxplot of slope of LDG in different climate states--------
# Filter LDG slope
slope_data_filtered <- slope_data %>%
  filter(climate_state %in% climate_levels) %>%
  mutate(climate_state = factor(climate_state, levels = climate_levels))

# Compute the total sample count for each climate state (combining all slope types)
sample_counts <- slope_cli_df_filter %>%
  select(climate_state, bin_midpoint) %>%
  group_by(climate_state) %>%  # Group only by climate state
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", count))  # Format label for display

# Draw boxplot
boxplot <- ggplot(slope_data_filtered, aes(x = climate_state, y = slope_value, fill = slope_type)) +
  # Boxplot layer
  geom_boxplot(outlier.size = 2, outlier.shape = 21, position = position_dodge(width = 0.75)) + 
  # Jittered points layer - ensuring they are placed within the correct box
  geom_jitter(aes(fill = slope_type),  # Ensure jitter uses correct fill colors
              color = "gray",  # Use black border for points to differentiate them
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1, alpha = 0.6) + 
  # Add total sample count labels (centered above the violins)
  geom_text(
    data = sample_counts,
    aes(x = climate_state, y = y_min_val, label = label),  # Use a fixed Y position
    inherit.aes = FALSE,  # Do not inherit fill mapping from Slope_Type
    size = 4,
    color = "black",
    vjust = -0.5  # Slight vertical adjustment
  ) +
  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) + 
  
  # Customizing colors for boxplot and jitter points
  scale_fill_manual(
    values = c("Northern" = "#0072B2",  # Blue for Northern slope
               "Southern" = "#E69F00"),  # Orange for Southern slope
    labels = c("Northern", "Southern")
  ) +
  
  # Labels and titles
  labs(
    title = "Boxplot of LDG Slopes in Different Climate States",
    x = "Climate state",
    y = "Slope value",
    fill = "Slope type"
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

# Save and display
print(boxplot)
bT_path <- sprintf("./figures/%skm boxplot Mrate %s and col is %s.jpg", 
                   params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)
ggsave(bT_path, boxplot, width = 8, height = 5, dpi = 300)


