#geowendywen@outlook.com

library(dplyr)
library(tidyr)
library(ggplot2)
source("./R/options.R")
# 02 calculating the slope of LDG in time bins (N,S,combined)---------------

# Read dataset--------------------------------------------------------------
rich_df <- read.csv("./results/LDG/250_cell_0.7_richness.csv") #_nonunique_occurrence
time_bins <- readRDS("./data/time_bins.RDS")

# Step 1 : Data filtering --------------------------------------------------
# Filter the dataset to keep only rows where collection_no > 5 and Mrate > 0.3
rich_df$stage <- time_bins$interval_name[match(rich_df$bin_midpoint, time_bins$mid_ma)]
rich_df <- rich_df %>% 
  filter(bin_midpoint < 486.8500)%>%
  mutate(bin_midpoint = factor(bin_midpoint, 
                               levels = sort(unique(bin_midpoint), decreasing = TRUE)))
rich_df <- rich_df %>%
  filter(collection_no >= 5 & Mrate >= 0.3 )

# Convert latitude to absolute values and classify hemispheres
rich_df <- rich_df %>%
  mutate(
    abs_lat = abs(cell_lat),  # Convert latitude to absolute value
    hemisphere = case_when(
      cell_lat >= 0 ~ "Northern",  # good cells & cell_lat >= 0 → Northern
      cell_lat < 0 ~ "Southern"  # good cells & cell_lat < 0 → Southern
    ))
hemisphere_counts <- rich_df %>%
  group_by(bin_midpoint, hemisphere) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = hemisphere, values_from = count, values_fill = 0)

rich_df$Northern <- hemisphere_counts$Northern[match(rich_df$bin_midpoint, hemisphere_counts$bin_midpoint)]
rich_df$Southern <- hemisphere_counts$Southern[match(rich_df$bin_midpoint, hemisphere_counts$bin_midpoint)]

rich_df <- rich_df %>%
  mutate(label = ifelse(Northern  >= rich_params$threshold & Southern >= rich_params$threshold, "good", "bad"))
rm(hemisphere_counts)

# Step 2: Calculate LDG slopes-------------------------------------------
combined_slope <- rich_df %>%
  group_by(bin_midpoint) %>%
  summarise(combined_slope = coef(lm(qD ~ abs_lat))[2], .groups = "drop")

northern_slope <- rich_df %>%
  filter(hemisphere == "Northern") %>%
  group_by(bin_midpoint) %>%
  summarise(northern_slope = coef(lm(qD ~ abs_lat))[2], .groups = "drop")

southern_slope <- rich_df %>%
  filter(hemisphere == "Southern") %>%
  group_by(bin_midpoint) %>%
  summarise(southern_slope = coef(lm(qD ~ abs_lat))[2], .groups = "drop")

# Step 3: Combine the results--------------------------------------------
LDG_slope <- combined_slope
LDG_slope$northern_slope <- northern_slope$northern_slope[match(LDG_slope$bin_midpoint, northern_slope$bin_midpoint)]
LDG_slope$southern_slope <- southern_slope$southern_slope[match(LDG_slope$bin_midpoint, southern_slope$bin_midpoint)]
LDG_slope$label <- rich_df$label[match(LDG_slope$bin_midpoint, rich_df$bin_midpoint)]
rm(combined_slope, northern_slope, southern_slope)
# View the LDG slope
View(LDG_slope)

# Step 4: Create the scatter plot with LDG slopes------------------------
LDG_s_plot <- ggplot(rich_df, aes(x = abs_lat, y = qD, 
                    color = ifelse(label == "bad", "bad", hemisphere),
                    shape = hemisphere,
                    linetype = ifelse(label == "bad", "bad", "normal")
                    )) +
  # Scatter plot points (Northern & Southern hemispheres get automatic colors)
  geom_point(alpha = 0.7, size = 1) +
  # Overlay precomputed LDG slopes with increased line thickness
  geom_smooth(method = "lm", se = FALSE, aes(x = abs_lat, y = qD, linetype = "Combined Slope",
                                             color = ifelse(label == "bad", "bad", "red")),
              linewidth = 1, data = rich_df, inherit.aes = FALSE) +  # Combined slope
  geom_smooth(method = "lm", se = FALSE, aes(x = abs_lat, y = qD, linetype = "Northern Slope", 
                                             color = ifelse(label == "bad", "bad", hemisphere)),
              linewidth = 1, data = filter(rich_df, hemisphere == "Northern"), inherit.aes = FALSE) +  # Northern slope
  geom_smooth(method = "lm", se = FALSE, aes(x = abs_lat, y = qD, linetype = "Southern Slope", 
                                             color = ifelse(label == "bad", "bad", hemisphere)),
              linewidth = 1, data = filter(rich_df, hemisphere == "Southern"), inherit.aes = FALSE) +  # Southern slope
  # Define legend for slope lines
  # Use viridis color palette (same for points and lines)
  scale_color_manual(name = "Legend", 
                     values = c("bad" = "#D3D3D3",
                                "Northern" = "#0072B2",  
                                "Southern" = "#E69F00")) +
  scale_shape_manual(name = "Legend", 
                     values = c("Northern" = 16, "Southern" = 17, "bad" = 16)) +
  scale_linetype_manual(name = "Legend", values = c("Combined Slope" = "solid", 
                                                        "Northern Slope" = "solid", 
                                                        "Southern Slope" = "solid"))+
  guides(
    color = "none",
    linetype = guide_legend(override.aes = list(color = c("black", "#0072B2", "#E69F00"), shape = c(16, 17, 16)))
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
             labeller = as_labeller(function(x) paste0(rich_df$stage[match(x, rich_df$bin_midpoint)], "\n", 
                                                       rich_df$label[match(x, rich_df$bin_midpoint)])), 
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

ggsave("./figures/LDG slopes figure.jpg", LDG_s_plot, width = 8, height = 9, dpi = 300)
