library(ggplot2)
library(dplyr)

# 过滤掉 bad 数据，并确保 bin_midpoint 为 numeric
slope_data <- LDG_slope %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint)),
         slope = ifelse(label == "bad", NA, slope))
# 定义颜色映射
quantile_colors <- c(
  "q25" = "#1F77B4",   # 蓝色
  "q50" = "#2CA02C",   # 绿色
  "q75" = "#FF7F0E",   # 橙色
  "q95" = "#D62728"    # 红色
)

# 绘制 slope 变化曲线
slope_vTime_plot <- ggplot(slope_data, aes(x = bin_midpoint, y = slope, color = quantile)) +
  geom_line(linewidth = 1) +  # 绘制曲线
  geom_point(size = 2) +  # 绘制散点
  
  # 添加气候状态颜色条
  geom_rect(data = climate_states, 
            aes(xmin = top, xmax = bottom, ymin = -Inf, ymax = Inf, fill = climate_color),
            alpha = 0.3, inherit.aes = FALSE) +
  
  # Y 轴 & X 轴
  scale_x_reverse(name = "Geological time (Ma)", breaks = seq(500, 0, -50), expand = c(0, 0)) +
  scale_y_continuous(name = "Slope Value", expand = c(0, 0)) +
  
  # 颜色映射到 quantile
  scale_color_manual(name = "Quantile Slopes", values = quantile_colors) +
  
  # 其他修饰
  facet_wrap(~ hemisphere, scales = "free_y", ncol=1) +
  labs(
    x = "Geological Time (Ma)",
    y = "Slope Value"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = "right"
  )

# 显示图像
print(slope_vTime_plot)

library(ggplot2)
library(dplyr)

# 过滤掉 bad 数据，并确保 bin_midpoint 为 numeric
slope_data <- LDG_slope %>%
  filter(label != "bad") %>%
  mutate(bin_midpoint = as.numeric(as.character(bin_midpoint)))

# 合并不同 quantile，并确保 climate_state 顺序
slope_boxplot_data <- slope_data %>%
  left_join(climate_states, by = c("bin_midpoint" = "mid")) %>%
  mutate(
    climate_state = factor(climate_state, levels = c("Coldhouse", "Coolhouse", "Transitional", "Warmhouse", "Hothouse"))
  )

# 定义 quantile 颜色
quantile_colors <- c(
  "q25" = "#1F77B4",   # 蓝色
  "q50" = "#2CA02C",   # 绿色
  "q75" = "#FF7F0E",   # 橙色
  "q95" = "#D62728"    # 红色
)

# 绘制 Boxplot
slope_boxplot <- ggplot(slope_boxplot_data, aes(x = climate_state, y = slope, fill = quantile)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # 绘制箱线图
  geom_jitter(aes(color = quantile), width = 0.2, size = 1, alpha = 0.5) +  # 添加数据散点
  scale_fill_manual(name = "Quantile", values = quantile_colors) +  # 填充颜色
  scale_color_manual(name = "Quantile", values = quantile_colors) +  # 点颜色
  facet_wrap(~ hemisphere, ncol=1) +  # 按照南北半球分面
  labs(
    x = "Climate State",
    y = "Slope Value",
    title = "LDG Slope Distributions Across Climate States"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

# 显示图像
print(slope_boxplot)

sb_path <- sprintf("./figures/test/%s km LDG boxplot figure Mrate is %s and col is %s.jpg", 
                   params$spacing, rich_params$Mrate_filter_no, rich_params$col_filter_no)

ggsave(sb_path,slope_boxplot, width = 6, height = 8, dpi = 300)
