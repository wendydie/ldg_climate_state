suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(GGally)
  library(palaeoverse)
})

source("./R/options.R")
source("./R/functions/check_hemisphere_good.R")

set.seed(123)

B <- 100
k_min <- 1
occurrence_min <- 5
out_dir <- "./figures/test"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hemi_cols <- c("N" = "#0072B2", "S" = "#E69F00")
slope_names <- c("direct_norm", "resample_norm", "direct_raw", "resample_raw")

rich_df <- read.csv(sprintf("./results/LDG/%s_cell_%s_richness.csv",
                            params$spacing, params$level))

time_bins <- readRDS("./data/time_bins.RDS")

lat_bins <- palaeoverse::lat_bins_area(n = rich_params$n_lat_bins) %>%
  arrange(min)

rich_df <- rich_df %>%
  filter(nT >= occurrence_min, t <= 2 * nT, bin_midpoint <= 486.8500) %>%
  mutate(
    stage = time_bins$interval_name[match(bin_midpoint, time_bins$mid_ma)],
    bin = lat_bins$bin[findInterval(cell_lat, vec = c(lat_bins$min, Inf))],
    abs_lat = abs(cell_lat),
    hemisphere = ifelse(cell_lat >= 0, "Northern", "Southern")
  ) %>%
  left_join(lat_bins %>% transmute(bin, abs_lat_bin_mid = abs(mid)), by = "bin") %>%
  filter(!is.na(qD), !is.na(abs_lat), !is.na(abs_lat_bin_mid),
         !is.na(bin_midpoint), !is.na(hemisphere)) %>%
  group_by(bin_midpoint) %>%
  mutate(qD_normalized = qD * 100 / max(qD, na.rm = TRUE)) %>%
  ungroup()

if (rich_params$n_lat_bins == 6) {
  south_bin_no <- lat_bins$bin[lat_bins$mid >= -40 & lat_bins$mid <= 0]
  north_bin_no <- lat_bins$bin[lat_bins$mid <= 40 & lat_bins$mid >= 0]
  
  rich_df <- rich_df %>%
    group_by(bin_midpoint, hemisphere) %>%
    mutate(label = ifelse(
      (hemisphere == "Southern" & all(south_bin_no %in% bin)) |
        (hemisphere == "Northern" & all(north_bin_no %in% bin)),
      "good", "bad"
    )) %>%
    ungroup()
} else {
  rich_df <- has_adjacent_bins(rich_df, lat_bins)
}

rich_df <- rich_df %>% filter(label == "good")

fit_slope <- function(df, y_col) {
  df <- df %>% filter(!is.na(.data[[y_col]]), !is.na(abs_lat))
  if (nrow(df) < 3 || n_distinct(df$abs_lat) < 2) return(NA_real_)
  coef(lm(reformulate("abs_lat", y_col), data = df))[["abs_lat"]]
}

calc_slope <- function(df, y_col, slope_name, balanced = FALSE, B = 100, k_min = 1) {
  df %>%
    group_by(bin_midpoint, stage, hemisphere) %>%
    group_modify(~ {
      dat <- .x %>% filter(!is.na(.data[[y_col]]), !is.na(abs_lat), !is.na(abs_lat_bin_mid))
      
      if (!balanced) {
        return(tibble(slope = fit_slope(dat, y_col)))
      }
      
      lat_counts <- dat %>% count(abs_lat_bin_mid, name = "n_cells")
      if (nrow(lat_counts) < 2) return(tibble(slope = NA_real_))
      
      k_sample <- min(lat_counts$n_cells, na.rm = TRUE)
      if (is.na(k_sample) || k_sample < k_min) return(tibble(slope = NA_real_))
      
      split_dat <- split(dat, dat$abs_lat_bin_mid)
      
      tibble(slope = median(map_dbl(seq_len(B), function(i) {
        sampled <- bind_rows(lapply(split_dat, \(z)
                                    z[sample(seq_len(nrow(z)), k_sample), , drop = FALSE]
        ))
        fit_slope(sampled, y_col)
      }), na.rm = TRUE))
    }) %>%
    ungroup() %>%
    mutate(slope_type = slope_name)
}

slope_all <- bind_rows(
  calc_slope(rich_df, "qD_normalized", slope_names[1], balanced = FALSE),
  calc_slope(rich_df, "qD_normalized", slope_names[2], balanced = TRUE, B = B, k_min = k_min),
  calc_slope(rich_df, "qD", slope_names[3], balanced = FALSE),
  calc_slope(rich_df, "qD", slope_names[4], balanced = TRUE, B = B, k_min = k_min)
)

slope_wide <- slope_all %>%
  pivot_wider(
    id_cols = c(bin_midpoint, stage, hemisphere),
    names_from = slope_type,
    values_from = slope
  ) %>%
  mutate(hemisphere_short = recode(hemisphere, Northern = "N", Southern = "S"))

slope_matrix <- slope_wide %>%
  select(hemisphere = hemisphere_short, all_of(slope_names))

correlation_plot <- ggpairs(
  slope_matrix,
  columns = 2:ncol(slope_matrix),
  aes(color = hemisphere, fill = hemisphere),
  upper = list(continuous = wrap("cor", size = 3)),
  lower = list(continuous = wrap("points", alpha = 0.7, size = 2)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5))
) +
  scale_color_manual(values = hemi_cols) +
  scale_fill_manual(values = hemi_cols) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    legend.position = "bottom",
    plot.background = element_blank()
  )

print(correlation_plot)

cor_path <- sprintf(
  "%s/%skm_%squota_%slatbins_per_cell_cor_k%s_B%s.jpg",
  out_dir, params$spacing, params$level, rich_params$n_lat_bins, k_min, B
)

ggsave(cor_path, correlation_plot, width = 6, height = 6, dpi = 300)

print(
  slope_wide %>%
    group_by(hemisphere) %>%
    summarise(
      n = n(),
      across(all_of(slope_names), ~ sum(!is.na(.x)), .names = "n_{.col}"),
      .groups = "drop"
    )
)

print(cor_path)