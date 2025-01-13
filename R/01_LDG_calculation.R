# Wendy geowendywen@outlook.com

data_clean <- as.data.frame(data_clean)
dat_list <- lapply(unique(data_clean$stbin), function(i) subset(data_clean, stbin == i))
# 1.1 Information in stages

# Start timer
tic("Running buffers function")
# Apply the compute_richness_summary function to each element in dat_list
richness_results <- lapply(dat_list, function(dat) {
  compute_richness_summary(
    dat = dat,
    xy = xy,
    iter = iter,
    nSite = nSite,
    r = r,
    crs = crs,
    q = q, 
    datatype = datatype, 
    base = base,
    level = level, 
    nboot = nboot
  )
})

# Extract and combine the detailed richness results (rich_df) across all datasets
combined_rich_df <- do.call(rbind, lapply(richness_results, `[[`, "combined_rich"))

# Extract and combine the summary richness results (rich_info_summary) across all datasets
combined_rich_summary <- do.call(rbind, lapply(richness_results, `[[`, "rich_info_summary"))

# Stop timer and display elapsed time
toc()
# 1.2 LDG using coverage-based rarefaction in buffer area
