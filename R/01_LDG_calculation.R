# Wendy geowendywen@outlook.com

data_clean <- as.data.frame(data_clean)
dat_list <- lapply(unique(data_clean$stbin), function(i) subset(data_clean, stbin == i))
dat_list <- Filter(function(x) length(x) >= nSite, dat_list)

# 1.1 Information in stages
cl <- parallel::makeCluster(parallel::detectCores() - 2) # Use all but one core
parallel::clusterEvalQ(cl, {
  source('./R/functions/buffer_subsampling.R')
  source('./R/functions/calculate_Info.R')
  
  # Load necessary libraries
  library(divvy)
  library(tictoc)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(iNEXT)
  library(pbapply)     # Ensure dependencies available in workers
})
parallel::clusterExport(cl, varlist = c("compute_richness_summary", "xy", "nSite", 
                                        "r", "crs", "q", "datatype", "base", "level", "nboot"))
# Start timer
tic("Running buffers function")

# Parallelized computation of richness results for each dataset in dat_list
richness_results <- pblapply(dat_list, function(dat) {
  # Step 1: Generate spatial buffers (locations)
  locs <- buffers(dat, xy, nSite = nSite, r = r, crs = crs, output = 'locs')
  stbin <- unique(dat$stbin)
  locs['stbin'] <- stbin
  
  # Handle cases where buffers return NULL
  if (is.null(locs)) {
    warning("Skipping dataset: No viable buffers created for the given parameters.")
    return(list(locs = NULL, richness_summary = NULL))
  }
  
  # Step 2: Compute richness metrics
  richness_summary <- compute_richness_summary(
    dat = dat,
    xy = xy,
    nSite = nSite,
    r = r,
    crs = crs,
    q = q, 
    datatype = datatype, 
    base = base,
    level = level, 
    nboot = nboot
  )
  
  # Return both `locs` and `richness_summary`
  list(locs = locs, richness_summary = richness_summary)
  # tryCatch({
  #   
  # }, error = function(e) {
  #   # Log and skip the dataset if an error occurs
  #   warning(sprintf("Error processing dataset: %s", e$message))
  #   return(list(locs = NULL, richness_summary = NULL))
  # })
}, cl = cl)


# Stop the cluster
parallel::stopCluster(cl)
# Step 1: Filter out NULL elements from richness_results
filtered_richness_results <- Filter(function(x) !is.null(x$richness_summary), richness_results)

# Step 2: Extract and combine `locs` results
combined_locs <- do.call(rbind, lapply(filtered_richness_results, `[[`, "locs"))

# Step 3: Extract and combine the detailed richness results (`combined_rich`) across all datasets
combined_rich_df <- do.call(rbind, lapply(filtered_richness_results, function(x) x$richness_summary$combined_rich))

# Step 4: Extract and combine the summary richness results (`rich_info_summary`) across all datasets
combined_rich_info_summary <- do.call(rbind, lapply(filtered_richness_results, function(x) x$richness_summary$rich_info_df))

# Optional: Save combined results to files
saveRDS(combined_locs, file = paste0(r, "km2_", "combined_locs.rds"))
saveRDS(combined_rich_df, file = paste0(r, "km2_",level,  "sqs quota combined_rich_df.rds"))
saveRDS(combined_rich_info_summary, file = paste0(r, "km2_", level, "sqs quota combined_rich_info_summary.rds"))

# Stop timer and display elapsed time
toc()
# 1.2 LDG using coverage-based rarefaction in buffer area
