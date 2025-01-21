
dat <- dat_list[[90]]
incefreq2 <- buffers(dat, xy, nSite=1, r=150,
                     crs = 'epsg:4326', output = 'incidence_freq')

locs_inbuffer <- buffers(dat, xy, nSite = nSite, r = r, crs = crs, output = 'locs')
# Compute richness for each buffer using iNEXT::estimateD
# lapply is used to iterate over each buffer in incefreq2
richness_list <- lapply(incefreq2, function(buffer) {
  iNEXT::estimateD(x = buffer, q = c(0), datatype = "incidence_freq", 
                   base = "coverage", level = 0.7, nboot = 50)
})

# Combine the richness results from all iterations into a single data frame
# Iteration number is added as a column for tracking purposes
combined_richness <- do.call(rbind, lapply(seq_along(richness_list), function(i) {
  richness_list[[i]] %>%
    mutate(Iteration = i)  # Add iteration ID
}))

# Calculate summary statistics (mean and 95% confidence intervals) for each Assemblage
richness_summary <- combined_richness %>%
  group_by(Assemblage) %>%
  summarise(
    # Mean and confidence interval for `t`
    mean_t = mean(t),
    ci_t_lower = mean(t) - 1.96 * sd(t) / sqrt(n()),
    ci_t_upper = mean(t) + 1.96 * sd(t) / sqrt(n()),
    
    # Mean and confidence interval for `SC`
    mean_SC = mean(SC),
    ci_SC_lower = mean(SC) - 1.96 * sd(SC) / sqrt(n()),
    ci_SC_upper = mean(SC) + 1.96 * sd(SC) / sqrt(n()),
    
    # Mean and confidence interval for `qD`
    mean_qD = mean(qD),
    ci_qD_lower = mean(qD) - 1.96 * sd(qD) / sqrt(n()),
    ci_qD_upper = mean(qD) + 1.96 * sd(qD) / sqrt(n())
  )

# Apply Fun.ince to all buffers in all iterations of incefreq2
rich.info.all <- lapply(incefreq2, function(iteration) {
  # Apply Fun.ince to each buffer in the current iteration
  lapply(names(iteration), function(assembly_name) {
    result <- Fun.ince(iteration[[assembly_name]]) # Apply Fun.ince
    # Add Assemblage name as the first column
    c(Assemblage = assembly_name, result)
  })
})

# Combine results across all iterations into a single data frame
rich.infodf.all <- do.call(rbind, lapply(seq_along(rich.info.all), function(i) {
  iteration_results <- do.call(rbind, rich.info.all[[i]])
  # Add Iteration ID
  iteration_results <- cbind(Iteration = i, iteration_results)
  return(iteration_results)
}))

rich.infodf.all <- as.data.frame(rich.infodf.all)
# Set column names for the combined data frame
colnames(rich.infodf.all) <- c("Iteration", "Assemblage", "nT", "Number of occurrence",
                               "Sobs", "Chat", "Good.u", "multition.rate",
                               paste0("Q", 1:10))
# Columns to convert to numeric
numeric_cols <- c("nT", "Number of occurrence", "Sobs", "Chat", "Good.u", 
                  "multition.rate", paste0("Q", 1:10))

# Convert specified columns to numeric
rich.infodf.all[numeric_cols] <- lapply(rich.infodf.all[numeric_cols], as.numeric)
# Summarize results by buffer (Assemblage) across iterations
rich_info_summary <- rich.infodf.all %>%
  group_by(Assemblage) %>% 
  summarise(
    # Mean and confidence intervals for nT
    mean_nT = mean(nT),
    ci_nT_lower = mean(nT) - 1.96 * sd(nT) / sqrt(n()),
    ci_nT_upper = mean(nT) + 1.96 * sd(nT) / sqrt(n()),
    
    # Mean and confidence intervals for "Number of occurrence"
    mean_occurrence = mean(`Number of occurrence`),
    ci_occurrence_lower = mean(`Number of occurrence`) - 1.96 * sd(`Number of occurrence`) / sqrt(n()),
    ci_occurrence_upper = mean(`Number of occurrence`) + 1.96 * sd(`Number of occurrence`) / sqrt(n()),
    
    # Mean and confidence intervals for Sobs
    mean_Sobs = mean(Sobs),
    ci_Sobs_lower = mean(Sobs) - 1.96 * sd(Sobs) / sqrt(n()),
    ci_Sobs_upper = mean(Sobs) + 1.96 * sd(Sobs) / sqrt(n()),
    
    # Mean and confidence intervals for Chat
    mean_Chat = mean(Chat),
    ci_Chat_lower = mean(Chat) - 1.96 * sd(Chat) / sqrt(n()),
    ci_Chat_upper = mean(Chat) + 1.96 * sd(Chat) / sqrt(n()),
    
    # Mean and confidence intervals for Good.u
    mean_Good_u = mean(Good.u),
    ci_Good_u_lower = mean(Good.u) - 1.96 * sd(Good.u) / sqrt(n()),
    ci_Good_u_upper = mean(Good.u) + 1.96 * sd(Good.u) / sqrt(n()),
    
    # Mean and confidence intervals for multition.rate
    mean_multition_rate = mean(multition.rate),
    ci_multition_rate_lower = mean(multition.rate) - 1.96 * sd(multition.rate) / sqrt(n()),
    ci_multition_rate_upper = mean(multition.rate) + 1.96 * sd(multition.rate) / sqrt(n()),
    
    # Mean and confidence intervals for Q1-Q10
    across(starts_with("Q"), list(
      mean = ~ mean(.),
      ci_lower = ~ mean(.) - 1.96 * sd(.) / sqrt(n()),
      ci_upper = ~ mean(.) + 1.96 * sd(.) / sqrt(n())
    ), .names = "{.col}_{.fn}")
  )