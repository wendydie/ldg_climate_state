# function 1. calculating Good'U index
incfreq <- function(dat){
  freq <- dat %>%
    count(genus, name = "frequency") %>%
    arrange(desc(frequency))
  nT <- nrow(dat)
  y <- c(nT, freq$frequency)
  y <- list(as.numeric(y))
  return(y)
}


Fun.ince <- function(x) {
  nT <- x[1]  # The total sample size, which is the first element in the vector x
  x <- x[-1]  # The remaining elements in x are the counts of species or occurrences
  U <- sum(x)  # Total number of occurrences (sum of all remaining elements)
  # Calculate Qk: a vector where Qk[k] represents the number of species occurring exactly k times
  Qk <- sapply(1:10, function(k) sum(x == k))
  Q1 <- Qk[1]  # Number of species that occur exactly once
  Q2 <- Qk[2]  # Number of species that occur exactly twice
  Sobs <- sum(x > 0)  # Number of observed species (species with at least one occurrence)
  # Estimate the number of undetected species (Q0.hat)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2, 
                   (nT - 1)/nT * Q1^2/2/Q2)
  # Calculate A, which is used to adjust the estimate of sample coverage
  A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1), 1)
  # Calculate the estimated sample coverage (Chat), rounded to four decimal places
  Chat <- round(1 - Q1/U * A, 4)
  Good.u <- 1 - Q1/U
  multition.rate <- (Sobs - Q1) / Sobs # cited from John Alroy on 21 March 2017
  # Return a vector of results including nT, U, Sobs, Chat, and Qk (Q1 to Q10)
  out <- c(nT, U, Sobs, Chat, Good.u, multition.rate, Qk)
}


# function 3. Common Function to count references, collections and occurrences
count_occ_col_refs <- function(dat) {
  # Ensure required columns exist
  required_columns <- c("occurrence_no", "collection_no", "reference_no")
  if (!all(required_columns %in% colnames(dat))) {
    stop("Missing required columns.")
  }
  # Calculate unique counts for each column
  sapply(dat[required_columns], function(col) length(unique(col)))
}


#function 4.
compute_richness_summary <- function(dat, xy, iter = 50, nSite = 10, r = 1500, 
                                     crs = 'epsg:4326', q = c(0), 
                                     datatype = "incidence_freq", 
                                     base = "coverage", level = 0.7, nboot = 0,
                                     stage_name = NULL, 
                                     stage_mid = NULL, stbin = NULL) {
  # Step 1: Check for stage information in dat if not provided
  if (is.null(stage_name) && "stage_name" %in% colnames(dat)) {
    stage_name <- unique(dat$stage)
    if (length(stage_name) > 1) {
      stop("Multiple stage names found in dat. Please provide a specific stage_name.")
    }
  }
  
  if (is.null(stage_mid) && "stage_mid" %in% colnames(dat)) {
    stage_mid <- unique(dat$stage_mid)
    if (length(stage_mid) > 1) {
      stop("Multiple stage times found in dat. Please provide a specific stage_mid.")
    }
  }
  
  if (is.null(stbin) && "stbin" %in% colnames(dat)) {
    stbin <- unique(dat$stbin)
    if (length(stbin) > 1) {
      stop("Multiple stbin values found in dat. Please provide a specific stbin.")
    }
  }
  
  # Step 2: Generate incidence frequencies using buffers
  incfreq <- buffers(dat, xy, iter = iter, nSite = nSite, r = r, crs = crs, 
                     output = 'incidence_freq')
  
  # Step 3: Compute richness metrics for each buffer using iNEXT::estimateD
  rich_list <- lapply(incfreq, function(buf) {
    iNEXT::estimateD(x = buf, q = q, datatype = datatype, 
                     base = base, level = level, nboot = nboot)
  })
  
  # Combine richness results across all iterations
  combined_rich <- do.call(rbind, lapply(seq_along(rich_list), function(i) {
    rich_list[[i]] %>% mutate(iter = i) # Add iteration ID
  }))
  
  # Add stage information to combined_rich
  combined_rich <- combined_rich %>%
    mutate(stage_name = stage_name, stage_mid = stage_mid, stbin = stbin)
  
  # Step 4: Apply Fun.ince to all buffers
  fun_results <- lapply(incfreq, function(iter) {
    lapply(names(iter), function(assem) {
      res <- Fun.ince(iter[[assem]])
      c(assem = assem, res) # Add Assemblage name
    })
  })
  
  # Combine results across all iterations into a single data frame
  rich_df <- do.call(rbind, lapply(seq_along(fun_results), function(i) {
    iter_res <- do.call(rbind, fun_results[[i]])
    cbind(iter = i, iter_res) # Add iteration ID
  }))
  
  # Convert to data frame and ensure numeric columns
  rich_df <- as.data.frame(rich_df)
  colnames(rich_df) <- c("Iteration", "assem", "nT", "n_occ",
                                 "Sobs", "Chat", "Good.u", "m.rate",
                                 paste0("Q", 1:10))
  numeric_cols <- c("nT", "n_occ", "Sobs", "Chat", "Good.u", "m.rate", paste0("Q", 1:10))
  rich_df[numeric_cols] <- lapply(rich_df[numeric_cols], as.numeric)
  
  # Add stage information to rich_df
  rich_df <- rich_df %>%
    mutate(stage_name = stage_name, stage_mid = stage_mid, stbin = stbin)
  
  # Step 5: Summarize results by Assemblage across iterations
  rich_info_summary <- rich_df %>%
    group_by(assem) %>%
    summarise(
      mean_nT = mean(nT, na.rm = TRUE),
      ci_nT_lower = mean(nT, na.rm = TRUE) - 1.96 * sd(nT, na.rm = TRUE) / sqrt(sum(!is.na(nT))),
      ci_nT_upper = mean(nT, na.rm = TRUE) + 1.96 * sd(nT, na.rm = TRUE) / sqrt(sum(!is.na(nT))),
      
      mean_n_occ = mean(n_occ, na.rm = TRUE),
      ci_n_occ_lower = mean(n_occ, na.rm = TRUE) - 1.96 * sd(n_occ, na.rm = TRUE) / sqrt(sum(!is.na(n_occ))),
      ci_n_occ_upper = mean(n_occ, na.rm = TRUE) + 1.96 * sd(n_occ, na.rm = TRUE) / sqrt(sum(!is.na(n_occ))),
      
      mean_Sobs = mean(Sobs, na.rm = TRUE),
      ci_Sobs_lower = mean(Sobs, na.rm = TRUE) - 1.96 * sd(Sobs, na.rm = TRUE) / sqrt(sum(!is.na(Sobs))),
      ci_Sobs_upper = mean(Sobs, na.rm = TRUE) + 1.96 * sd(Sobs, na.rm = TRUE) / sqrt(sum(!is.na(Sobs))),
      
      mean_Chat = mean(Chat, na.rm = TRUE),
      ci_Chat_lower = mean(Chat, na.rm = TRUE) - 1.96 * sd(Chat, na.rm = TRUE) / sqrt(sum(!is.na(Chat))),
      ci_Chat_upper = mean(Chat, na.rm = TRUE) + 1.96 * sd(Chat, na.rm = TRUE) / sqrt(sum(!is.na(Chat))),
      
      mean_Good_u = mean(Good.u, na.rm = TRUE),
      ci_Good_u_lower = mean(Good.u, na.rm = TRUE) - 1.96 * sd(Good.u, na.rm = TRUE) / sqrt(sum(!is.na(Good.u))),
      ci_Good_u_upper = mean(Good.u, na.rm = TRUE) + 1.96 * sd(Good.u, na.rm = TRUE) / sqrt(sum(!is.na(Good.u))),
      
      mean_m_rate = mean(m.rate, na.rm = TRUE),
      ci_m_rate_lower = mean(m.rate, na.rm = TRUE) - 1.96 * sd(m.rate, na.rm = TRUE) / sqrt(sum(!is.na(m.rate))),
      ci_m_rate_upper = mean(m.rate, na.rm = TRUE) + 1.96 * sd(m.rate, na.rm = TRUE) / sqrt(sum(!is.na(m.rate))),
      
      across(starts_with("Q"), list(
        mean = ~ mean(., na.rm = TRUE),
        ci_lower = ~ mean(., na.rm = TRUE) - 1.96 * sd(., na.rm = TRUE) / sqrt(sum(!is.na(.))),
        ci_upper = ~ mean(., na.rm = TRUE) + 1.96 * sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))
      ), .names = "{.col}_{.fn}")
    ) %>%
    mutate(stage_name = stage_name, stage_mid = stage_mid, stbin = stbin)
  
  # Return combined richness results and the summary
  list(combined_rich = combined_rich, rich_info_summary = rich_info_summary)
}