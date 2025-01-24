# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 00_data_preparation.R
# Last updated: 2025-01-21
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; die.wen@ucl.ac.uk
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(palaeoverse)
library(dplyr)
library(dggridR)
source("./R/options.R")

# Data downloading from PBDB --------------------------------------------
if (params$download || !file.exists("./data/raw/pbdb_data.RDS")) {
  # Use for fresh downloads
  library(RCurl)
  library(httr)
  RCurl::curlSetOpt(3000)
  # Read data 
  url <- httr::modify_url(url = params$base_url, query = params$query)
  occdf <- RCurl::getURL(url = url, ssl.verifypeer = FALSE)
  occdf <- read.csv(textConnection(occdf))
  # Save raw data (do not modify to ensure complete reproducibility)
  saveRDS(occdf, "./data/raw/pbdb_data.RDS")
} else {
  # Read data
  occdf <- readRDS("./data/raw/pbdb_data.RDS")
}

# Set up time bins -----------------------------------------------------
bins <- time_bins(interval = "Phanerozoic",
                  rank = params$rank,
                  scale = params$GTS)
# Collapse Holocene equivalent bins
vec <- which(bins$interval_name == "Greenlandian")
bins$interval_name[vec] <- "Holocene"
bins$abbr[vec] <- "H"
# Update min_ma
bins$min_ma[vec] <- 0.0000
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Drop rows
bins <- bins[-which(bins$interval_name %in% c("Meghalayan", "Northgrippian")), ]
# Collapse Pleistocene equivalent bins
# Drop bins
pleis <- c("Late Pleistocene", "Chibanian", "Calabrian")
bins <- bins[-which(bins$interval_name %in% pleis), ]
# update Gelasian to be all of the Pleistocene
vec <- which(bins$interval_name == "Gelasian")
bins$interval_name[vec] <- "Pleistocene"
# Update min_ma
bins$min_ma[vec] <- bins[which(bins$interval_name == "Holocene"), "max_ma"]
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Update bin numbers
bins$bin <- 1:nrow(bins)
row.names(bins) <- 1:nrow(bins)
# Save time bins
saveRDS(object = bins, file = "./data/time_bins.RDS")

# Data cleaning and processing -----------------------------------------
# Clean up age assignments
# Remove any regular prefixes (Early, Middle, Late)
occdf$early_interval <- gsub(pattern = "Early ",
                             replacement = "",
                             x = occdf$early_interval)
occdf$late_interval <- gsub(pattern = "Early ",
                            replacement = "",
                            x = occdf$late_interval)
occdf$early_interval <- gsub(pattern = "Middle ",
                             replacement = "",
                             x = occdf$early_interval)
occdf$late_interval <- gsub(pattern = "Middle ",
                            replacement = "",
                            x = occdf$late_interval)
occdf$early_interval <- gsub(pattern = "Late ",
                             replacement = "",
                             x = occdf$early_interval)
occdf$late_interval <- gsub(pattern = "Late ",
                            replacement = "",
                            x = occdf$late_interval)

# Add required columns for look_up()
bins$early_stage <- bins$interval_name
bins$late_stage <- bins$interval_name
# Look up ages for intervals names using GTS2023/08
# Add ages
occdf <- look_up(occdf = occdf,
                 early_interval = "early_interval",
                 late_interval = "late_interval",
                 int_key = bins,
                 assign_with_GTS = FALSE)
# Which intervals could not be looked up?
vec_max <- which(is.na(occdf$interval_max_ma))
vec_min <- which(is.na(occdf$interval_min_ma))
# Use original input ages
occdf$interval_max_ma[vec_max] <- occdf$max_ma[vec_max]
occdf$interval_min_ma[vec_min] <- occdf$min_ma[vec_min]
occdf$early_stage[vec_max] <- occdf$early_interval[vec_max]
occdf$late_stage[vec_min] <- occdf$late_interval[vec_min]
# Replace max_ma and min_ma ages
occdf$max_ma <- occdf$interval_max_ma
occdf$min_ma <- occdf$interval_min_ma
# Calculate interval_mid_ma
occdf$interval_mid_ma <- (occdf$interval_max_ma + occdf$interval_min_ma) / 2
# Remove any collections with a large age range (> 50 Myr)
occdf <- occdf[-which(abs(occdf$max_ma - occdf$min_ma) > 50), ]

# Remove suffixes from genus names
occdf$genus <- sub(" .*", "", occdf$genus)

# Round off coordinates to stack collections 
occdf$lng <- round(occdf$lng, digits = params$n_decs)
occdf$lat <- round(occdf$lat, digits = params$n_decs)

# Temporal binning -----------------------------------------------------
# Use collections to speed up binning and palaeogeographic reconstruction
colldf <- unique(occdf[, c("collection_no", "lng", "lat", "max_ma", "min_ma")])
# Use the majority method
colldf <- bin_time(occdf = colldf, bins = bins, method = params$method)
# Remove data which do not hit the majority threshold (params$threshold)
colldf <- colldf[-which(colldf$overlap_percentage < params$threshold), ]

# Palaeorotate collections ---------------------------------------------
colldf <- palaeorotate(occdf = colldf,
                       lng = params$lng,
                       lat = params$lat,
                       age = params$age,
                       model = params$models,
                       method = "point",
                       uncertainty = FALSE,
                       round = NULL)
# Exclude collections which palaeocoordinates could not be estimated for
colldf <- subset(colldf, !is.na(colldf$p_lat))

# Spatial binning -------------------------------------------------------
# Construct a global grid with cells with approx params$spacing km
dggs <- dgconstruct(spacing = params$spacing)
# Get cells
colldf$cell <- dgGEO_to_SEQNUM(dggs = dggs, 
                         in_lon_deg = colldf[, params$p_lng], 
                         in_lat_deg = colldf[, params$p_lat])$seqnum
# Get coordinates from cells
xy <- dgcellstogrid(dggs = dggs, cells = colldf$cell, return_sf = FALSE)
# Rename columns
colnames(xy) <- c("cell_lng", "cell_lat", "cell_index")
# Merge dataframes
# Join datasets
m <- match(x = colldf$cell, table = xy$cell_index)
# Add data
colldf[, colnames(xy)] <- xy[m, colnames(xy)]
# Drop column
colldf <- colldf[, -ncol(colldf)]
# Join datasets --------------------------------------------------------
# Retain collections present in colldf
occdf <- occdf[which(occdf$collection_no %in% colldf$collection_no), ]
# Join datasets
m <- match(x = occdf$collection_no, table = colldf$collection_no)
# Add data
occdf[, colnames(colldf)] <- colldf[m, colnames(colldf)]
# Filter for unique occurrences from stacked collections
occdf <- distinct(occdf, lat, lng, family, genus, bin_assignment,
                  .keep_all = TRUE)
# Save processed data
saveRDS(object = occdf, file = "./data/processed/pbdb_data.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}
