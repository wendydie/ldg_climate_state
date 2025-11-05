# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: 00_data_preparation.R
# Last updated: 2025-10-15
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; geowendywen@outlook.com
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Load libraries and options --------------------------------------------
library(palaeoverse)
library(dplyr)
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
  # Read data 576848
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
bins$abbr[vec] <- "Ple"
# Update min_ma
bins$min_ma[vec] <- bins[which(bins$interval_name == "Holocene"), "max_ma"]
# Update mid_ma
bins$mid_ma[vec] <- (bins$min_ma[vec] + bins$max_ma[vec]) / 2
# Update duration
bins$duration_myr[vec] <- (bins$max_ma[vec] - bins$min_ma[vec])
# Update bin numbers
bins$bin <- 1:nrow(bins)
row.names(bins) <- 1:nrow(bins)

## GTS 2023
GTS_2023 <- read.csv('./data/GTS_2023.csv')
bins <- bins %>%
  left_join(GTS_2023, by="bin") %>%
  mutate(
    max_ma = bottom,
    mid_ma = mid,
    min_ma = top,
    duration_myr = dur
  ) %>%
  select(bin, interval_name, rank, max_ma, mid_ma, min_ma,
         duration_myr, short, font, sys, system, series,
         systemCol, seriesCol, stageCol, stageRGB)
# bins <- bins %>% filter(mid_ma < 485.4)
# Save time bins
saveRDS(object = bins, file = "./data/time_bins.RDS")

# # Data cleaning and processing -----------------------------------------

# Remove suffixes from genus names
occdf$genus <- sub(" .*", "", occdf$genus)

## Remove rows where the 'genus' column contains the value "NO_GENUS_SPECIFIED" 575880
occdf <- occdf %>%
filter(genus != "NO_GENUS_SPECIFIED")
## Remove rows where any of the columns 'genus', 'lat', or 'lng' have NA values
occdf <- occdf %>%
  filter(!is.na(genus) & !is.na(lng) & !is.na(lat))

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

# Join datasets --------------------------------------------------------
# Retain collections present in colldf 529008
occdf <- occdf[which(occdf$collection_no %in% colldf$collection_no), ]
# Join datasets
m <- match(x = occdf$collection_no, table = colldf$collection_no)
# Add data
occdf[, colnames(colldf)] <- colldf[m, colnames(colldf)]
# Filter for unique occurrences from stacked collections 461054 
occdf <- distinct(occdf, lat, lng, family, genus, bin_assignment, collection_no,
                  .keep_all = TRUE)

# Save processed data
saveRDS(object = occdf, file = "./data/processed/pbdb_data.RDS")
# Notify
if (params$notify) {
  beepr::beep(4)
}
