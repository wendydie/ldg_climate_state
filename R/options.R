# options
# Required packages
required_packages <- c("httr", 
                       "readr", 
                       "stringdist", 
                       "divvy", 
                       "tictoc", 
                       "dplyr", 
                       "tidyr", 
                       "stringr", 
                       "iNEXT",
                       "ggplot2",
                       "pbapply",
                       "parallel")
invisible(lapply(required_packages, 
                 function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) 
                   install.packages(pkg)))

# Load libraries
invisible(lapply(required_packages, library, character.only = TRUE))

data_clean <- readRDS("data/data_clean.rds")

# Source functions
source('./R/functions/buffer_subsampling.R')
source("./R/functions/calculate_Info.R")

# 0.1 data downloading from PBDB
taxa <- c("Bivalvia", 
          "Brachiopoda", 
          "Gastropoda", 
          "Trilobita",
          "Bryozoa", 
          "Echinoidea")
max_ma <- 500
min_ma <- 0
base_url <- "https://paleobiodb.org/data1.2/occs/list.csv"
query <- list(
  datainfo = "",
  rowcount = "",
  base_name = paste(taxa, collapse = ","),
  taxon_reso = "genus",
  max_ma = max_ma,
  min_ma = min_ma,
  pgm = "gplates,scotese,seton",
  show = "full"
)

# 0.2 data cleaning
# Filter out missing values
selected_columns <- c("occurrence_no", 
                      "reference_no",
                      "collection_no", 
                      "genus", 
                      "lat", 
                      "lng",
                      "max_ma",
                      "min_ma", 
                      "early_interval", 
                      "late_interval",
                      "paleomodel2",
                      "paleolng2", 
                      "paleolat2")

# 1.1 LDG calculation
# Parameters
xy = c("paleolat2", "paleolng2") # Columns representing geographic coordinates (latitude and longitude)
nSite = 1                        # Number of sites to sample within each buffer
r = 200                          # Buffer radius in kilometers
crs = 'epsg:4326'                # Coordinate reference system, default is WGS84
q = c(0)                         # Diversity order(s) to calculate; 0 for species richness
datatype = "incidence_freq"      # Data type; 'incidence_freq' for presence-absence data aggregated into frequencies
base = "coverage"                # Standardization base; 'coverage' ensures comparisons at the same coverage level
level = 0.7                      # Desired sample coverage level (e.g., 0.7 = 70% of species expected to be detected)
nboot = 50                        # Number of bootstrap replicates for confidence intervals; 0 means no bootstrapping
stage_name = NULL 
stage_mid = NULL
stbin = NULL
source("./R/01_LDG_calculation.R")
