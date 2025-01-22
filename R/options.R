# Header ----------------------------------------------------------------
# Project: LDG_climate_state
# File name: options.R
# Last updated: 2025-01-21
# Author: Lewis A. Jones; Die (Wendy) Wen
# Email: lewis.jones@ucl.ac.uk; die.wen@ucl.ac.uk
# Repository: https://github.com/wendydie/LDG_climate_state
# -----------------------------------------------------------------------
# Parameters and options for analyses
params <- list(
  # Download updated dataset
  download = FALSE,
  # PBDB base url
  base_url = "https://paleobiodb.org/data1.2/occs/list.csv",
  # PBDB Download options
  query = list(
    base_name = paste(c("Bivalvia", "Brachiopoda", "Gastropoda", 
                        "Trilobita", "Bryozoa", "Echinoidea"), 
                      collapse = ","),
    taxon_reso = "genus",
    ident = "latest",
    taxon_status = "valid",
    idqual = "genus_certain",
    pres = "regular",
    interval = "Fortunian,Holocene",
    envtype = "marine",
    show = "genus,pres,strat,coll,coords,loc,class"),
  # The geological rank for conducting analyses
  rank = "stage",
  # The Geological Time Scale to be used
  GTS = "international ages",
  # Naming convention for temporal bin midpoint
  age = "bin_midpoint",
  # How should occurrences be temporally binned?
  method = "majority",
  # Threshold for majority binning rule
  threshold = 50,
  # Naming convention for longitude
  lng = "lng",
  # Naming convention for latitude
  lat = "lat",
  # Naming convention for palaeolongitude
  p_lng = "p_lng",
  # Naming convention for palaeolatitude
  p_lat = "p_lat",
  # Number of lat/lng decimal places to define stacked collections
  n_decs = 2,
  # Which Global Plate Models should be used?
  models = c("PALEOMAP"),
  # Columns representing geographic coordinates (latitude and longitude)
  xy = c("paleolat2", "paleolng2"),
  # Number of sites to sample within each buffer
  nSite = 1, 
  # Buffer radius in kilometers
  r = 200, 
  # Coordinate reference system, default is WGS84
  crs = 'epsg:4326', 
  # Diversity order(s) to calculate; 0 for species richness
  q = c(0),
  # Data type; 'incidence_freq' for presence-absence data aggregated into frequencies
  datatype = "incidence_freq", 
  # Standardization base; 'coverage' ensures comparisons at the same coverage level
  base = "coverage", 
  # Desired sample coverage level (e.g., 0.7 = 70% of species expected to be detected)
  level = 0.7, 
  # Number of bootstrap replicates for confidence intervals; 0 means no bootstrapping
  nboot = 50, 
  stage_name = NULL, 
  stage_mid = NULL,
  stbin = NULL,
  notify = TRUE
)


# # Required packages
# required_packages <- c("httr", 
#                        "readr", 
#                        "stringdist", 
#                        "divvy", 
#                        "tictoc", 
#                        "dplyr", 
#                        "tidyr", 
#                        "stringr", 
#                        "iNEXT",
#                        "ggplot2",
#                        "pbapply",
#                        "parallel")
# invisible(lapply(required_packages, 
#                  function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) 
#                    install.packages(pkg)))
# 
# # Load libraries
# invisible(lapply(required_packages, library, character.only = TRUE))
# 
# data_clean <- readRDS("data/data_clean.rds")
# 
# # Source functions
# source('./R/functions/buffer_subsampling.R')
# source("./R/functions/calculate_Info.R")
# 
# 
# # 0.2 data cleaning
# # Filter out missing values
# selected_columns <- c("occurrence_no", 
#                       "reference_no",
#                       "collection_no", 
#                       "genus", 
#                       "lat", 
#                       "lng",
#                       "max_ma",
#                       "min_ma", 
#                       "early_interval", 
#                       "late_interval",
#                       "paleomodel2",
#                       "paleolng2", 
#                       "paleolat2")
# 
# source("./R/01_LDG_calculation.R")
