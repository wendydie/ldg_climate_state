# LDG_climate_state
# LDG Climate State Analysis

This repository contains scripts and data for analyzing latitudinal diversity gradients (LDG) in relation to climate states, with a focus on Paleozoic periods. The project employs a combination of subsampling techniques and comparative analysis to explore the relationships between climate, biodiversity, and geographical distribution.

## Repository Structure

### Functions
- `buffer_subsampling.R`: Implements a buffer subsampling method for data analysis.
- `calculate_Info.R`: Contains functions to calculate relevant metrics and information for LDG analysis.
- `00_data_preparation.R`: Prepares raw data for analysis, including cleaning and formatting.
- `01_LDG_calculation.R`: Calculates latitudinal diversity gradients using various methods.
- `02_LDG_slope.R`: Analyzes the slopes of LDG curves to assess trends.
- `03_LDG_compared_in_climate_state.R`: Compares LDG results across different climate states.
- `test_code.R`: A script for testing and debugging functions.
- `options.R`: Configuration and global options for the project.

### Data
- `CGMW_ICS_colour_codes.xlsx`: International stratigraphic chart color codes used for visualizations.
- `data.rds`: Raw data used in the analysis.
- `data_clean.rds`: Cleaned and processed data ready for analysis.
- `stages.csv`: Information about geological stages and time periods.

### Additional Files
- `.gitignore`: Specifies files to be excluded from version control.
- `LICENSE`: License information for the project.

## Getting Started

### Prerequisites
- R (version 4.0 or above)
- R packages:
  - `tidyverse`
  - `ggplot2`
  - `dplyr`
  - `readr`
  - `data.table`

### Installation
Clone the repository to your local machine:
```bash
git clone <repository_url>


## Current State and Pending Issues

### 1. Paleolatitude Reconstruction Annotation
In the script `00_data_preparation.R`, the section labeled **# 0.3 To-Do: Paleolatitude Reconstruction Annotation** is currently undefined. The parameters and methods for this part of the analysis need to be determined.

### 2. Richness Results Calculation
The following code snippet in `00_data_preparation.R` calculates richness results using a set of parameters, but some aspects require further discussion:

```R
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