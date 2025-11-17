### LDG_climate_state
Latitudinal Diversity Gradient (LDG) analysis workflow for deep-time biodiversity and climate-state reconstruction.
This repository provides a fully reproducible R pipeline for computing LDG curves, LDG slopes, and climate-state comparisons based on global fossil occurrence data.


### Main Functions
LDG_climate_state/
│
├── 000_Main.R                     # Master pipeline controller
│
├── 00_data_preparation.R          # Data cleaning, paleolat, time-bin assignment
├── 00_data_distribution_map.R     # Fossil distribution maps
│
├── 01_LDG_calculation.R           # Richness calculations (dggridR + iNEXT)
│
├── 02_LDG_slope.R                 # LDG slope estimation
├── 02_LDG_slope_fig.R             # Slope figure (main)
├── 02_LDG_slope_fig2.R            # Slope figure (variations)
├── 02_LDG_slope_fig3.R            # Slope figure (sensitivity)
├── 02_LDG_slope_sensitivity_test.R# Slope robustness tests
│
├── 03_LDG_compared_in_climate_state.R   # LDG × climate-state comparison
├── 03_LDG_sensitivity_test.R            # Additional sensitivity tests
│
├── 04_LDG_completeness_estimate.R       # Cell completeness analysis
├── 04_LDG_histogram.R                   # Richness distribution plots
│
├── calculate_Info.R               # Utility functions
├── calculate_LDG_slope.R          # Core slope-calculation function
├── check_hemisphere_good.R        # QC for hemisphere-level data sufficiency
│
├── options.R                      # Global settings and parameters
├── test_code.R                    # Testing & debugging
│
├── data/                          # PBDB, climate state, time bins, SC16
└── LICENSE


### Data
- `climate_states.csv`: climate state information from Judd et al., 2021.
- `time_bins.rds` : International stratigraphic chart color codes used for visualizations.
- `./data/raw/pbdb_data.RDS`: Raw data from PBDB used in the analysis.
- `./data/processed/data_clean.rds`: Cleaned and processed data ready for analysis.
- `./SC16`: paleogeographic surfaces from Scotese (2016).

### Additional Files
- `.gitignore`: Specifies files to be excluded from version control.
- `LICENSE`: License information for the project.

## Getting Started

### Prerequisites
- R (version 4.0 or above)
- R packages:
  `tidyverse`
  `data.table`
  `ggplot2`
  `dplyr`
  `iNEXT`
  `dggridR`
  `sf`
  `terra`
  `deeptime`

### Installation
Clone the repository to your local machine:
```bash
git clone <repository_url>
# 
# Usage
# 
# Simply execute the master script:
# source("000_Main.R")
# 
# This will automatically run all steps:
# Data cleaning & paleolatitude assignment
# Equal-area richness estimation
# LDG curve generation
# LDG slope analysis
# Climate-state comparison
# Completeness and sensitivity tests

### Pending items:

standardize the default iNEXT quorum (q)

confirm whether incidence or abundance should be used for PBDB occurrences

finalize the equal-area grid spacing (500 km vs 250 km etc sensitivity)

finalize bootstrap settings (nboot)