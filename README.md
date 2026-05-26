[![DOI](https://zenodo.org/badge/895013248.svg)](https://doi.org/10.5281/zenodo.17631059)


### LDG_climate_state
Latitudinal Diversity Gradient (LDG) analysis workflow for deep-time biodiversity and climate-state reconstruction.
This repository provides a fully reproducible R pipeline for computing LDG curves, LDG slopes, and climate-state comparisons based on global fossil occurrence data.


## Overview

The main goals of this project are to:

1. reconstruct sampling-standardised marine invertebrate richness across palaeolatitude;
2. estimate LDG slopes through geological time;
3. compare LDG strength among coldhouse, coolhouse, transitional, warmhouse, and hothouse climate states;
4. evaluate the robustness of LDG estimates to sampling, spatial coverage, grid size, and slope-calculation method;
5. test whether hemisphere-specific LDG differences are associated with sampling structure.

The current main slope-estimation approach is based on **per-cell balanced resampling OLS**, in which richness values from equal-area grid cells are resampled within palaeolatitudinal bands before fitting hemisphere-specific LDG slopes.


### Main Functions
LDG_climate_state/
- `000_Main_new.R`                     # Master pipeline controller
- `00_data_preparation.R`          # Data cleaning, paleolat, time-bin assignment
- `00_data_distribution_map.R`     # Fossil distribution maps
- `01_LDG_calculation.R`           # Richness calculations (dggridR + iNEXT)
- `02_LDG_slope_per_cell.R`                 # LDG slope estimation
- `02_LDG_slope_fig_per_cell.R`             # Slope figure (main)
- `02_LDG_slope_fig2_per_cell.R`            # Slope figure (variations)
- `02_LDG_slope_fig3_per_cell.R`            # Slope figure (sensitivity)
- `02b_LDG_slope_QC_sensitivity.R`    # Quality control (QC) robustness tests
- `02_03_LDG_slope_per_cell_allcells_climate_state.R`    # Slope calculation methods robustness tests
- `03_LDG_compared_in_climate_state_per_cell.R`   # LDG × climate-state comparison
- `04_LDG_completeness_estimate.R`       # Cell completeness analysis
- `04_LDG_histogram.R`                   # Richness distribution plots
- `05_high_latitude_coverage_summary.R`  # High latitude coverage test
  - `calculate_Info.R`               # Utility functions
  - `calculate_LDG_slope.R `         # Core slope-calculation function
  - `check_hemisphere_good.R`        # QC for hemisphere-level data sufficiency
- `options.R`                      # Global settings and parameters
- `test_code.R`                    # Testing & debugging
- `data/`                          # PBDB, climate state, time bins, SC16
- `LICENSE`


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
`source("000_Main.R")`
# 
This will automatically run all steps:
Data cleaning & paleolatitude assignment
Equal-area richness estimation
LDG curve generation
LDG slope analysis
Climate-state comparison
Completeness and sensitivity tests
```

### Pending items:

standardize the default iNEXT quorum (q)

confirm whether incidence or abundance should be used for PBDB occurrences

finalize the equal-area grid spacing (500 km vs 250 km etc sensitivity)

finalize bootstrap settings (nboot)