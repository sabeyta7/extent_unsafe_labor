# Geographic Variation in Labor vs. Wage/Hour Violations - Hurdle Model

Analysis code for "Prevalence and Extent of Unsafe Labor in the United States"

## Data Access

The analysis dataset is available on Zenodo: **[10.5281/zenodo.17335395]**

**To run this analysis:**

1. Download the dataset from Zenodo (link above)
2. Extract all shapefile components (.shp, .shx, .dbf, .prj, .cpg)
3. Place all files in the `data/` folder of this repository
4. Run the main analysis script

## System Requirements

This analysis involves computationally intensive spatial models using INLA. 

**Successfully tested on:**
- MacBook Pro with Apple M2 Pro chip (10-core)
- 32 GB RAM
- R version 4.4.1
- macOS 15.6.1
- Runtime: Approximately 2-4 hours for full analysis

**Minimum recommended:**
- 16 GB RAM
- Multi-core processor
- R 4.0+

**Note:** The continuous interaction models (Section 7) may fail on university-managed laptops with artificial processing constraints.

## Required R Packages

dplyr, ggplot2, sf, tidyr, INLA, spdep, spatstat, viridis, RColorBrewer, corrplot, cowplot, patchwork, gridExtra, knitr, kableExtra, pROC, DCluster, survey

## Running the Analysis

```r
source("code/extent_unsafe_labor.R")
```

## Project Structure

```
├── code/
│   ├── extent_unsafe_labor.R  # Main analysis script
│   └── functions.R            # Functions for main analyses
├── data/                      # Place downloaded shapefiles here
│   ├── labor_data.shp
│   ├── labor_data.shx
│   ├── labor_data.dbf
│   ├── labor_data.prj
│   └── labor_data.cpg
└── README.md
```

## Analysis Overview

1. **Data Preparation** - Creates spatial neighborhoods and subsets
2. **Descriptives** - Summary statistics by industry and region
3. **Regional Profiles** - Industry-specific rates by region
4. **Relative Risk Analysis** - Industry-outcome correlations
5. **Bivariate Risk Maps** - Spatial visualization
6. **Prevalence Estimation** - SMRs with Empirical Bayes smoothing
7. **Regional Variation** - Industry-region interaction models (1-2 hours runtime)