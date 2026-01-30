# HJA Dynamic Storage Analysis

Analysis of dynamic and mobile water storage across watersheds at the H.J. Andrews Experimental Forest, Oregon.

## Project Overview

This repository contains R code for calculating storage metrics from hydrometric, chemistry, and isotope data, and analyzing relationships between storage, catchment characteristics, and ecological responses (stream temperature, low-flow).

**Timeline:** Water Years 1997-2020 (Oct 1, 1996 - Sep 30, 2020)

### Storage Metrics

| Type | Metric | Description | Data Source |
|------|--------|-------------|-------------|
| **Dynamic** | RBI | Richards-Baker Flashiness Index | Daily discharge |
| **Dynamic** | RCS | Recession Curve Slope | Daily discharge |
| **Dynamic** | FDC | Flow Duration Curve Slope | Daily discharge |
| **Dynamic** | SD | Storage-Discharge (S_annual_mm) | Daily P, Q, ET |
| **Mobile** | CHS | Chemical Hydrograph Separation (mean_bf) | Specific conductance |
| **Mobile** | MTT | Mean Transit Time | Isotopes |
| **Mobile** | Fyw | Young Water Fraction | Isotopes |
| **Mobile** | DR | Isotopic Damping Ratio | Isotopes |
| **Extended Dynamic** | WB | Water Balance Drawdown (DS_sum) | Daily P, Q, ET |

### Study Sites

- **Hydrometric sites (10):** GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC
- **Isotope-only sites (6):** MR (McRae), NC (Nostoc), LC (Longer), LO1, LO2, CC (Cold Creek)

---

## Repository Structure

```
hja-dynamic-storage/
├── config.R                        # Central configuration (paths, sites, constants)
├── data/                           # Data folder (files not tracked due to size)
│   └── README.md                   # Data file requirements
├── 00_Data_Preprocessing/          # Meteorological data preparation
│   └── Create_Master_Hydrometric_Dataset.R
├── 01_Dynamic_Storage/             # Dynamic storage metrics
│   ├── 01_RBI_Recession.R          # RBI and recession curve slope
│   └── 02_Storage_Discharge_FDC.R  # Storage-discharge and FDC
├── 02_Mobile_Storage/              # Mobile storage metrics
│   ├── 01_Load_Isotope_Metrics.R   # MTT, Fyw, DR from isotope data
│   └── 02_Chemical_Hydrograph_Separation.R  # CHS from EC data
├── 03_Extended_Dynamic/            # Extended dynamic storage
│   └── 01_Water_Balance_Drawdown.R # Water balance drawdown
├── 04_Response_Variables/          # Ecological response variables
│   └── 01_Stream_Temperature_LowFlow.R
├── 05_Aggregate_Metrics/           # Combine all metrics
│   ├── 00_Data_Availability_Summary.R
│   └── 01_Aggregate_All_Metrics.R
├── 06_Statistical_Analyses/        # Statistical analyses
│   ├── 01_Correlations.R
│   ├── 02_PCA.R
│   ├── 03_ANOVA_Tukey.R
│   ├── 04_Catchment_Storage_MLR.R
│   └── 06_Storage_Predicts_Thermal_LowFlow.R
└── 07_Plots/                       # Publication figures
    ├── 00_Run_All_Plots.R          # Master script to generate all figures
    ├── Hydrometric_Plots.R         # Time series and summary plots
    ├── Recession_Curves.R          # Recession curve log-log plots
    └── Publication_Figures.R       # Final publication-ready figures
```

---

## Getting Started

### 1. Configure Paths

Edit `config.R` to set your data paths:

```r
# Option 1: Local data in repo data/ folder
USE_LOCAL_DATA <- TRUE

# Option 2: Box cloud storage (update path for your system)
USE_LOCAL_DATA <- FALSE
BOX_BASE_DIR <- "/path/to/your/Box/05_Storage_Manuscript"
```

### 2. Required Data Files

See `data/README.md` for complete list. Key files:

| Folder | File | Description |
|--------|------|-------------|
| Q/ | `HF00402_v14.csv` | Daily mean discharge for all watersheds |
| Q/ | `drainage_area.csv` | Drainage areas for each site |
| DynamicStorage/ | `daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv` | Daily P, Q, ET |
| DynamicStorage/ | `Catchment_Charc.csv` | Catchment characteristics |
| EC/ | `CF01201_v3.txt` | Continuous specific conductance |
| Isotopes/ | `MTT_FYW.csv` | Mean transit times and young water fractions |
| Isotopes/ | `DampingRatios_2025-07-07.csv` | Isotopic damping ratios |
| Stream_T/ | `HT002*.csv` | Stream temperature files |

Data available from the [H.J. Andrews LTER Data Portal](https://andrewsforest.oregonstate.edu/data).

### 3. Run Analysis Workflow

Execute scripts in order:

```r
# 1. Preprocess meteorological data (if needed)
source("00_Data_Preprocessing/Create_Master_Hydrometric_Dataset.R")

# 2. Calculate storage metrics
source("01_Dynamic_Storage/01_RBI_Recession.R")
source("01_Dynamic_Storage/02_Storage_Discharge_FDC.R")
source("02_Mobile_Storage/01_Load_Isotope_Metrics.R")
source("02_Mobile_Storage/02_Chemical_Hydrograph_Separation.R")
source("03_Extended_Dynamic/01_Water_Balance_Drawdown.R")

# 3. Calculate response variables
source("04_Response_Variables/01_Stream_Temperature_LowFlow.R")

# 4. Aggregate all metrics
source("05_Aggregate_Metrics/01_Aggregate_All_Metrics.R")

# 5. Run statistical analyses
source("06_Statistical_Analyses/01_Correlations.R")
source("06_Statistical_Analyses/02_PCA.R")
source("06_Statistical_Analyses/03_ANOVA_Tukey.R")
source("06_Statistical_Analyses/04_Catchment_Storage_MLR.R")
source("06_Statistical_Analyses/06_Storage_Predicts_Thermal_LowFlow.R")

# 6. Generate publication figures
source("07_Plots/00_Run_All_Plots.R")
```

---

## Publication Figures

The `07_Plots/` folder generates all figures for the manuscript:

| Script | Figures Generated |
|--------|-------------------|
| `Hydrometric_Plots.R` | Time series by site, summary plots with Tukey letters, grid of all metrics |
| `Recession_Curves.R` | Log-log recession curves (-dQ/dt vs Q), recession slope comparison |
| `Publication_Figures.R` | Correlation matrix, PCA biplot, storage-response relationships |

Run all plots with:
```r
source("07_Plots/00_Run_All_Plots.R")
```

---

## Required R Packages

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "lubridate", "ggplot2",
  "zoo", "pracma", "colorspace", "scales", "patchwork",
  "GGally", "ggcorrplot", "vegan", "MASS", "car", "ggrepel"
))
```

---

## Key Outputs

| File | Description |
|------|-------------|
| `HJA_StorageMetrics_Annual_All.csv` | Annual storage metrics by site |
| `HJA_Ave_StorageMetrics_CatCharacter.csv` | Site-averaged metrics with catchment characteristics |
| `Sample_Size_by_Metric.csv` | Data availability summary |
| `RBI_RecessionCurve_Annual.csv` | Annual RBI and recession slopes |
| `StorageDischarge_FDC_Annual.csv` | Storage-discharge and FDC metrics |
| `DS_drawdown_annual.csv` | Annual dynamic storage drawdown |
| `Annual_GW_Prop.csv` | Annual mean baseflow proportion |
| `stream_thermal_lowflow_metrics_annual.csv` | Thermal and low-flow response variables |

---

## Data Availability Summary

### Storage Metrics by Type

| Storage Type | Abbreviation | Sites |
|--------------|--------------|-------|
| **Dynamic** | RBI | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC (10) |
| | RCS | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC (10) |
| | FDC | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC (10) |
| | SD | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC (10) |
| **Mobile** | MTT | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS07, GSWS08, GSWSMC, MR (10) |
| | Fyw | GSWS01, GSWS02, GSWS07, GSWS08, GSWSMC, GSLOOK, MR (7) |
| | CHS | GSWS10, GSWS01, GSWS02, GSWS03, GSWS06*, GSWS07, GSWS08, GSWSMC (8) |
| | DR | 15 sites (all hydrometric + MR, NC, LC, LO2, CC, LO1) |
| **Extended Dynamic** | WB | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC (10) |

*GSWS06 has limited CHS data (2017-2019); CHS excludes GSWS09 and GSLOOK (no EC data)

### Mean Transit Time (MTT) and Young Water Fraction (Fyw)

| Site | MTT1 (yr) 2001-2003* | MTT2 (yr) 2015-2018** | Fyw 2015-2018** |
|--------|----------------------|------------------------|-----------------|
| GSWS09 | 0.8 (0.18) | - | - |
| GSWS10 | 1.2 (0.29) | - | - |
| GSWS01 | - | 1.99 (0.44) – 4.14 (1.1) | 0.14 (0.07) – 0.32 (0.21) |
| GSLOOK | 2.0 (1.00) | 1.11 (0.15) – 5.01 (0.1) | 0.08 (0.03) – 0.54 (0.15) |
| GSWS02 | 2.2 (0.56) | 1.3 (0.11) – 4.41 (1.79) | 0.04 (0.01) – 0.39 (0.2) |
| GSWS03 | 1.3 (0.31) | - | - |
| GSWS06 | - | - | - |
| GSWS07 | - | 0.95 (0.08) – 5.14 (0.07) | 0.04 (0) – 0.91 (0.2) |
| GSWS08 | 3.3 (1.28) | 0.95 (0.04) – 3.28 (1.37) | 0.05 (0.01) – 0.61 (0.16) |
| GSWSMC | 2.0 (0.49) | 0.89 (0.06) – 5.55 (0.67) | 0.06 (0.02) – 0.47 (0.13) |
| MR | - | 1.45 (0.32) – 5.03 (1.14) | 0.05 (0) – 0.2 (0.02) |

*McGuire, 2005; **Segura, 2021. Values show mean (SE); ranges indicate different estimation methods.

### Temporal Coverage

| Metric | Sites | Years | Notes |
|--------|-------|-------|-------|
| RBI, RCS | 10 | 1950-2020 (varies) | Full hydrometric record |
| FDC, SD | 10 | 1998-2020 | Water balance period |
| CHS | 8 | 2013-2019 (varies) | EC data availability |
| WB | 10 | 1998-2019 | Water balance period |
| MTT | 11 | Site-level | Single value per site |
| Fyw | 8 | Site-level | Not all sites have data |
| DR | 15 | Site-level | Includes isotope-only sites |

---

