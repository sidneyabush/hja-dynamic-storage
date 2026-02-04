# HJA Dynamic Storage Analysis

Analysis of dynamic, mobile, and extended dynamic storage across watersheds at the H.J. Andrews Experimental Forest, Oregon.

## Project Overview

This repository contains R code for calculating storage metrics from streamflow, stream chemistry, and stream isotope data, and analyzing relationships between storage, catchment characteristics, and ecological responses (stream temperature, low-flow).

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

### Required Data Files

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

Raw data available from the [H.J. Andrews LTER Data Portal](https://andrewsforest.oregonstate.edu/data).

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

## Data Availability Summary

### Storage Metrics by Type

| Storage Type | Abbreviation | n | Sites | Years |
|--------------|--------------|---|-------|-------|
| **Dynamic** | RBI | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1950-2020 (varies) |
| | RCS | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1950-2020 (varies) |
| | FDC | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2020 |
| | SD | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2020 |
| **Mobile** | MTT | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS07, GSWS08, GSWSMC, MR | Site-level |
| | Fyw | 7 | GSWS01, GSWS02, GSWS07, GSWS08, GSWSMC, GSLOOK, MR | Site-level |
| | CHS | 8 | GSWS10, GSWS01, GSWS02, GSWS03, GSWS06*, GSWS07, GSWS08, GSWSMC | 2013-2019 (varies) |
| | DR | 15 | All hydrometric + MR, NC, LC, LO2, CC, LO1 | Site-level |
| **Extended Dynamic** | WB | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2019 |

*GSWS06 has limited CHS data (2017-2019); CHS excludes GSWS09 and GSLOOK (no EC data)

---

