# HJA Dynamic Storage Analysis

Analysis of dynamic, mobile, and extended dynamic storage across watersheds at the H.J. Andrews Experimental Forest, Oregon.

## Project Overview

This repository contains R code for calculating storage metrics from streamflow, stream chemistry, and stream isotope data, and analyzing relationships between storage, catchment characteristics, and ecological responses (stream temperature, low-flow).

**Storage metrics calculated:**
- **Dynamic** — RBI (flashiness), recession curve slope, flow duration curve slope, storage-discharge
- **Mobile** — chemical hydrograph separation (EC-based), mean transit time, young water fraction, damping ratio
- **Extended dynamic** — water balance drawdown

**Sites:** 10 hydrometric watersheds (GSWS01–10, GSLOOK, GSWSMC) + 6 isotope-only sites


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

## Requirements

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "lubridate", "ggplot2",
  "zoo", "pracma", "colorspace", "scales", "patchwork",
  "GGally", "ggcorrplot", "vegan", "MASS", "car", "ggrepel"
))
```
---