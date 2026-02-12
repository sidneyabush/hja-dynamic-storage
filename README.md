# HJA Dynamic Storage Analysis
Analysis of dynamic, mobile, and extended dynamic storage across watersheds at the H.J. Andrews Experimental Forest, Oregon.

## Project Overview

Quantifying storage metrics from streamflow, stream chemistry, and stream isotope data, and analyzing relationships between storage, watershed characteristics, and ecological responses (stream temperature, low-flow).

**Storage metrics calculated:**
- **Dynamic** — RBI (flashiness), recession curve slope, flow duration curve slope, storage-discharge
- **Mobile** — chemical hydrograph separation (EC-based), mean transit time, young water fraction, damping ratio
- **Extended dynamic** — water balance drawdown

## Requirements

```r
install.packages(c(
  "dplyr", "readr", "tidyr", "lubridate", "ggplot2",
  "zoo", "pracma", "colorspace", "scales", "patchwork",
  "GGally", "ggcorrplot", "vegan", "MASS", "car", "ggrepel"
))
```
---
