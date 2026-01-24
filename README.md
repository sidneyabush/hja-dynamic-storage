# HJA Storage Manuscript - Code Repository

Code and workflow for the H.J. Andrews Experimental Forest hydrometric storage analysis manuscript.

**Timeline:** Water Years 1997-2020
**Sites:** 10 HJA watersheds (GSLOOK, GSWS01-10, GSWSMC)

---

## Repository Structure

### [`Create_Master_Hydrometric_Dataset/`](Create_Master_Hydrometric_Dataset/)
Comprehensive workflow harmonizing meteorology and streamflow for ET calculations and downstream analyses.

### [`ET_calculations/`](ET_calculations/)
Evapotranspiration estimation using Hamon method with Zhang-style calibrated coefficient.

### [`Discharge_metrics/`](Discharge_metrics/)
Discharge-based storage indicators and ecological response variables:
- **RBI (flashiness)** - Richards-Baker Flashiness Index
- **Recession curve slope** - log-log slope of -dQ/dt vs Q
- **Dynamic storage** - Kirchner-Staudinger storage-discharge approach
- **Flow Duration Curves (FDC)** - Q99, Q50, Q01 percentiles
- **Chemical Hydrograph Separation** - Baseflow from specific conductance
- **Stream Temperature & Low-Flow Metrics** - 7-day moving averages for ecological stress indicators

### [`Statistical_analyses/`](Statistical_analyses/)
Statistical models and multivariate analyses:
- **Correlations** - Correlation matrices of storage metrics and catchment attributes
- **PCA** - Principal component analysis of storage metrics
- **MLR** - Multiple linear regression (catchment → storage)
- **RDA** - Redundancy analysis and variance partitioning
- **Storage → Thermal/Low-Flow** - Prediction models testing storage effects on ecological responses

### [`Plots/`](Plots/)
Figures for exploration and publication.

### [`Prior_analyses/`](Prior_analyses/)
Archived scripts from early project stages.

### [`deprecated/`](deprecated/)
Old or replaced scripts retained for reference.

---

## Workflow Overview

### Phase 1: Data Harmonization & ET
**Scripts:** `Create_Master_Hydrometric_Dataset/`, `ET_calculations/`

1. **Meteorology Harmonization**
   - Combine multi-station data using OLS (pairs) or multiple regression (triplets)
   - Minimum overlap: pairs ≥5 days, triplets ≥10 days
   - Output: `watersheds_met_data_q.csv`

2. **ET Calculation**
   - Hamon PET with Zhang calibration coefficient
   - Join P, Q, ET into daily water balance
   - Output: `daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv`

### Phase 2: Storage Metrics Calculation
**Scripts:** `Discharge_metrics/`

Run in order:
1. `01_RBI_Recession.R` - Calculate RBI and recession slope
2. `02_Storage_Discharge_FDC.R` - Kirchner-Staudinger storage-discharge, FDC
3. `03_Dynamic_Storage_Drawdown.R` - Dynamic storage drawdown from water balance
4. `04_Chemical_Hydrograph_Separation.R` - Baseflow from specific conductance
5. `05_Stream_Temperature_LowFlow.R` - Thermal and low-flow ecological indicators
6. `06_Aggregate_All_Metrics.R` - Combine all metrics into master tables

**Outputs:**
- `RBI_RecessionCurve_Annual.csv` - Annual RBI and recession slopes
- `StorageDischarge_FDC_Annual.csv` - Storage-discharge and FDC metrics
- `DS_drawdown_annual.csv` - Annual dynamic storage drawdown
- `Annual_GW_Prop.csv` - Annual mean baseflow proportion
- `stream_thermal_lowflow_metrics_annual.csv` - Max 7-day temp, min 7-day Q, temp at min Q, Q5_CV
- `HJA_StorageMetrics_Annual_All.csv` - All metrics, annual
- `HJA_Ave_StorageMetrics_CatCharacter.csv` - Site-averaged metrics

### Phase 3: Statistical Analyses
**Scripts:** `Statistical_analyses/`

Run in order:
1. `01_Correlations.R` - Site-averaged metrics + correlation matrices
2. `02_PCA.R` - PCA on storage metrics
3. `03_ANOVA_Tukey.R` - Test for site differences in storage metrics
4. `04_Catchment_Storage_MLR.R` - MLR predicting storage from catchment attributes
5. `05_RDA.R` - Redundancy analysis and variance partitioning
6. `06_Storage_Predicts_Thermal_LowFlow.R` - **KEY ANALYSIS** - Test storage → ecological response

**Outputs:**
- `HJA_Ave_StorageMetrics_CatCharacter.csv` - Site-averaged data
- `ANOVA_results.csv` - ANOVA F-statistics and p-values
- `Tukey_HSD_results.csv` - Pairwise site comparisons
- `MLR_Storage_Catchment_Results.csv` - Catchment → storage models
- `RDA_variance_explained.csv` - Variance partitioning results
- `Storage_Thermal_LowFlow_Models.csv` - **Storage → thermal/low-flow models**
- QA plots for all analyses

---

## Key Analyses for Manuscript

### 1. Storage Metrics (Hydrometrics)
- **RBI**: Flashiness index
- **Recession slope**: Drainage rate
- **Q5**: Low-flow magnitude
- **CV(Q5)**: Low-flow variability
- **Mean baseflow**: Chemical hydrograph separation
- **FDC slope**: Flow duration curve slope
- **S (annual)**: Annual storage (Kirchner-Staudinger)
- **ΔS (sum)**: Cumulative storage change

### 2. Ecological Response Variables (WY 1997-2020)
- **Max 7-day avg stream temp (°C)**: Indicator of thermal stress
- **Min 7-day avg discharge (mm/d)**: Indicator of drought stress
- **Temp at min Q (°C)**: Combined thermal-hydrologic stress

### 3. Hypotheses Tested
**H1:** Greater storage capacity → Lower maximum stream temperatures
**H2:** Greater storage capacity → Higher minimum discharge
**H3:** Greater storage capacity → Lower temperature at minimum discharge timing

---

## Data Requirements

### Inputs
All data reside in: `/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/03_Data/`

- **Discharge:** `Q/HF00402_v14.csv`
- **Stream Temperature:** `Stream_T/HT00201_results_*.csv`
- **Specific Conductance:** `EC/CF01201_v3.txt`
- **Catchment Attributes:** `DynamicStorage/Catchment_Charc.csv`
- **Drainage Areas:** `Q/drainage_area.csv`

### Outputs
All outputs saved to: `/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/`

---

## Data Coverage & Timeline

**Water Years:** 1997-2020 (Oct 1, 1996 - Sep 30, 2020)

**Sites:**
- GSLOOK (Lookout Creek)
- GSWS01, GSWS02, GSWS03 (WS01-03)
- GSWS06, GSWS07, GSWS08, GSWS09, GSWS10 (WS06-10)
- GSWSMC (Mack Creek)

**Data availability varies by site and metric** - see individual scripts for details.

---

## Software Requirements

**R ≥ 4.3.x**

### Core Packages
- Data manipulation: `dplyr`, `tidyr`, `readr`, `lubridate`
- Statistics: `MASS` (stepAIC), `car` (VIF), `vegan` (RDA)
- Visualization: `ggplot2`, `patchwork`, `ggcorrplot`, `ggrepel`, `GGally`
- Hydrology: `EflowStats`, `zoo`

---

## Notes

- **Paths:** All scripts use Box paths. Update `base_dir` and `output_dir` if running locally.
- **.gitignore:** Excludes `.DS_Store`, `.Rhistory`, `.RData`, and other R temp files.
- **Message statements:** Removed from all scripts for cleaner execution.
- **Water year definition:** Oct-Dec belongs to next calendar year (e.g., Oct 1996 = WY 1997).

---
