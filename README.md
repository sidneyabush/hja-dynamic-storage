# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews watersheds.

## run

1. Set paths in `config.R` (`USE_LOCAL_DATA`, Box/local directories).
2. Run the full workflow:
   - Terminal: `Rscript run_all.R`
   - Or set `HJA_REPO_DIR` if needed: `HJA_REPO_DIR="/path/to/hja-dynamic-storage" Rscript run_all.R`

## what it runs

- `01_metrics`: build master datasets and metric calculations
- `02_stats`: ANOVA/Tukey, PCA, MLR
- `03_plots`: manuscript and supplemental figures
- `verify_outputs.R`: required output checks

## mlr setup

- MLR is strict-only.
- Predictors are centered/scaled.
- Stepwise selection uses backward AIC/AICc workflow.
- Iterative VIF filtering is applied until all retained predictors are `<= 10`.
- Known high-correlation predictor pairs are constrained in final models.
- LOOCV statistics are exported.

## output location

By default (Box), outputs go to:

`/Users/sidneybush/Library/CloudStorage/Box-Box/05_Storage_Manuscript/05_Outputs/final_workflow`

Main subfolders:

- `master`
- `metrics`
- `stats`
- `tables`
- `figs`
