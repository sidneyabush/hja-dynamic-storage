# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews watersheds.

## run

1. Set paths in `config.R` (`USE_LOCAL_DATA`, Box/local directories).
2. Run the full workflow:
   - Terminal: `Rscript run_all.R`
   - Or set `HJA_REPO_DIR` if needed: `HJA_REPO_DIR="/path/to/hja-dynamic-storage" Rscript run_all.R`

## scripts

- `01_metrics`: build master datasets and metric calculations
- `02_stats`: ANOVA/Tukey, PCA, MLR
- `03_plots`: manuscript and supplemental figures
- `verify_outputs.R`: required output checks


Main output subfolders:

- `master`
- `metrics`
- `stats`
- `tables`
- `figs`
