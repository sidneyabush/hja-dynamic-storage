# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews watersheds.

Eco-response variable notes:
- `WS09` has no stream-temperature records in `HT00451_v10.txt`, so `T_7DMax` and `T_Q7Q5` are unavailable (`NA`) for `WS09`.

CHS (EC baseflow) notes:
- CHS keeps a water year only when there are at least `300` days with overlapping EC + discharge data (`CHS_MIN_DAYS_PER_WY` in `config.R`).
- EC input file is `inputs/CF01201_v4.txt`.
- Meteorological input files are in `inputs/all_hydromet/`.

## scripts

- `01_metrics`: build master datasets and metric calculations
- `02_stats`: ANOVA/Tukey, PCA, MLR
- `03_plots`: manuscript and supplemental figures

Main output subfolders:

- `master`
- `metrics`
- `models`
- `tables`
- `figs`
