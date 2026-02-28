# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews catchments.

Eco-response variable notes:
- `WS09` has no stream-temperature records in `HT00451_v10.txt`, so `T_7DMax` and `T_Q7Q5` are unavailable (`NA`) for `WS09`.

CHS (EC baseflow) notes:
- CHS keeps a water year only when there are at least `300` days with overlapping EC + discharge data (`CHS_MIN_DAYS_PER_WY` in `config.R`).
- EC input file is `inputs/CF01201_v4.txt`.
- Meteorological input files are in `inputs/all_hydromet/`.

## scripts

- `00_data_preprocessing`: hydromet and ET preprocessing
- `01_storage_calcs`: master datasets and metric calculations
- `02_analysis`: ANOVA/Tukey, PCA, MLR, and manuscript table assembly
- `03_plots`: manuscript and supplemental figures

Main output subfolders:

- `master`
- `metrics`
- `models`
- `exploratory_plots`

Manuscript-ready exports:

- `ms_materials/main`
- `ms_materials/supp`
