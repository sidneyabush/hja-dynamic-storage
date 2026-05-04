# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews catchments.

Eco-response variable notes:
- `WS09` has no stream-temperature records in `HT00451_v10.txt`, so `T_7DMax` and `T_Q7Q5` are unavailable (`NA`) for `WS09`.

BF notes:
- BF keeps a water year only when there are at least `300` days with overlapping EC + discharge data (`CHS_MIN_DAYS_PER_WY` in `config.R`).
- EC input file is `inputs/CF01201_v4.txt`.
- Ca/chemistry comparison uses `inputs/CF00201_v7.csv` with `CA` and `COND` and keeps years with at least `10` chemistry observations (`CHS_MIN_OBS_PER_WY_CHEM`).
- Main-workflow BF is Ca-based (`annual_gw_prop_ca.csv`).
- EC-vs-Ca supplemental comparison uses BF from daily EC (`annual_gw_prop.csv`) versus BF from Ca (`annual_gw_prop_ca.csv`).
- Comparison outputs are written to `outputs/metrics/mobile/annual_gw_prop_ec_ca_comparison.csv`, plus site-year paired values (`annual_gw_prop_ec_ca_site_year_pairs.csv`) and per-site stats (`annual_gw_prop_ec_ca_site_by_site_stats.csv`) in the same folder.
- Meteorological input files are in `inputs/all_hydromet/`.

Isotope notes:
- Final workflow uses site-level mean isotope metrics only (`MTT`, `Fyw`, `DR`).
- Annual isotope-ingestion paths are removed from the core workflow.

WB (extended dynamic storage) note:
- `WB` now uses a max within-year depletion workflow (`max(cummax(cumsum(P-Q-ET)) - cumsum(P-Q-ET))`), not a simple annual range; the legacy peak-anchored WB workflow is removed from the core pipeline.

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
