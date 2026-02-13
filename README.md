# hja_dynamic_storage

Dynamic, mobile, and extended dynamic storage workflow for H.J. Andrews watersheds.

Eco-response notes:
- Stream-temperature responses are computed from `final_workflow/inputs/stream_t/HT00451_v10.txt` (`WATERTEMP_MEAN`).
- Eco responses are computed on water years (`WY`) `1997`-`2020`.
- `WS09` has no stream-temperature records in `HT00451_v10.txt`, so `T_7DMax` and `T_Q7Q5` are unavailable (`NA`) for `WS09`.

## scripts

- `01_metrics`: build master datasets and metric calculations
- `02_stats`: ANOVA/Tukey, PCA, MLR
- `03_plots`: manuscript and supplemental figures
- `verify_outputs.R`: required output checks

Main output subfolders:

- `master`
- `metrics`
- `models`
- `tables`
- `figs`
