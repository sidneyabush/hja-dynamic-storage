# Input Files

This folder is intentionally empty in the code archive

To rerun the paper, place the required files here using the filenames below. `helpers/check_inputs.R` looks for these names directly

After the files are in place, `Rscript helpers/check_inputs.R` gives a quick check before the full run

## Public Study Exports

Download these from the official Andrews Forest study pages listed in the main [README](../README.md):

- `all_hydromet/MS00102_v9.csv`
- `all_hydromet/MS05025_v3.csv`
- `all_hydromet/MS00403_v2.csv`
- `HF00402_v14.csv`
- `HT00451_v10.txt`
- `CF01201_v4.txt`
- `CF00201_v7.csv`

## Prepared Input Files

These are the prepared files used in the paper. The code assumes they already exist in this form. For a Zenodo copy, these files should be archived with the code or linked clearly:

- `all_hydromet/Temperature_original_&_filled_1979_2023_v2.csv`
- `all_hydromet/Precipitation_original_&_filled_1979_2023.csv`
- `all_hydromet/SWE_original_&_filled_1997_2023_v5.csv`
- `catchment_char.csv`
- `drainage_area.csv`
- `MTT_FYW.csv`
- `DampingRatios_2025-07-07.csv`

## Expected Tree

```text
inputs/
├── all_hydromet/
│   ├── MS00102_v9.csv
│   ├── MS00403_v2.csv
│   ├── MS05025_v3.csv
│   ├── Precipitation_original_&_filled_1979_2023.csv
│   ├── SWE_original_&_filled_1997_2023_v5.csv
│   └── Temperature_original_&_filled_1979_2023_v2.csv
├── catchment_char.csv
├── drainage_area.csv
├── CF00201_v7.csv
├── CF01201_v4.txt
├── DampingRatios_2025-07-07.csv
├── HF00402_v14.csv
├── HT00451_v10.txt
└── MTT_FYW.csv
```
