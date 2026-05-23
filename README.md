# hja_dynamic_storage

This repository reruns the analysis used in the paper on dynamic, extended-dynamic, and mobile storage in H.J. Andrews catchments

It recreates the paper figures, tables, and supporting outputs from the prepared inputs in `inputs/`. It does not rebuild every dataset from raw Andrews downloads, and it is not meant to be a general hydrology package

To run it, put the required files in `inputs/`, install packages with `Rscript install_packages.R`, and run `Rscript run_all.R`

## Software

- Tested with `R 4.5.2`
- Required packages are listed in [install_packages.R](install_packages.R)

## What This Repository Produces

Main-text figures:

- Figure 2: dynamic and extended-dynamic storage summaries
- Figure 3: PCA of dynamic and extended-dynamic metrics
- Figure 4: mobile storage summaries
- Figure 5: dynamic-mobile correlation matrix
- Figure 6: catchment controls on storage metrics
- Figure 7: storage controls on ecological response variables
- Figure 8: observed versus predicted ecological responses
- Figure 9: geology-landslide-storage framework

Main-text tables:

- Table 4
- Table 5

Supplement items written by the code:

- Figure S1
- Table S5
- Table S6
- Table S7
- Table S8

The code does not write Tables S1-S4 or Text S1-S4 as separate files. Those stay outside the code-generated outputs

## Before You Run

This code uses two kinds of inputs

Public Andrews files:

- `MS001`: benchmark and secondary meteorology
- `MS004`: precipitation network
- `MS050`: solar radiation / net radiation support used in ET processing
- `HF004`: stream discharge
- `HT004`: stream temperature
- `CF012`: stream specific conductance
- `CF002`: stream chemistry

Prepared paper inputs:

- gap-filled temperature, precipitation, and SWE files used in the paper
- `catchment_char.csv`
- `drainage_area.csv`
- `MTT_FYW.csv`
- `DampingRatios_2025-07-07.csv`

Those prepared files are part of the reproducible analysis. This repository does not rebuild every input from raw public downloads alone. If you post this code on Zenodo, those files should go with it or be linked clearly

The exact filenames expected by the code are listed in [inputs/README.md](inputs/README.md)

The code will stop if any required file is missing or misnamed

## Where To Get The Public Data

Start with the Andrews Forest data catalog:

- <https://andrewsforest.oregonstate.edu/data>

Study pages used here:

- `MS001` meteorology: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=MS001>
- `MS004` precipitation: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=MS004>
- `MS050` solar radiation: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=MS050>
- `HF004` discharge: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=HF004>
- `HT004` stream temperature: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=HT004>
- `CF012` specific conductance: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=CF012>
- `CF002` stream chemistry: <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=CF002>

## Default Folder Layout

By default, the repository looks for inputs and writes outputs in this simple local layout:

```text
hja_dynamic_storage/
├── inputs/
│   └── ...
├── outputs/
├── ms_materials/
│   ├── main/
│   └── supp/
├── run_all.R
└── install_packages.R
```

You only need to provide `inputs/`

`outputs/` and `ms_materials/` are created automatically when the code runs

If you want different folders, set:

- `HJA_REPO_DIR`
- `HJA_FINAL_WORKFLOW_ROOT`
- `HJA_BASE_DATA_DIR`
- `HJA_OUTPUT_DIR`
- `HJA_MS_MATERIALS_DIR`

## What `run_all.R` Does

In order, the full run:

1. Checks that the required input files are present.
2. Builds catchment-scale daily meteorological forcing and ET inputs.
3. Calculates annual storage metrics and ecological response metrics.
4. Runs ANOVA, PCA, regression models, and the MTT sensitivity analysis.
5. Makes the main-text figures and supplement figures.
6. Writes the main-text and supplement tables.
7. Verifies that the expected outputs were written.

## Run The Code

1. Put the required input files in `inputs/`.
2. Install packages if needed:

```bash
Rscript install_packages.R
```

3. Optional input check:

```bash
Rscript helpers/check_inputs.R
```

4. Run the code:

```bash
Rscript run_all.R
```

5. Optional output check:

```bash
Rscript helpers/verify_outputs.R
```

`run_all.R` already calls both checks

## Repository Layout

- `00_data_preprocessing/`: builds catchment-scale daily meteorological and ET inputs
- `01_storage_calcs/`: calculates annual storage and ecological response metrics
- `02_analysis/`: runs ANOVA, PCA, regression models, and the MTT sensitivity analysis
- `03_plots/`: makes the paper figures and the supplement figures
- `04_tables/`: writes the main-text and supplement tables
- `helpers/`: input checks, unit checks, helper functions, and final output checks
- `legacy_code/`: older development and comparison scripts kept for reference, but not required by `run_all.R`

## Files Written By The Code

- `outputs/metrics/`: annual and site-level hydrologic metrics
- `outputs/models/`: ANOVA, PCA, regression, and sensitivity-analysis outputs
- `outputs/master/`: combined annual and site summary tables used later
- `ms_materials/main/`: paper-ready main-text figures and tables
- `ms_materials/supp/`: supplement figure and tables made by the code

The MTT sensitivity step also writes a fuller set of outputs to `outputs/MTT_sensitivity/` and copies `TableS5_MTT_sensitivity.csv` into `ms_materials/supp/`

Main-text figures are written as both `.png` and `.pdf`. Supplement figure exports also include `.png` and `.pdf` where needed

## What Counts As A Successful Run

- `Rscript run_all.R` exits without error
- `outputs/master/master_annual.csv` and `outputs/master/master_site.csv` are created
- `ms_materials/main/` contains Figures 2-9
- `ms_materials/supp/` contains `FigS1` and `TableS5` through `TableS8`
- `Rscript helpers/verify_outputs.R` returns without error

## A Few Notes

- The analysis period is water years 1997-2020.
- `WB` is treated as a maximum within-year water-balance deficit, not as absolute storage.
- `DR`, `Fyw`, and `MTT` are used as site-level isotope metrics in this run
- `WS09` has no usable stream-temperature record in `HT00451_v10.txt`, so `T7DMax` and `Q7Q5`-based thermal metrics are not available there.
- The MTT sensitivity step is the slowest part of the run and may sit quiet for a while before finishing
- `legacy_code/` is not required to rerun the paper results
