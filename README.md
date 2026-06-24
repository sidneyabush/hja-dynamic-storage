# A Unified Framework for Dynamic and Mobile Storage Across Volcanic Montane Catchments

Bush, S. A., Johnson, K., Segura, C., Perry, Z., Duffy, S., and Sullivan, P. L. (submitted). A Unified Framework for Dynamic and Mobile Storage Across Volcanic Montane Catchments. *Water Resources Research*.

## Abstract

Subsurface storage regulates streamflow variability, sustains low flows, and moderates stream temperatures, but is difficult to measure directly and is therefore typically inferred from stream flow and chemistry data. We analyzed long-term hydrometric, chemical, and isotopic records from Lookout Creek and nine subcatchments in the H.J. Andrews Experimental Forest, Oregon, USA, to evaluate how commonly used storage metrics covary and relate to ecologically relevant hydrologic and thermal responses. We developed a unified framework combining dynamic storage inferred from streamflow, extended-dynamic storage represented by catchment water-balance deficit, and mobile storage characterized by calcium-derived baseflow fraction and isotope-derived mixing and transport metrics. Across this heterogeneous volcanic landscape, storage metrics varied among catchments but did not produce a single consistent high-to-low storage ranking of catchments. Relationships among storage metrics followed several clear patterns but varied in strength: isotope-derived mobile metrics were most strongly related to flashiness and storage-discharge behavior, whereas baseflow fraction was higher in catchments with greater storage-discharge estimates and lower in flashier catchments. Geology, slope, and landslide extent further organized variation among catchments. Storage metrics retained in predictive models differed between stream low-flow and thermal responses. Annual seven-day low-flow discharge was most strongly associated with dynamic and extended-dynamic storage metrics and wet-season precipitation, whereas annual maximum seven-day stream temperature retained both flow-based and tracer-based storage metrics, indicating that thermal responses reflect both discharge regulation and source-water mixing. These results support a multi-metric framework that interprets subsurface storage using complementary storage indicators, rather than any single metric.

## Quick Start

From the repository folder, run:

```r
Rscript install_packages.R
Rscript run_all.R
```

`install_packages.R` is a one-time setup step. `run_all.R` runs the full workflow.

The workflow reads inputs from `HJA_FINAL_WORKFLOW_ROOT/inputs` and writes outputs to `HJA_FINAL_WORKFLOW_ROOT/outputs` and `HJA_FINAL_WORKFLOW_ROOT/figs_tables_pub`. Set `HJA_FINAL_WORKFLOW_ROOT` to the folder that stores the workflow data and publication figure and table files.

## Software

- `R 4.5.2`
- Packages are listed in `install_packages.R`
- Software citation metadata are provided in `CITATION.cff` and `.zenodo.json`
- License: MIT

## Code And Data Availability

The R code used to calculate storage metrics, run the analyses, and produce the manuscript figures and tables is developed at <https://github.com/sidneyabush/hja-dynamic-storage> and should be cited using the archived Zenodo software DOI after release. Data used by the workflow are publicly available from the sources cited below. The gap-filled meteorological inputs used by the workflow are derived from public H.J. Andrews Forest records.

## Data References

- Bell, D. M. (2023). Time-series of lidar data for the H. J. Andrews Experimental Forest and Willamette National Forest, 2008-2021 [Data set]. Long-Term Ecological Research; Forest Science Data Bank. Retrieved June 18, 2026, from <https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=GI013>
- Daly, C., Schulze, M. D., & McKee, W. A. (2025). Air temperature, relative humidity, dewpoint temperature, water vapor pressure deficit, and atmospheric pressure data from benchmark stations at the HJ Andrews Experimental Forest, 1957 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/5d4ab4b210165d6e860ebe58e0579e4e>
- H.J. Andrews Experimental Forest. (2019, May 27). Harvest sites [Data set]. H.J. Andrews Experimental Forest and LTER Open Data Hub. Retrieved June 18, 2026, from <https://data-osugisci.opendata.arcgis.com/datasets/2d8986828d3d4d07af77784bcea3845f/explore?layer=0>
- Johnson, S. L., Gregory, S., Henshaw, D., & Kennedy, A. (2026). Stream and air temperature data from stream gages and stream confluences in the Andrews Experimental Forest, 1950 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/4d0e647de6e461c29ff13456849ea328>
- Johnson, S. L., Henshaw, D., Nash, B., Remillard, S., & Rothacher, J. (2025). Stream discharge in gaged watersheds at the HJ Andrews Experimental Forest, 1949 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/9826b79c5bb9e80e37c08c02b2ee13f6>
- Johnson, S. L., Henshaw, D., Remillard, S., & Fredriksen, R. L. (2026). Stream chemistry concentrations and fluxes using proportional sampling in the Andrews Experimental Forest, 1968 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/42cb9e9d33248e56386454c89cdba777>
- Schulze, M. D., Daly, C., & McKee, W. A. (2026). Solar radiation data from benchmark stations at the HJ Andrews Experimental Forest, 1973 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/6183856f5450ad6be6745de8581ced63>
- Schulze, M. D., Daly, C., & Rothacher, J. (2026). Precipitation measurements from historic and current standard, storage and recording rain gauges at the Andrews Experimental Forest, 1951 to present [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/c30028e97bf7ec5be4fb72d08ab64bd2>
- Segura, C., Perry, Z., & Ortega, J. (2025). Water stable isotopes for streams and precipitation samples in the HJ Andrews Experimental Forest, 2014-2023 [Data set]. Long-Term Ecological Research; Forest Science Data Bank. <https://doi.org/10.6073/pasta/460268c31591c065921bcb7426d987cc>
