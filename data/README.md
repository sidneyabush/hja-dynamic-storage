# Data Files Required

## Data Sources

All data are available from the H.J. Andrews Experimental Forest Long-Term Ecological Research program:
- **Website**: https://andrewsforest.oregonstate.edu/data
- **Data Portal**: https://portal.edirepository.org/

## Required Files

### Q/ (Discharge Data)
- `HF00402_v14.csv` - Daily mean discharge for all watersheds
- `drainage_area.csv` - Drainage areas for each site

### DynamicStorage/ (Water Balance Data)
- `daily_water_balance_ET_Hamon-Zhang_coeff_interp.csv` - Daily P, Q, ET
- `Catchment_Charc.csv` - Catchment characteristics (slope, aspect, geology, etc.)

### EC/ (Specific Conductance)
- `CF01201_v3.txt` - Continuous specific conductance measurements

### Isotopes/ (Isotope Data)
- `MTT_FYW.csv` - Mean transit times and young water fractions
- `DampingRatios_2025-07-07.csv` - Isotopic damping ratios

### Stream_T/ (Stream Temperature)
- `HT002*.csv` - Stream temperature files (multiple files by site)

### MET/ (Meteorology)
- Various meteorological data files

## File Format Notes

- Most files are comma-separated (CSV) or tab-delimited (TXT)
- Date formats vary: `YYYY-MM-DD` or `MM/DD/YYYY`
- Site codes use HJA naming convention (e.g., GSWS01, GSLOOK, GSWSMC)
