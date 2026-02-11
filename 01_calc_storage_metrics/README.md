# Calculate Storage Metrics

**Storage metrics calculated:**
- **Dynamic (DS)** — RBI (flashiness), recession curve slope, flow duration curve slope, storage-discharge
- **Mobile (MS)** — chemical hydrograph separation (EC-based), mean transit time, young water fraction, damping ratio
- **Extended dynamic (EDS)** — water balance drawdown

# Script Overview:
- 01_DS_RBI_RCS - calculates RBI (flashiness) and RCS (recession curve slope)
- 02_DS_FDC_SD - calculates storage-discharge relationship and flow duration curve slope
- 03a_MS_load_isotope_metrics - loads isotope-derived mobile storage metrics -- here is where we can determine how we want to combine metrics (i.e., based on time or site)
- 04_EDS_water_balance - calculates water balance drawdown

**Sites:** 10 hydrometric watersheds (GSWS01–10, GSLOOK, GSWSMC) + 6 isotope-only sites

### Storage Metrics by Type

| Storage Type | Abbreviation | n | Sites | Years |
|--------------|--------------|---|-------|-------|
| **Dynamic** | RBI | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1950-2020 (varies) |
| | RCS | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1950-2020 (varies) |
| | FDC | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2020 |
| | SD | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2020 |
| **Mobile** | MTT | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS07, GSWS08, GSWSMC, MR | Site-level |
| | Fyw | 7 | GSWS01, GSWS02, GSWS07, GSWS08, GSWSMC, GSLOOK, MR | Site-level |
| | CHS | 8 | GSWS10, GSWS01, GSWS02, GSWS03, GSWS06*, GSWS07, GSWS08, GSWSMC | 2013-2019 (varies) |
| | DR | 15 | All hydrometric + MR, NC, LC, LO2, CC, LO1 | Site-level |
| **Extended Dynamic** | WB | 10 | GSWS09, GSWS10, GSWS01, GSLOOK, GSWS02, GSWS03, GSWS06, GSWS07, GSWS08, GSWSMC | 1998-2019 |

*GSWS06 has limited CHS data (2017-2019); CHS excludes GSWS09 and GSLOOK (no EC data)

---