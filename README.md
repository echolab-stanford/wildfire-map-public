# wildfire-map-public
Repo supporting "The changing risk and burden of wildfire in the United States"


# Using results

The main estimates generated in the paper are available at `clean/results_all.RDS`. This file is a dataframe, with an estimate for each year, for each grid cell. The column `preds` is the prediction for the cell overall, and `preds_0` is the prediction when the smoke input value is artificially set to 0. For spatial coordinates, the data can be merged with the grid shapefile `clean/national_grid.RDS` on the variable `id`.  

There are several settings that can be changed in the `work/00_functions.R`. 

# Scripts

Description of script structuring, etc.

# Issues, Constraints, Choices

* Lots of file paths are used that are Mac specific
* Code currently uses the velox package for many raster extract operations, which isn't compatible with R 4.0.0+.
* RSelenium requires some pretty specific setup, including downloading the dev version of chrome. This is only reqired for downloading the airport data, which is included in the repo and could be updated manually for new years of data.


## Niche Issues

* download for the airport data will fail if you have any files in your download folder that fulfill the regex "\_T\_.+zip$"

# Data

Input data from freely distributed sources are included in the repo in a processed form. Directions for accessing the raw data are included along with the download links provided below. Examples of the processed data format include:
1) the TIGER line files are included in a simplified form (same polygons, just passed through gSimplify so that they're small enough to include in a GitHub repo) 
2) the EPA PM2.5 data are included in one big rds file whereas it can only be downloaded on a state-year basis (so the raw data is ±600 separate files)
3) the smoke polygon data are included in one big rds file, whereas they can only be downloaded on a daily basis (4000+ separate files) 

The last two datasets in that list do fall into the cateogry of "freely distributed sources" but are too large for Github. Those two files can be found at: [XYZZY]

Additionally, some input datasets are not included here in any form because they require registration prior to access and thus cannot be redistributed. Those files, along with links for access, are listed below.

## Not included in repo

* data/boundaries/uszips.csv: proprietary data from simplemaps are available for [purchase](https://simplemaps.com/data/us-zips). The `Pro` version of the zips database is utilized for the analysis (filename.csv).
* data/improve: Data requires [registration](http://views.cira.colostate.edu/fed/Auth/Register.aspx) prior to [downloading](http://vista.cira.colostate.edu/Improve/improve-data/). We use files `IMPROVE_1988-2006.txt`, `IMPROVE_2007.txt`, `IMPROVE_2008-2016`, `IMPROVE_2017.txt`, and `IMPROVE_2018.txt`.
* data/pop: SEDAC population data requires [registration](https://sedac.ciesin.columbia.edu/user-registration) prior to [downloading](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download). We use the 2.5 minute data for 2005, 2010, 2015, and 2020. Files are of the form: `gpw_v4_population_count_rev11_[Year]_2pt5_min.tif`.
* data/traffic-darte: DARTE data requires [registration](https://urs.earthdata.nasa.gov/users/new?client_id=YQOhivHfMTau88rjbMOVyg&redirect_uri=https%3A%2F%2Fdaac.ornl.gov%2Fcgi-bin%2Furs%2Furs_logon_proc.pl&response_type=code&state=https%3A%2F%2Fdaac.ornl.gov%2Fcgi-bin%2Fdsviewer.pl%3Fds_id%3D1735) prior to [downloading](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1735) the files `onroad_[YYYY].tif` for years 2006-2017.

## Included in repo

* data/airport: Airport data come from [here](https://www.transtats.bts.gov/ONTIME/) including airport locations, quarterly ticketing, and monthly on-time status.

* data/boundaries/GACC: The [GACC boundaries](https://hub.arcgis.com/datasets/nifc::national-gacc-boundaries) can be downloaded [here](https://opendata.arcgis.com/datasets/7dc5f4a286bd47e0aaafa0ab05302fe9_0.gdb)

* data/boundaries/tl_2019_us_state: [state boundaries](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html) are TIGER line files and can be downloaded [here](https://www2.census.gov/geo/tiger/TIGER2019/STATE/tl_2019_us_state.zip)

* data/boundaries/all_national_zips.rds: simplified version of TIGER line files

* data/boundaries/counties.RDS: simplified version of TIGER line files

* data/boundaries/state_fips_codes.csv: manually created

* data/census: all public access census and ag census files

* data/coal: all sourced from https://www.eia.gov/electricity/data/eia860m/

* data/emissions/abatzoglou_data.csv: Dataset S01 from https://www.pnas.org/content/suppl/2016/10/06/1607171113.DCSupplemental

* data/emissions/supression_costs.csv: manually annotated from https://www.nifc.gov/fireInfo/fireInfo_documents/SuppCosts.pdf

* data/emissions/USAGDPDEFAISMEI.csv: https://fred.stlouisfed.org/series/USAGDPDEFAISMEI

* data/EPA_trend: each region downloaded from https://www.epa.gov/air-trends/particulate-matter-pm25-trends

* data/fire/hms_fires: adapted from data downloaded (individually by day) from https://www.ospo.noaa.gov/Products/land/hms.html

* data/fire/prescribed_burn_acres.csv: Manually annotated from Wildland Fire Summaries at https://www.nifc.gov/fireInfo/fireInfo_statistics.html

* data/natural_gas: field level natural gas estimates, need to be aggregated to county level. https://www.eia.gov/naturalgas/ngqs/#?report=RP7&year1=2005&year2=2019&company=Name

* data/physio_shp: from https://water.usgs.gov/GIS/dsdl/physio_shp.zip

* data/pm: sourced from https://www.epa.gov/outdoor-air-quality-data/download-daily-data and slightly processed. Processing file is provided in  `work/supplemental` in case you would like to update the data in the future.

* data/powerplants/emissions_[YYYY].csv: from the data browser, you need to manually select the year of interest and download from the button above the table https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2009&tab=annual_emissions

* data/powerplants/overview_[YYYY].csv: same as above https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2006&tab=overview&start=200101&end=201710 

* data/powerplants/Plant_Y[YYYY].xlsx: https://www.eia.gov/electricity/data/eia860/

* data/smoke: adapted from data downloaded (individually by day) from https://www.ospo.noaa.gov/Products/land/hms.html. Processing file is provided in `work/supplemental` in case you would like to update the data in the future.

* data/WUI: state level estimates of number of homes in the wildland urban interface. Generated using National Land Cover Database and proprietary CoreLogic data including the locations of all homes in the US.

### Necessary starting data folder structure to reproduce

```
data
 ├── airport
 │	 ├── airport_locations.csv
 │	 ├── Q [Quarter]_[Year]_tickets.csv (for 'Q 1' to 'Q 4' 2006-2018)
 │   └── [Month]_[Year]_ontime.csv (for Jan-Dec 2006-2018)
 ├── boundaries
 │   ├── GACC
 │   │	 └── National_GACC_Current_20200226.shp (and associated dbf, prj) 
 │	 ├── tl_2010_us_state
 │	 ├── all_national_zips.rds
 │	 ├── counties.RDS
 │	 ├── state_fips_codes.csv
 │	 └── uszips.csv
 ├── census
 │	 ├── acres_conventional_tillage.csv 
 │	 ├── ACS_17_1YR_S0101.csv
 │	 ├── BP_[Year]_00A1_with_ann.csv (for 2006-2016)
 │	 ├── CAINC6N__ALL_STATES_2001_2017.csv
 │	 ├── cattle_number.csv
 │	 ├── county_ag_sales.csv
 │	 └── fertilizer_totals.csv
 ├── coal
 │	 ├── existing_gen_units_2006.xls (2006 - 2014)
 │	 └── existing_gen_units_2015.xlsx (2015 - 2018)
 ├── emissions
 │	 ├── abatzoglou_data.csv
 │	 ├── suppression_costs.csv 
 │	 └── USAGDPDEFAISMEI.csv
 ├── EPA_trend
 │	 └── PM25[Region].csv (for Central, National, Northeast, NorthernRockies, Northwest, 
 │						   South, Southeast, Southwest, UpperMidwest, West)
 ├── fire
 │	 ├── hms_fires.RDS
 │	 └── prescribed_burn_acres.csv
 ├── improve
 │	 ├── IMPROVE_1998-2006.txt
 │	 ├── IMPROVE_2007.txt
 │	 ├── IMPROVE_2008-2016.txt
 │	 ├── IMPROVE_2017.txt
 │	 └── IMPROVE_2018.txt
 ├── model_inputs
 ├── natural_gas
 │ 	 └── county_level_natural_gas_estimates_2006_2016.rds
 ├── physio_shp
 │ 	 └── physio (shapefile of physiographic divisions) 
 ├── pm
 │ 	 └── epa_station_level_pm25_data.rds
 ├── pop
 │	 └── gpw_v4_population_count_rev11_[Year]_2pt5_min.tif (2005,  2010, 2015, 2020)
 ├── powerplants
 │	 ├── Plant_Y[Year].xls (2006 - 2010)
 │	 ├── Plant_Y[Year].xlsx (2011 -  2018)
 │	 ├── overview_[Year].csv (2006 - 2018)
 │	 └── emissions_[Year].csv  (2006 - 2018)
 ├── smoke
 │	 └──  smoke_plumes.rds
 ├── traffic-darte
 │	 └── onroad_[Year].tif (2006 - 2017)
 └── WUI
 	 └── [Year]_combined_wui_hh_data.csv (2001, 2004, 2006, 2008, 2011, 2013, 2016)
```
# R Packages needed

R packages required for replications are:

- data.table
- dplyr
- gdata 
- geosphere
- imputeTS
- ncdf4
- openxlsx 
- raster
- readr
- rgdal
- rgeos
- RSelenium
- sf
- signal 
- sp
- splines
- stringr
- tidyr
- velox

Users can run the following one-off command to install the most recent versions of these packages:
```
install.packages(c('BAMMtools','data.table','devtools','dplyr','gdata', 'geosphere','Hmisc','imputeTS','ncdf4','openxlsx', 'raster','readr','rgdal','rgeos','RSelenium','sf','signal', 'sp','stringr','tidyr','velox'), dependencies = T)
```

Finally, some of the scraping requires the [dev version](https://www.google.com/chrome/dev/) of Google Chrome.


All scripts were written in R version 3.6.1 (2019-07-05) -- "Action of the Toes".

