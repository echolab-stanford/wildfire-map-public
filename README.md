# wildfire-map-public
Repo supporting "The changing risk and burden of wildfire in the United States"


# Using results

The estimates generated in the paper can be found in the folder `clean/results_all.RDS`. This file is a dataframe, with an estimate for each year, for each grid cell. The column `preds` is the prediction for the cell overall, and `preds_0` is the prediction when the smoke value is artificially set to 0. The grid cells are only identified by id, and can be merged with a shapefile for any spatial work. The grid shapefile can be found at `clean/national_grid.RDS`.  

There are several settings that can be changed in the `work/00_functions.R`. 

# Packages needed

A custom package used in the lab that includes several plotting and scraping functions is used throughout, and can be downloaded at the github `burkelab/census.tools`. 

Additionally, some of the scraping requires RSelenium, along with the dev version of chrome.


# Data

Several datasets are included in the repo in a processed form, and the raw data can be downloaded at that link provided below. For example:
1) the TIGER line files are included in a simplified form (same polygons, just passed through gSimplify so that they're small enough to include in GitHub) 
2) the EPA PM2.5 data is included in one big rds file, it can only be downloaded on a state-year basis (so the raw data is ±600 separate files)
3) the smoke polygon data is included in one big rds file, it can only be downloaded on a daily basis (4000+ separate files) 

Additionally, some datasets are not included at all, since they require log in to access and cannot be distributed. Those files are listed below.

## Not included in repo

boundaries/uszips.csv: roprietary data from simplemaps. http://simplemaps.com/data/updates/order/FLA201009-8749-21137/
improve: Sourced from http://vista.cira.colostate.edu/Improve/improve-data/, requires log in to access
pop: sourced from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download, requires log in to download
traffic-darte: Can be found at https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1735, requires log in to download

## Included in repo

airport: from https://www.transtats.bts.gov/ONTIME/
boundaries/GACC: boundaries are from https://hub.arcgis.com/datasets/nifc::national-gacc-boundaries
boundaries/tl_2019_us_state: state boundaries are TIGER line files https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
boundaries/all_national_zips.rds: simplified version of TIGER line files
boundaries/counties.RDS: simplified version of TIGER line files
boundaries/state_fips_codes.csv: manually created
census: all public access census and ag census files
coal: all sourced from https://www.eia.gov/electricity/data/eia860m/
emissions/abatzoglou_data.csv: Dataset S01 from https://www.pnas.org/content/suppl/2016/10/06/1607171113.DCSupplemental
emissions/supression_costs.csv: manually annotated from https://www.nifc.gov/fireInfo/fireInfo_documents/SuppCosts.pdf
emissions/USAGDPDEFAISMEI.csv: https://fred.stlouisfed.org/series/USAGDPDEFAISMEI
EPA_trend: each region downloaded from https://www.epa.gov/air-trends/particulate-matter-pm25-trends
fire/hms_fires: adapted from data downloaded (individually by day) from https://www.ospo.noaa.gov/Products/land/hms.html
fire/prescribed_burn_acres.csv: Manually annotated from Wildland Fire Summaries at https://www.nifc.gov/fireInfo/fireInfo_statistics.html
natural_gas: field level natural gas estimates, need to be aggregated to county level. https://www.eia.gov/naturalgas/ngqs/#?report=RP7&year1=2005&year2=2019&company=Name
physio_shp: from https://water.usgs.gov/GIS/dsdl/physio_shp.zip
pm: sourced from https://www.epa.gov/outdoor-air-quality-data/download-daily-data and slightly processed. Processing file is provided in  `work/supplemental` in case you would like to update the data in the future.
powerplants/emissions_[YYYY].csv: from the data browser, you need to manually select the year of interest and download from the button above the table https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2009&tab=annual_emissions
powerplants/overview_[YYYY].csv: same as above https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2006&tab=overview&start=200101&end=201710 
powerplants/Plant_Y[YYYY].xlsx: https://www.eia.gov/electricity/data/eia860/
smoke: adapted from data downloaded (individually by day) from https://www.ospo.noaa.gov/Products/land/hms.html. Processing file is provided in `work/supplemental` in case you would like to update the data in the future.
WUI: state level estimates of number of homes in the wildland urban interface. Generated using National Land Cover Database and proprietary CoreLogic data including the locations of all homes in the US.

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
