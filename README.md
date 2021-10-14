# wildfire-map-public

Repo supporting ["The changing risk and burden of wildfire in the United States"](https://www-pnas-org.stanford.idm.oclc.org/content/118/2/e2011048118).


# Using results

The main estimates generated in the paper are available at `clean/results_all.RDS`. This file is a dataframe, with an estimate for each year, for each grid cell. The column `preds` is the prediction for the cell overall, and `preds_0` is the prediction when the smoke input value is artificially set to 0. For spatial coordinates, the data can be merged with the grid shapefile `clean/national_grid.RDS` on the variable `id`.  

There are several settings that can be changed in the `work/00_functions.R`. 

The entire repo is ~1.4 GB zipped and ~1.7 GB unzipped.


# Scripts

Description of script structuring, etc.


# Issues, Constraints, Choices

* Lots of file paths are used that are Mac specific
* Code currently uses the velox package for many raster extract operations, which isn't compatible with R 4.0.0+.
* RSelenium requires some pretty specific setup, including downloading the dev version of Chrome. This is only reqired for downloading the airport data, which is included in the repo and could be updated manually for new years of data.

## Niche Issues

* download for the airport data will fail if you have any files in your download folder that fulfill the regex "\_T\_.+zip$"


# Data

Input data from freely distributed sources are included in the repo in a processed form. Directions for accessing the raw data are included along with the download links provided below. Examples of the processed data format include:
1) the TIGER line files are included in a simplified form (same polygons, just passed through gSimplify so that they're small enough to include in a GitHub repo) 
2) the EPA PM2.5 data are provided in one big rds file whereas it can only be downloaded on a state-year basis (so the raw data is ±600 separate files)
3) the smoke polygon data are included in one big rds file, whereas they can only be downloaded on a daily basis (4000+ separate files) 

The last two datasets in that list do fall into the cateogry of "freely distributed sources" but are too large for Github. Those two files can be found at the links below.

Additionally, some input datasets are not included here in any form because they require registration prior to access and thus cannot be redistributed. Those files, along with links for access, are listed below.

## Not included in repo

* data/boundaries/tl_2019_us_county: [county boundaries](https://www2.census.gov/geo/tiger/TIGER2019/COUNTY/) are TIGER line files and are too large for github but can be downloaded [here](https://www2.census.gov/geo/tiger/TIGER2019/COUNTY/tl_2019_us_county.zip)

* data/fire/hms_fires.RDS: File too big for Github but can be downloaded [here](https://www.dropbox.com/s/hv9qb0l7no3ec19/hms_fires.RDS?dl=0). Adapted from data downloaded (individually by day) from [here](https://www.ospo.noaa.gov/Products/land/hms.html).

* data/improve: Data requires [registration](http://views.cira.colostate.edu/fed/Auth/Register.aspx) prior to [downloading](http://vista.cira.colostate.edu/Improve/improve-data/). We use the OCf parameter from the IMPROVE Aerosol daily data set for all sites 1988-2018. Files were downloaded by year groups and the data component of the files were extracted to create the files referenced in the script: `IMPROVE_1988-2006.txt`, `IMPROVE_2007.txt`, `IMPROVE_2008-2016`, `IMPROVE_2017.txt`, and `IMPROVE_2018.txt`.

* data/pm: File is too large for Github but can be downloaded [here](https://www.dropbox.com/s/uueqfjixp74fxh7/epa_station_level_pm25_data.rds?dl=0). Raw data come from the [EPA download portal](https://www.epa.gov/outdoor-air-quality-data/download-daily-data) and are slightly processed to create the file linked above using the script  `work/supplemental/create_simplified_epa_data.R` in case you would like to update the data in the future. 

* data/pop: SEDAC population data requires [registration](https://sedac.ciesin.columbia.edu/user-registration) prior to [downloading](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download). We use the 2.5 minute data for 2005, 2010, 2015, and 2020. Files are of the form: `gpw_v4_population_count_rev11_[Year]_2pt5_min.tif`.

* data/traffic-darte: DARTE data requires [registration](https://urs.earthdata.nasa.gov/users/new?client_id=YQOhivHfMTau88rjbMOVyg&redirect_uri=https%3A%2F%2Fdaac.ornl.gov%2Fcgi-bin%2Furs%2Furs_logon_proc.pl&response_type=code&state=https%3A%2F%2Fdaac.ornl.gov%2Fcgi-bin%2Fdsviewer.pl%3Fds_id%3D1735) prior to [downloading](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1735) the files `onroad_[YYYY].tif` for years 2006-2017.

* data/vanD/vanD-[Year].asc: Too large to include on Github, can be downloaded [here](http://fizz.phys.dal.ca/~atmos/martin/?page_id=140#V4.NA.03). Files have to be renamed as when downloaded they follow the naming convention `V4NA03_PM25_NA_[Year]01_[Year]12-RH35-NoNegs.asc.zip`


## Included in repo

* data/airport: Airport data come from [here](https://www.transtats.bts.gov/ONTIME/) including airport locations, quarterly ticketing, and monthly on-time status.

* data/boundaries/GACC: The [GACC boundaries](https://hub.arcgis.com/datasets/nifc::national-gacc-boundaries) can be downloaded [here](https://opendata.arcgis.com/datasets/7dc5f4a286bd47e0aaafa0ab05302fe9_0.gdb)

* data/boundaries/tl_2019_us_state: [state boundaries](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html) are TIGER line files and can be downloaded [here](https://www2.census.gov/geo/tiger/TIGER2019/STATE/tl_2019_us_state.zip)

* data/boundaries/all_national_zips.rds: simplified version of TIGER line files

* data/boundaries/counties.RDS: simplified version of TIGER line files

* data/boundaries/state_fips_codes.csv: manually created

* data/Brey: region specific estimates of smoke source from [Brey et al 2018](https://d-nb.info/1162355557/34) [repo](https://salix.atmos.colostate.edu/svn/smokeSource/).

* data/census: all public access census and ag census files

* data/coal: Annual data on coal generation. Data come from Form EIA-860-Schedule 3, "Generator Data" originally sourced from [here](https://www.eia.gov/electricity/data/eia860/). The files included in the repo are these annual files (`3_1_Generator[YEAR]`) subset to only Coal plants with the files renamed.

* data/emissions/abatzoglou_data.csv: Dataset S01 from [Abatzoglou and Williams 2016](https://www.pnas.org/content/113/42/11770) was downloaded from the supplement [here](https://www.pnas.org/highwire/filestream/623620/field_highwire_adjunct_files/0/pnas.1607171113.sd01.csv)

* data/emissions/supression_costs.csv: manually annotated from [reported NIFC numbers](https://www.nifc.gov/fireInfo/fireInfo_documents/SuppCosts.pdf). 

* data/emissions/USAGDPDEFAISMEI.csv: Data from the [St. Louis Fed Economic Data](https://fred.stlouisfed.org/series/USAGDPDEFAISMEI) downloaded [here](https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=USAGDPDEFAISMEI&scale=left&cosd=1960-01-01&coed=2019-01-01&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Annual&fam=avg&fgst=lin&fgsnd=2019-01-01&line_index=1&transformation=lin&vintage_date=2021-01-04&revision_date=2021-01-04&nd=1960-01-01)

* data/EPA_trend: each region-time period combination downloaded from [here](https://www.epa.gov/air-trends/particulate-matter-pm25-trends).

* data/fire/prescribed_burn_acres.csv: Manually annotated from [Wildland Fire Summaries](https://www.nifc.gov/fireInfo/fireInfo_statistics.html).

* data/natural_gas: field level natural gas estimates, aggregated to county level, from the [EIA Natural Gas Annual Respondent Query System](https://www.eia.gov/naturalgas/ngqs/#?report=RP7&year1=2005&year2=2019&company=Name).

* data/physio_shp: USGS shapefile of physio divisions can be downloaded [here](https://water.usgs.gov/GIS/dsdl/physio_shp.zip).

* data/powerplants/emissions_[YYYY].csv: from the [EIA Electricity Data Browser](https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2009&tab=annual_emissions), you need to manually select the year of interest and download from the button above the table. 

* data/powerplants/overview_[YYYY].csv: same process for emissions data above starting from [here](https://www.eia.gov/beta/electricity/data/browser/#/topic/1?agg=2,0,1&fuel=vtvv&sec=g&geo=g&freq=A&datecode=2006&tab=overview&start=200101&end=201710). 

* data/powerplants/Plant_Y[YYYY].xlsx: Plant level data from Form EIA-860 can be downloaded from the [here](https://www.eia.gov/electricity/data/eia860/).

* data/smoke: adapted from data downloaded (individually by day) from [NOAA's Hazard Mapping System](https://www.ospo.noaa.gov/Products/land/hms.html). Processing file is provided in `work/supplemental` in case you would like to update the data in the future.

* data/WUI: state level estimates of number of homes in the wildland urban interface. Generated using National Land Cover Database and proprietary CoreLogic data including the locations of all homes in the US.

### Necessary starting data folder structure to reproduce

```
data
 ├── airport
 │	 ├── airport_locations.csv
 │	 ├── airport_ids.csv
 │	 ├── Q [Quarter]_[Year]_tickets.csv (for 'Q 1' to 'Q 4' 2006-2018)
 │   └── [Month]_[Year]_ontime.csv (for Jan-Dec 2006-2018)
 ├── boundaries
 │   ├── GACC
 │   │	 └── National_GACC_Current_20200226.shp (and associated dbf, prj) 
 │	 ├── tl_2019_us_county
 │	 ├── tl_2010_us_state
 │	 ├── all_national_zips.rds
 │	 ├── counties.RDS
 │	 └── state_fips_codes.csv
 ├── Brey
 │	 ├── [Region]_smokeHourSummary_2007-2014_season=6-9_gdas1_2day+6day_daily_plume=TRUE.RData
 │	 └── BreyStatesRegions.csv
 ├── census
 │	 ├── acres_conventional_tillage.csv 
 │	 ├── ACS_17_1YR_S0101.csv
 │	 ├── BP_[Year]_00A1_with_ann.csv (for 2006-2016)
 │	 ├── CAINC6N__ALL_STATES_2001_2017.csv
 │	 ├── cattle_number.csv
 │	 ├── county_ag_sales.csv
 │	 ├── fertilizer_totals.csv
 │	 └── us_county_average_pm25_income_race.rds
 ├── clean
 │	 └── [Year]_by_county_with_wildfire.csv (for 2008, 2011, 2014, 2017)
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
 ├── vanD
 │	 └── vanD-[Year].asc (2006 - 2018)
 └── WUI
 	 └── [Year]_combined_wui_hh_data.csv (2001, 2004, 2006, 2008, 2011, 2013, 2016)
```
# R Packages needed

R packages required for replications are:

- BAMMtools
- caret
- cleangeo
- data.table
- dplyr
- gdata 
- geosphere
- ggplot2
- ggpubr
- ggthemes
- gridExtra
- hutils
- Hmisc
- imputeTS
- latticeExtra
- mapproj
- maptools
- mapview
- ncdf4
- openxlsx 
- plyr
- raster
- readr
- RColorBrewer
- rgdal
- rgeos
- rnaturalearth
- rnaturalearthdata
- rworldmap
- RSelenium
- sf
- signal 
- sp
- splines
- stringr
- svMisc
- tidyr
- triangle
- truncnorm
- velox
- zoo


Users can run the following one-off command to install the most recent versions of these packages:
```
install.packages(c('BAMMtools','caret','cleangeo','data.table','devtools','dplyr','gdata', 'geosphere','ggplot2','ggpubr','ggthemes','gridExtra',`hutils`,'Hmisc','imputeTS','latticeExtra','mapproj','maptools','mapview','ncdf4','openxlsx','plyr', 'raster','readr','RColorBrewer','rgdal','rgeos','rnaturalearth','rnaturalearthdata','RSelenium','rworldmap','sf','signal', 'sp','stringr','svMisc','tidyr','triangle','truncnorm','velox','zoo'), dependencies = T)
```

Finally, some of the scraping requires the [dev version](https://www.google.com/chrome/dev/) of Google Chrome.


## sessionInfo()

```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] splines   datasets  stats     graphics  grDevices utils     methods   base     

other attached packages:
 [1] RSelenium_1.7.5   BAMMtools_2.1.7   ape_5.3           Hmisc_4.3-0       ggplot2_3.3.2     Formula_1.2-3     survival_3.1-8   
 [8] lattice_0.20-38   imputeTS_3.0      signal_0.7-6      openxlsx_4.1.4    velox_0.2.0       gdata_2.18.0      raster_3.3-13    
[15] readr_1.3.1       stringr_1.4.0     tidyr_1.1.2       dplyr_1.0.2       ncdf4_1.17        geosphere_1.5-10  rgdal_1.5-18     
[22] rgeos_0.5-2       sp_1.4-4          sf_0.8-0          data.table_1.13.2

loaded via a namespace (and not attached):
  [1] backports_1.2.0     plyr_1.8.5          usethis_1.5.1       digest_0.6.27       htmltools_0.4.0     fansi_0.4.1        
  [7] magrittr_1.5        checkmate_1.9.4     memoise_1.1.0       cluster_2.1.0       remotes_2.1.0       xts_0.11-2         
 [13] askpass_1.1         forecast_8.10       tseries_0.10-47     prettyunits_1.1.1   colorspace_1.4-1    rvest_0.3.5        
 [19] pan_1.6             xfun_0.11           callr_3.5.1         crayon_1.3.4        jsonlite_1.7.1      lme4_1.1-25        
 [25] zoo_1.8-6           glue_1.4.2          gtable_0.3.0        pkgbuild_1.1.0      weights_1.0         semver_0.2.0       
 [31] quantmod_0.4-15     jomo_2.6-10         scales_1.1.1        stinepack_1.4       DBI_1.1.0           ggthemes_4.2.0     
 [37] Rcpp_1.0.5          htmlTable_1.13.3    units_0.6-5         foreign_0.8-72      htmlwidgets_1.5.1   httr_1.4.1         
 [43] gplots_3.0.1.1      RColorBrewer_1.1-2  acepack_1.4.1       ellipsis_0.3.1      mice_3.6.0          pkgconfig_2.0.3    
 [49] XML_3.98-1.20       nnet_7.3-12         tidyselect_1.1.0    rlang_0.4.8         munsell_0.5.0       tools_3.6.1        
 [55] cli_2.1.0           generics_0.1.0      audio_0.1-7         devtools_2.2.1      broom_0.7.2         yaml_2.2.0         
 [61] binman_0.1.1        processx_3.4.4      knitr_1.26          fs_1.3.1            zip_2.0.4           caTools_1.17.1.3   
 [67] purrr_0.3.4         mitml_0.3-7         nlme_3.1-142        xml2_1.2.2          compiler_3.6.1      rstudioapi_0.11    
 [73] curl_4.3            e1071_1.7-3         testthat_3.0.0      tibble_3.0.4        statmod_1.4.35      stringi_1.5.3      
 [79] ps_1.4.0            desc_1.2.0          Matrix_1.2-18       classInt_0.4-2      nloptr_1.2.2.2      urca_1.3-0         
 [85] vctrs_0.3.4         pillar_1.4.6        lifecycle_0.2.0     lmtest_0.9-37       bitops_1.0-6        wdman_0.2.4        
 [91] R6_2.5.0            latticeExtra_0.6-28 KernSmooth_2.23-16  gridExtra_2.3       sessioninfo_1.1.1   codetools_0.2-16   
 [97] boot_1.3-23         MASS_7.3-51.4       gtools_3.8.1        assertthat_0.2.1    pkgload_1.1.0       openssl_1.4.1      
[103] rprojroot_1.3-2     withr_2.3.0         fracdiff_1.5-1      parallel_3.6.1      hms_0.5.2           quadprog_1.5-8     
[109] grid_3.6.1          rpart_4.1-15        timeDate_3043.102   class_7.3-15        minqa_1.2.4         TTR_0.23-5         
[115] base64enc_0.1-3 
```