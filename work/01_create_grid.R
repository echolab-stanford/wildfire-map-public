source("work/00_functions.R")
library(Hmisc)

###############################################################
# Read in data
###############################################################

raster = raster("data/pop/gpw_v4_population_count_rev11_2005_2pt5_min.tif")
physio = readOGR("data/physio_shp", "physio", stringsAsFactors=F)

###############################################################
# Create simple county file
###############################################################

if (!file.exists("data/boundaries/counties.RDS")) {
    # The Census TIGER line file is not provided in the repo - if you'd like to 
    # download it to recreate from base data it can be downloaded at:
    # https://www2.census.gov/geo/tiger/TIGER2019/COUNTY/

    counties = readOGR("data/boundaries/tl_2019_us_county", "tl_2019_us_county")
    counties$AWATER = as.numeric(as.character(counties$AWATER))
    counties$ALAND = as.numeric(as.character(counties$ALAND))
    counties = counties[!counties$STATEFP %in% c("02", "15", "60", "66", "69", "72", "78"), ]
    counties = spTransform(counties, CRS(crs_using))
    counties_simp = gSimplify(counties, 0.05, topologyPreserve=T)
    counties_simp = SpatialPolygonsDataFrame(counties_simp, counties@data)
    saveRDS(counties_simp, "data/boundaries/counties.RDS")
} else {
    counties_simp = readRDS("data/boundaries/counties.RDS")
}


###############################################################
# Use raster to make grid, add state, county, fips 
###############################################################

# create the grid cells 
raster = crop(raster, extent(counties_simp))
km = round((km/111/2)/res(raster)[1])
raster = aggregate(raster, fact=km, fun=sum)
saveRDS(raster, "data/clean/national_grid_raster.RDS")

# remove grid cells that don't overlap a state
grid = rasterToPolygons(raster)
keep = over(grid, counties_simp)
keep = which(!is.na(keep$STATEFP))
data_ll = grid[keep, ]

# create IDs for the grid cells 
data_ll$id = str_pad(1:nrow(data_ll), 4, "0", side="left")
row.names(data_ll) = data_ll$id
data_ll@data[, c("lon", "lat")] = coordinates(data_ll)
data_ll = data_ll[, -1]

# allocate grid cells to counties
c = over(data_ll, counties_simp)
data_ll$fips = c$GEOID
data_ll$county = c$NAME
data_ll$state = c$STATEFP
data_ll = data_ll[!is.na(data_ll$county), ]
data_ll$county = as.character(data_ll$county)

# fix a few county names for merging later
data_ll$county[data_ll$county == "Richmond"] = "City of Richmond"
data_ll$county[data_ll$county == "Denali"] = "Denali Borough"
data_ll$county[data_ll$county == "LaPorte"] = "La Porte"
data_ll$county[data_ll$county == "St. Clair"] = "St Clair"
data_ll$county[data_ll$county == "St. Charles"] = "St Charles"
data_ll$county[data_ll$county == "St. Joseph"] = "St Joseph"
data_ll$county[data_ll$county == "St. Louis"] = "St Louis"

# get state names and finalize column names
fips = read_csv("data/boundaries/state_fips_codes.csv")[, 1:2]
data_ll = merge(data_ll, fips, by.x="state", by.y="fips", all.x=T)
data_ll = data_ll[, c("id", "lon", "lat", "fips", "county", "state.y")]
names(data_ll) = c("id", "lon", "lat", "fips", "county", "state")


###############################################################
# Create physio zones across US and add to raster
###############################################################

# Physiographic regions are in: Division, Province, Section (where sections is largest)
# Creating physiographic divisions that are similar sized by subdividing in certain areas.

physio = physio[, c("AREA", "PROVINCE", "SECTION", "DIVISION")]
physio@data[is.na(physio$SECTION), "SECTION"] = physio@data[is.na(physio$SECTION), "PROVINCE"]
physio = physio[!is.na(physio$DIVISION), ]
physio = spTransform(physio, CRS(crs_using))

# count grid cells in each province (range 2 - 827), in areas where less than 200 cells in
# each province, don't use the lower level of sections, no need to further subdivide
counts = table(over(data_ll, physio)$PROVINCE)
physio@data[physio$PROVINCE %in% names(counts[counts <= 200]), "SECTION"] = 
    physio@data[physio$PROVINCE %in% names(counts[counts <= 200]), "PROVINCE"]

# for small sections, merge with their smaller neighboring section
physio[physio$SECTION == "OLYMPIC MOUNTAINS", "SECTION"] = "OREGON COAST RANGE"
physio[physio$SECTION == "OUACHITA", "SECTION"] = "OZARK PLATEAUS"
physio[physio$SECTION == "ST. LAWRENCE VALLEY", "SECTION"] = "ADIRONDACK"
physio[physio$SECTION == "RATON", "SECTION"] = "COLORADO PIEDMONT"
physio[physio$SECTION == "BLACK HILLS", "SECTION"] = "MISSOURI PLATEAU, UNGLACIATED"
physio[physio$SECTION == "SACRAMENTO", "SECTION"] = "PECOS VALLEY"
physio[physio$SECTION == "BLUE RIDGE", "SECTION"] = "VALLEY AND RIDGE"
physio[physio$SECTION %in% c("EDWARDS PLATEAU", "CENTRAL TEXAS"), "SECTION"] = "WEST GULF COASTAL PLAIN"
physio[physio$SECTION %in% c("SALTON TROUGH", "LOS ANGELES RANGES"), "SECTION"] = "LOWER CALIFORNIAN"
physio[physio$DIVISION == "INTERIOR HIGHLANDS"] = "ATLANTIC PLAIN"
physio[physio$DIVISION == "LAURENTIAN UPLAND"] = "INTERIOR PLAINS"

# assign each grid cell their physiographic location and save
groups = over(data_ll, physio)
data_ll$physio_section = groups$SECTION
data_ll$physio_region = groups$DIVISION
data_ll = data_ll[!is.na(data_ll$physio_section), ]
data_ll$physio_section = as.factor(data_ll$physio_section)
data_ll$physio_region = as.factor(data_ll$physio_region)
saveRDS(data_ll, "data/clean/national_grid.RDS")


###############################################################
# Get the EPA data in to the grid form
###############################################################

# This data is derived from the EPA data downloaded at: 
# https://www.epa.gov/outdoor-air-quality-data/download-daily-data
# tools to help download that raw data can be found on GitHub in BurkeLab/census.tools

# read in EPA PM2.5 data
epa = readRDS("data/pm/epa_station_level_pm25_data.rds")
epa = dplyr::select(epa, year, month, day, lat, lon, pm25, aqi)
epa = epa %>% dplyr::filter(year %in% years)
epa$lat = round(epa$lat, 3)
epa$lon = round(epa$lon, 3)
epa$id = epa %>% dplyr::group_indices(lat, lon)
epa_ll =  epa[!duplicated(epa$id), c("lon", "lat", "id")]
epa_ll = SpatialPointsDataFrame(epa_ll[, c("lon", "lat")], data=epa_ll)
proj4string(epa_ll) = CRS(crs_using)

# get the pm from stations for each grid and take mean of all stations within grid to get 
# yearly PM number for each cell
matches = over(data_ll, epa_ll, returnList=T) 
data = expand.grid(id=data_ll$id, year=years)
data$pm = NA
data$obs = NA

#runs slowly (about 3 minutes)
prog = txtProgressBar(min=0, max=length(matches), initial=0, char="-", style=3)
for (i in 1:length(matches)) {
    cur = matches[[i]]
    if (nrow(cur) > 0) {
        
        # filter by epa data from correct stations
        cur = epa %>% dplyr::filter(id %in% cur$id) %>% 
            # take daily mean by station
            group_by(year, month, day, id) %>%
            dplyr::summarise(pm = mean_na(pm25), n=n()) %>% 
            # take daily mean 
            group_by(year, month, day) %>%
            dplyr::summarise(pm = wtd.mean(pm, weights=n), n=sum(n)) %>% 
            # take monthly mean
            group_by(year, month) %>%
            dplyr::summarise(pm = wtd.mean(pm, weights=n), n=n()) %>%
            # overall year mean for that grid cell
            group_by(year) %>% 
            dplyr::summarise(pm = wtd.mean(pm, weights=n), obs=sum(n))
        id = names(matches)[[i]]
        cur$id = id
        cur = cur[, c("id", "year", "pm", "obs")]
        
        data[data$id == id & data$year %in% cur$year, ] = cur
        
        if (i %% 30 == 0) {setTxtProgressBar(prog, i)}
    }
}
saveRDS(data, "data/clean/gridded_pm_data.RDS")


###############################################################
# Create adjacency matrix
###############################################################

touching = gTouches(data_ll, byid=T)
adjacency = apply(touching, 1, which)
saveRDS(adjacency, "data/clean/grid_adjacency.RDS")

adjacency2 = adjacency
for (i in 1:length(adjacency)) {
    cur = adjacency[[i]]
    added = unique(unlist(adjacency[cur]))
    
    cur = c(cur, added)
    cur = cur[cur != i & !(cur %in% adjacency[[i]])]
    names(cur) = str_pad(cur, 4, "left", "0")
    adjacency2[[i]] = cur
}
saveRDS(adjacency2, "data/clean/grid_adjacency2.RDS")

adjacency3 = adjacency
for (i in 1:length(adjacency)) {
    cur = adjacency2[[i]]
    added = unique(unlist(adjacency[cur]))
    
    cur = c(cur, added)
    cur = cur[cur != i & !(cur %in% adjacency[[i]]) & !(cur %in% adjacency2[[i]])]
    names(cur) = str_pad(cur, 4, "left", "0")
    adjacency3[[i]] = cur
}
saveRDS(adjacency3, "data/clean/grid_adjacency3.RDS")
