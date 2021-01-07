source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Get coal, natural gas and other plant data interpolated at the grid cell level
########################################################################################

###############################################################
# Read in data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")
zips = readRDS("data/boundaries/all_national_zips.rds")
cents = coordinates(zips)
zips = SpatialPointsDataFrame(coords=cents, data=zips@data, 
                              proj4string=CRS(crs_using))


###############################################################
# Process coal and power plant locations
###############################################################

# loop through all the Excel files with information on coal and power plants
# to create a shapefile of their locations (georeferenced to centroid of zip code)
plants_out = list(length(years))
coal_out = list(length(years))

for (i in 1:length(years)) {
    year = years[i]
    if (year >= 2011) {
        if (year >= 2012) {k = c("Plant.Code", "State", "Zip")} else {
            k = c("PLANT_CODE", "STATE", "ZIP5")}
        file = paste0("data/powerplants/Plant_Y", year, ".xlsx")
        plants = read.xlsx(file, startRow=2)
        plants = plants[, k]
    } else {
        if (year == 2006) {k = c("PLNTCODE", "STATE", "PLNTZIP")} else if (year <= 2008) {
            k = c("PLNTCODE", "STATE", "ZIP5")} else {
                k = c("PLANT_CODE", "STATE", "ZIP5")}
        file = paste0("data/powerplants/Plant_Y", year, ".xls")
        plants = read.xls(file)
        plants = plants[, k]
    }
    plants_out[[i]] = plants
    
    #### do the same for hte location of coal plants specifically
    file = paste0("data/coal/existing_gen_units_", year, ".xls")
    if (year >= 2015) {file = paste0(file, "x")}
    if (!file.exists(file)) {warning(paste(year, "does not exist.")); next}
    
    if (year >= 2012) {skip=1} else {skip=3}
    if (year <= 2014) {
        temp = read.xls(file, skip=skip, header=T)
    } else {
        temp = read.xlsx(file, startRow=skip+1)
    }
    coal_out[[i]] = temp
}

# make a flat file of spatial points
plants = rbindlist(plants_out, use.names=F)
plants = unique(plants)
names(plants) = c("id", "state", "zip")
plants$zip = str_pad(plants$zip, 5, "left", "0")
plants = merge(plants, zips, by.y="ZCTA5CE10", by.x="zip")
plants$INTPTLAT10 = as.numeric(as.character(plants$INTPTLAT10))
plants$INTPTLON10 = as.numeric(as.character(plants$INTPTLON10))
plants = SpatialPointsDataFrame(plants[, c("INTPTLON10", "INTPTLAT10")], 
                                plants[, c("id", "state", "zip")], 
                                proj4string=CRS(crs_using))

# flatten data for coal plants 
coal_plants = rbindlist(coal_out, fill=T)
coal_plants = coal_plants[, c("Plant.ID", "Plant.Name")]
names(coal_plants) = c("id", "name")
coal_plants = unique(coal_plants)

###############################################################
# Process power plant emissions
###############################################################

#loop through years to count emissions
output = list(length(years))
prog = txtProgressBar(min=0, max=length(years), initial=0, char="-", style=3)
for (i in 1:length(years)) {
    
    year = years[i]
    
    temp_df = data.frame(id=as.character(data_ll@data[["id"]]), year=as.character(year))
    
    # read in overview to get hte zip code
    over = read.csv(paste0("data/powerplants/overview_", year, ".csv"), skip=4)
    
    # read in emissions to get the CO2 emissions fo the plant
    emissions = read.csv(paste0("data/powerplants/emissions_", year, ".csv"), skip=4)
    emissions = merge(over, emissions, by=c("Plant.Name", "State"))
    emissions = dplyr::select(emissions, State, Plant.Code, Sector.Name, CO2.Emissions..tons.)
    names(emissions) = c("state", "id", "sector", "co2")
    emissions = emissions %>% group_by(state, id) %>% 
        dplyr::summarize(co2=sum(as.numeric(co2), na.rm=T))
    
    emissions = merge(plants, emissions, by=c("id", "state"))
    emissions = emissions[!is.na(emissions$co2), ]
    coalemissions = emissions[emissions$id %in% coal_plants$id, ]
    emissions = emissions[!emissions$id %in% coal_plants$id, ]
    
    # allocate emissions to grid cells
    count = over(data_ll, emissions, fun=sum, na.rm=T)
    count$co2[is.na(count$co2)] = 0
    temp_df[, "plant"] = count$co2
    
    # same for coal
    count = over(data_ll, coalemissions, fun=sum, na.rm=T)
    count$co2[is.na(count$co2)] = 0
    temp_df[, "coalplant"] = count$co2
    
    output[[i]] = temp_df
    setTxtProgressBar(prog, i)
}

# flatten data
data = rbindlist(output)
saveRDS(data, "data/model_inputs/powerplants_processed.RDS")


###############################################################
# Process coal plants
###############################################################

output = list(length(years))

prog = txtProgressBar(min=0, max=length(years), initial=0, char="-", style=3)
for (i in 1:length(years)) {
    
    year = years[i]
    file = paste0("data/coal/existing_gen_units_", year, ".xls")
    if (year >= 2015) {file = paste0(file, "x")}
    if (!file.exists(file)) {warning(paste(year, "does not exist.")); next}
    
    if (year >= 2012) {skip=1} else {skip=3}
    if (year <= 2014) {
        temp = read.xls(file, skip=skip, header=T)
    } else {
        temp = read.xlsx(file, startRow=skip+1)
    }
    
    temp = temp %>% dplyr::group_by(County, State) %>% dplyr::summarize(n = n())
    names(temp) = c("county", "state", "coalplant")
    temp = merge(temp, data_ll@data, by=c("county", "state"))
    temp$lat = as.numeric(as.character(temp$lat))
    temp$lon = as.numeric(as.character(temp$lon))
    temp = SpatialPointsDataFrame(temp[, c("lon", "lat")], 
                                  data.frame(count=temp$coalplant), 
                                  proj4string=CRS(crs_using))
    count = over(data_ll, temp, returnList=T)
    
    temp_df = data.frame(id=data_ll@data[, c("id")])
    temp_df$year = year
    temp_df$coalcount = sapply(count, FUN=function(x){sum(x$count)})
    
    knn_results = RANN::nn2(data=coordinates(temp), query=coordinates(data_ll), k=2)
    temp_df$coaldist1 = as.vector(knn_results$nn.dists[,1]) * 111
    temp_df$coaldist2 = (as.vector(knn_results$nn.dists[,2]) * 111) - temp_df$coaldist1
    
    output[[i]] = temp_df
    setTxtProgressBar(prog, i)
}

# flatten data
data = rbindlist(output)
saveRDS(data, "data/model_inputs/coal_processed.RDS")
