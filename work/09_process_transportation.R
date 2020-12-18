source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Get transportation data interpolated at the grid cell level
########################################################################################

###############################################################
# Read in data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")


###############################################################
# Get traffic for every year
###############################################################

# this is really slow (one hourish)
output = as.list(rep(NA, length(years)))
prog = txtProgressBar(min=0, max=length(years), initial=0, char="-", style=3)
for (i in 1:length(years)) {
    
    year = years[i]
    
    # try to load in the raster and convert to velox
    traffic = tryCatch({raster(paste0("data/traffic-darte/onroad_", year, ".tif"))}, 
                       error = function(e) {next})
    traffic = projectRaster(traffic, crs=CRS(crs_using))
    traffic = velox(traffic)
    t = traffic$extract(data_ll, fun=sum_na, small=T)
    
    # make the extraction to a data frame of the right format
    df = data.frame(id=data_ll$id, year=year, traffic=t)
    output[[i]] = df
    setTxtProgressBar(prog, i)
}

# flatten and save the traffic data
data = rbindlist(output)
saveRDS(data, "data/model_inputs/traffic_processed.RDS")
