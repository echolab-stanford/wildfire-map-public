source("work/00_functions.R")

###############################################################
# Read in data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")


###############################################################
# Process population
###############################################################

pop_folder="data/pop"
years=c(2005, 2010, 2015, 2020)
    
output = list(length(years))
geo_df = data.frame(id = unique(data_ll$id), year=as.numeric(NA), pop=as.numeric(NA))

#loop through the pop files  (< 1 minute)
prog = txtProgressBar(min=0, max=length(years), initial=0, char="-", style=3)
for (i in 1:length(years)){
    
    year = years[i]
    file = paste0("data/pop/gpw_v4_population_count_rev11_", as.character(year), 
                  "_2pt5_min.tif")
    if (!file.exists(file)) {warning(paste(year, "does not exist.")); next}
    
    # extract the summed population numbers for each grid cell
    pop = raster(file)
    pop = crop(pop, extent(data_ll))
    pop = velox(pop)
    pop_ex = pop$extract(data_ll, fun=function(x) sum(x,na.rm=T))
    
    # create dataframe to save that info under
    temp_df = geo_df
    temp_df$year = as.character(year)
    temp_df$pop = as.vector(pop_ex)
    temp_df$id = as.character(temp_df$id)
    
    output[[i]] = temp_df
    
    setTxtProgressBar(prog, i)
}

# combine the data for every grid-cell year in one data file
data = rbindlist(output)
data$id = as.character(data$id)
data$year = as.character(data$year)
pad = expand.grid(id=unique(data$id), year=as.character(2006:2016))
data = merge(data, pad, by=c("id", "year"),  all=T)
data$year = as.numeric(as.character(data$year))

d_pop = lm_interpolate(data, "pop")

saveRDS(d_pop, "data/model_inputs/pop_processed.RDS")
