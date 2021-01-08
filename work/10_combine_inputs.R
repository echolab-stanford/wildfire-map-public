source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Combine the various co-variates that are needed for modelling and combine
########################################################################################

###############################################################
# read in geo data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")
data = readRDS("data/clean/gridded_pm_data.RDS")


###############################################################
# get input_data
###############################################################

# fire and smoke variables
d_smoke = readRDS("data/model_inputs/smoke_processed.RDS")
d_smoke$year = as.numeric(as.character(d_smoke$year))

# other yearly covariates
d_pop = readRDS("data/model_inputs/pop_processed.RDS")
d_plants = readRDS("data/model_inputs/powerplants_processed.RDS")
d_coal = readRDS("data/model_inputs/coal_processed.RDS")
d_ag = readRDS("data/model_inputs/ag_processed.RDS")
d_airports = readRDS("data/model_inputs/airport_interpolated.RDS")
d_traffic = readRDS("data/model_inputs/traffic_processed.RDS")
d_tourism = readRDS("data/model_inputs/tourism_processed.RDS")
d_entertainment = readRDS("data/model_inputs/entertainment_processed.RDS")


d_census = readRDS("data/model_inputs/census_interpolated.RDS")
d_census = merge(d_census, unique(data_ll@data[, c("id", "fips")]), 
                 by.x="county", by.y="fips", all.y=T)
d_census = dplyr::select(d_census, -county, -ag, -all, -transport)

# process natural gas data
d_gas = readRDS("data/natural_gas/county_level_natural_gas_estimates_2006_2016.rds")
d_gas[d_gas$state == "IL" & d_gas$county_name == "Mclean", "county_name"] = "McLean"
d_gas[d_gas$state == "IL" & d_gas$county_name == "Warren/Mercer", "county_name"] = "Mercer"
d_gas[d_gas$state == "MT" & d_gas$county_name == "Glacier National", "county_name"] = "Glacier"
d_gas[which(d_gas$state == "PA" & d_gas$county_name == "Mckean"), "county_name"] = "McKean"
d_gas = merge(d_gas, unique(data_ll@data[, c("id", "county", "state")]), 
              by.x=c("state", "county_name"), by.y=c("state", "county"))
pad = expand.grid(id=data_ll$id, year=years)
d_gas = merge(d_gas, pad, by=c("id", "year"), all=T) %>% 
    dplyr::select(-state, -county_name)
d_gas = lm_interpolate(d_gas, "tot_production")

#the ids with no data now is bcz there is no production wells there
d_gas[is.na(d_gas)] = 0
d_gas = d_gas %>% dplyr::select(id, year, tot_production)
names(d_gas)[3] = "gas"


###############################################################
# combine data and convert data to long format
###############################################################

data$id = as.character(data$id)
d_coal$id = as.character(d_coal$id)
d_coal$year = as.character(d_coal$year)
d_smoke$id = as.character(d_smoke$id)
d_pop$id = as.character(d_pop$id)
d_plants$id = as.character(d_plants$id)
d_traffic$id = as.character(d_traffic$id)

# combine data
full_data = Reduce(function(x, y) merge(x, y, by=c("id", "year"), all=TRUE), 
              list(data, d_smoke, d_pop, d_plants, d_coal, 
                    d_airports, d_gas, d_census, d_traffic, d_ag, d_tourism, 
                   d_entertainment)) 
full_data = full_data[full_data$year >= as.numeric(min(years)) &
                          full_data$year <= as.numeric(max(years)), ]

# merge in physio_regions
full_data = Reduce(function(x, y) merge(x, y, by=c("id"), all=TRUE), 
                   list(full_data, 
                        unique(data_ll@data[, c("id", "state", "physio_section", "physio_region")])))
full_data = full_data[!is.na(full_data$physio_section), ] #removes great lakes
saveRDS(full_data, "data/clean/epa_full_long_data.RDS")

###############################################################
# impute
###############################################################

# get list of variables that shouldn't be used for interpolation
vars = names(full_data)
vars = vars[grepl("fire", vars) | grepl("smoke", vars) | grepl("den", vars)]
vars = c("id", "pm", "physio_section", "state", "obs", "area", vars)

# find which columns have missing data, and sort from those missing the least to the most
missing = apply(full_data[, !names(full_data) %in% vars], 2, function(x){sum(is.na(x))})
missing = sort(missing)
missing = names(missing[missing > 0])

# loop through the missing variables and interpolate
for (miss in missing) {
    if (miss == "pm") {next}
    v = vars[vars != miss]
    m = interpolate_cases(full_data, miss, v)
    w = is.na(full_data[, miss])
    full_data[w, miss] = predict(m, full_data[w, ])
    full_data[w & full_data[, miss] < 0, miss] = 0
}

saveRDS(full_data, "data/clean/epa_full_long_data_imputed.RDS")

###############################################################
# get spatial lag
###############################################################

adjacency_vars = c("pop", "coalplant", "mining")
adjacency_vars = adjacency_vars[adjacency_vars %in% names(full_data)]

adj1 = readRDS("data/clean/grid_adjacency.RDS")
adj2 = readRDS("data/clean/grid_adjacency2.RDS")
adj3 = readRDS("data/clean/grid_adjacency3.RDS")


# run that function over each list of adjacent cells
temp = lapply(1:length(adj1), FUN=get_adj)
temp = rbindlist(temp)

d = merge(full_data, temp, by=c("id", "year"))
saveRDS(d, "data/clean/epa_full_long_data_spatial_lag.RDS")
