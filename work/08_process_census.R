source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Get census data interpolated at the grid cell level
########################################################################################

###############################################################
# Read in data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")


###############################################################
# Process industry data
###############################################################

out = as.list(rep(NA, sum(years<=2016)))

for (i in 1:length(years)) {
    
    # bring in census data on salary paid within certain categories
    year = years[i]
    file = paste0("data/census/BP_", year, "_00A1_with_ann.csv")
    if (!file.exists(file)) {warning(paste(year, "does not exist.")); next}
    
    # filter to relevant categories and rename cols
    census = read_csv(file)[-1, ]
    census = census[census$NAICS.id %in% c("00", "11", "21", "23", "48-49", "48", 
                                           "49", "53"), ]
    if (year == 2006) {p = "PAYANT"} else {p = "PAYANN"}
    census = census[, c("GEO.id2", "NAICS.display-label", p)]
    names(census) = c("county", "naics", "payroll")
    census$payroll = as.numeric(gsub("[a-zA-Z,.]", "", census$payroll))
    
    # aggregate to the county-naics level
    census = census %>% dplyr::group_by(county, naics) %>% 
        dplyr::summarize(payroll = sum_na(payroll))
    
    # normalize format for combination later
    c = expand.grid(county=unique(census$county), naics=unique(census$naics))
    census = merge(census, c, by=c("county", "naics"), all=T)
    census = spread(census, key="naics", value="payroll")
    names(census) = c("county", "ag", "construction", "mining", "real_estate", "all", "transport")
    census$year = year
    
    out[[i]] = census
}

# combine data and linearly interpolate between missing years
out = rbindlist(out)
out = as.data.frame(out)
out = lm_interpolate(out, "mining", id="county")
out = lm_interpolate(out, "ag", id="county")
out = lm_interpolate(out, "construction", id="county")
out = lm_interpolate(out, "real_estate", id="county")
out = lm_interpolate(out, "all", id="county")

saveRDS(out, "data/model_inputs/census_interpolated.RDS")


###############################################################
# get tourism numbers
###############################################################

tourism = read.csv("data/census/CAINC6N__ALL_STATES_2001_2017.csv")
tourism = tourism[, c("GeoFIPS", "Description", paste0("X", 2006:2017))]
tourism$Description = as.character(tourism$Description)
tourism$type = "" 
tourism[tourism$Description %in% 
            c("    Performing arts, spectator sports, and related industries",
              "    Food services and drinking places"), "type"] = "entertainment"
tourism[tourism$Description %in% 
            c("    Scenic and sightseeing transportation", 
              "    Museums, historical sites, and similar institutions", 
              "    Amusement, gambling, and recreation industries", 
              "    Accommodation"), "type"] = "tourism"
tourism = tourism[tourism$type != "", ]

# convert text to numeric, removing characters first
tourism[, paste0("X", 2006:2017)] = lapply(tourism[, paste0("X", 2006:2017)], 
                                      function(x){as.numeric(gsub("[A-Za-z().,]", "", 
                                                                  as.character(x)))})

# aggregate over the different classes of tourism and entertainment
tourism = tourism %>% group_by(GeoFIPS, type) %>% 
    dplyr::summarize(X2006 = sum_na(X2006), X2007 = sum_na(X2007), X2008 = sum_na(X2008), 
              X2009 = sum_na(X2009), X2010 = sum_na(X2010), X2011 = sum_na(X2011), 
              X2012 = sum_na(X2012), X2013 = sum_na(X2013), X2014 = sum_na(X2014), 
              X2015 = sum_na(X2015), X2016 = sum_na(X2016), X2017 = sum_na(X2017))
tourism$GeoFIPS = str_trim(as.character(tourism$GeoFIPS))

# convert from wide to long data
tourism = reshape2::melt(tourism, c("GeoFIPS", "type"))
tourism$variable = substr(tourism$variable, 2, 5)
tourism = merge(tourism, data_ll[, c("id", "fips")], by.x="GeoFIPS", by.y="fips", all.y=T, all.x=F)

entertainment = tourism[tourism$type == "entertainment" & !is.na(tourism$type), c("id", "variable", "value")]
names(entertainment) = c("id", "year", "entertainment")
entertainment = lm_interpolate(entertainment, "tourism")

tourism = tourism[tourism$type == "tourism" & !is.na(tourism$type), c("id", "variable", "value")]
names(tourism) = c("id", "year", "tourism")
tourism = lm_interpolate(tourism, "tourism")

saveRDS(tourism, "data/model_inputs/tourism_processed.RDS")
saveRDS(entertainment, "data/model_inputs/entertainment_processed.RDS")
