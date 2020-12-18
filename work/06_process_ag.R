source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Get ag data interpolated at the grid cell level
########################################################################################

###############################################################
# Read in data
###############################################################

data_ll = readRDS("data/clean/national_grid.RDS")
ag = read_csv("data/census/county_ag_sales.csv")
tillage = read_csv("data/census/acres_conventional_tillage.csv")
cattle = read_csv("data/census/cattle_number.csv")
fertilizer = read_csv("data/census/fertilizer_totals.csv")


###############################################################
# Process AG
###############################################################

# get the state + county identifiers all lined up
d_ag = ag[!is.na(ag$`County ANSI`), ]
d_ag$fips = paste0(d_ag$`State ANSI`, d_ag$`County ANSI`)
d_ag = d_ag[, c("Year", "fips", "Value")]
names(d_ag) = c("year", "fips", "value")
d_ag$value = as.numeric(gsub("[a-zA-z.,()]", "", d_ag$value))
d_ag = lm_interpolate(d_ag, "value", "fips")
names(d_ag) = c("id", "year", "ag")


###############################################################
# Process cattle
###############################################################

# get the state + county identifiers all lined up
cattle$fips = paste0(cattle$`State ANSI`, cattle$`County ANSI`)
cattle = cattle %>% dplyr::filter(!is.na(`County ANSI`) & Program == "CENSUS") %>% 
    dplyr::select(Year, fips, Value)
names(cattle) = c("year", "fips", "cattle")
cattle$cattle = as.numeric(gsub("[a-zA-Z,()?]", "", as.character(cattle$cattle)))
cattle = lm_interpolate(cattle, "cattle", "fips")


###############################################################
# Process fertilizer
###############################################################

# get the state + county identifiers all lined up
fertilizer$fips = paste0(fertilizer$`State ANSI`, fertilizer$`County ANSI`)
fertilizer = fertilizer %>% 
    dplyr::filter(!is.na(`County ANSI`) & 
           `Data Item` == "FERTILIZER TOTALS, INCL LIME & SOIL CONDITIONERS - EXPENSE, MEASURED IN $") %>% 
    dplyr::select(Year, fips, Value)
names(fertilizer) = c("year", "fips", "fertilizer")
fertilizer$fertilizer = as.numeric(gsub("[a-zA-Z,()?]", "", as.character(fertilizer$fertilizer)))
fertilizer = lm_interpolate(fertilizer, "fertilizer", "fips")


# combine all the different ag files
d_ag = merge(d_ag, fertilizer, by.x=c("year", "id"), by.y=c("year", "fips"), all=T)
d_ag = merge(d_ag, cattle, by.x=c("year", "id"), by.y=c("year", "fips"), all=T)

fips = expand.grid(fips=unique(data_ll$fips), year=years)
fips = merge(data_ll@data[, c("id", "fips")], fips, by="fips", all.x=T)
d_ag = merge(fips, d_ag, by.x=c("fips", "year"), by.y=c("id", "year"))
pad = expand.grid(id=unique(d_ag$id), year=years)
pad = merge(pad, unique(d_ag[, c("id", "fips")]), by="id", all=T)
d_ag = merge(d_ag, pad, by=c("id", "year", "fips"), all=T) %>% dplyr::select(-fips)

saveRDS(d_ag, "data/model_inputs/ag_processed.RDS")


###############################################################
# Process tillage
###############################################################

tillage$fips = paste0(tillage$`State ANSI`, tillage$`County ANSI`)
tillage = tillage %>% dplyr::filter(!is.na(`County ANSI`)) %>% dplyr::select(fips, Value)
names(tillage) = c("fips", "tillage")
tillage$tillage = as.numeric(gsub("[a-zA-Z,()?]", "", as.character(tillage$tillage)))

fips = unique(fips[, c("fips", "id")])
tillage = merge(fips, tillage, by="fips")
tillage = dplyr::select(tillage, id, tillage)

saveRDS(tillage, "data/model_inputs/tillage_processed.RDS")
