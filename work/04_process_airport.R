source("work/00_functions.R")

########################################################################################
# Written by: Anne Driscoll
# Get airport data interpolated at the grid cell level
########################################################################################

files = list.files("data/airport")

###############################################################
# Get number of flights by airport
###############################################################

files_tix = files[grepl("tickets", files)]
num_flight = as.list(rep(NA, length(files_tix)))
for (i in 1:length(files_tix)) {
    f = files_tix[i]
    cur = readRDS(paste0("data/airport/", f))
    cur = cur %>% group_by(YEAR, ORIGIN, ORIGIN_AIRPORT_ID) %>%
        dplyr::summarize(n = sum(PASSENGERS, na.rm=T))
    
    num_flight[[i]] = cur
}

num_flight = rbindlist(num_flight)
num_flight = num_flight %>% group_by(YEAR, ORIGIN) %>%
    dplyr::summarize(tickets = sum(n, na.rm=T)) #becuase each file is a quarter

###############################################################
# Get minutes of taxiing by airport
###############################################################

files = files[grepl("ontime", files)]
min_taxi = as.list(rep(NA, length(files)))
done = expand.grid(year=years, month=1:12, done=0)
for (i in 1:length(files)) {
    f = files[i]
    cur = readRDS(paste0("data/airport/", f))
    exists = done[done$month == cur$MONTH[1] & done$year == cur$YEAR[1], ]$done
    
    if(exists == 0) {
        cur = cur %>% 
            group_by(YEAR, MONTH, ORIGIN, DEST) %>%
            dplyr::summarize(taxi_in = sum(TAXI_IN, na.rm=T), 
                      taxi_out = sum(TAXI_OUT, na.rm=T))
        done[done$month == cur$MONTH[1] & done$year == cur$YEAR[1], ]$done = 1
    } else { next }

    min_taxi[[i]] = cur
}

w = sapply(min_taxi, function(x) {any(!is.na(x))})
min_taxi = rbindlist(min_taxi[w])

#combine the taxi in + taxi out, then aggregate
min_taxi = rbind(min_taxi[, c("YEAR", "MONTH", "ORIGIN", "taxi_out")], 
          min_taxi[, c("YEAR", "MONTH", "DEST", "taxi_in")], use.names=F)
min_taxi = min_taxi %>% group_by(YEAR, ORIGIN) %>%
    summarize(taxi = sum(taxi_out, na.rm=T)) #to combine taxi in and taxi out


###############################################################
# Go from airport ll to grid id
###############################################################

locations = read.csv("data/airport/airport_locations.csv", stringsAsFactors=F)
locations  = locations[locations$country == "United States", c("code", "lat", "lon")]
locations = rbind(locations, c("BIC", 46.7748198,- 100.7574257), 
                  c("BKX", 44.303667,-96.8136067), c("BQN", 18.4953957,-67.1377666), 
                  c("EAR", 40.7251123,-99.0111922), c("GUM", 13.4852976,144.7986232), 
                  c("IFP", 35.1562931,-114.5606222), c("IGM", 35.2589259,-113.9456463), 
                  c("ILE", 31.0873439,-97.6876582), c("IWD", 46.5250178,-90.133080), 
                  c("KKI", 60.8628897,-161.4807899), c("MAZ", 18.255084,-67.1506746), 
                  c("PLB", 44.6540516,-73.4638638), c("WDG", 36.3785883,-97.7928685), 
                  c("PPG", -14.331389,-170.71357), c("PSE", 18.0090714,-66.56476), 
                  c("ROP", 14.172008,145.2415723), c("SJU", 18.4383715,-66.002289), 
                  c("SPN", 15.1197429,145.7260902), c("STT", 18.3360608,-64.974461), 
                  c("STX", 17.6995292,-64.79965), c("SVC", 32.6480879,-108.15538), 
                  c("TSS", 40.7425701,-73.97413), c("VQS", 18.1352726,-65.493517), 
                  c("JQF", 35.3804727,-80.7161548), c("PJO", 68.348878,-166.800906), 
                  c("SXP", 62.52027,-164.85002))
#removed "WAS" and "NYC" bcz they double count the airports there

ll = data.frame(ORIGIN = unique(c(as.character(num_flight$ORIGIN), 
                                  as.character(min_taxi$ORIGIN))))
ll = merge(ll, locations, by.x="ORIGIN", by.y="code", all.x=T)
ll$lat = as.numeric(ll$lat)
ll$lon = as.numeric(ll$lon)

airport_ll = unique(ll[, c("ORIGIN", "lon", "lat")])
airport_ll = airport_ll[!is.na(airport_ll$lat), ] #remove some virgin islands, PR
airport_ll = SpatialPointsDataFrame(SpatialPoints(airport_ll[,c("lon", "lat")]), 
                                    data.frame(origin=airport_ll$ORIGIN))
proj4string(airport_ll) = CRS(crs_using)

data_ll = readRDS("data/clean/national_grid.RDS")
a = over(data_ll, airport_ll, returnList=T)
names(a) = data_ll$id

for (i in 1:length(a)) {
    cur = a[[i]]
    flights = num_flight %>% filter(ORIGIN %in% cur$origin) %>%
        group_by(YEAR) %>% summarize(tickets = sum(tickets))
    taxi = min_taxi %>% filter(ORIGIN %in% cur$origin) %>%
        group_by(YEAR) %>% summarize(taxi = sum(taxi))
    cur = merge(flights, taxi, by=c("YEAR"), all=T)
    if (nrow(cur) == 0) {cur = data.frame(year=years, tickets=0, taxi=0)}
    if (nrow(cur) == 1) {cur = data.frame(year=years, tickets=cur$tickets, taxi=cur$taxi)}
    if (nrow(cur) > 1) {
        names(cur) = c("year", "tickets", "taxi")
        cur$id = 1
        cur = kalman_mean_interpolate(cur, "tickets")
        cur = kalman_mean_interpolate(cur, "taxi")
        cur = cur %>% dplyr::select(-id)
    }
    cur$id = data_ll$id[i]
    a[[i]] = cur
}

d_airports = rbindlist(a)
names(d_airports) = c("year", "tickets", "taxi", "id")
saveRDS(d_airports, "data/model_inputs/airport_interpolated.RDS")
