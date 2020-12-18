source("work/00_functions.R")
library(Hmisc)
library(BAMMtools)

########################################################################################
# Written by: Anne Driscoll
# Gets clusters from the fire points to try to identify large fires (eg Camp Fire)
# Combines smoke plume data with fire to get a crude estimate of which fire the smoke 
#   plume originated from
########################################################################################

###############################################################
# Functions
###############################################################

get_km_bounds = function(v, k, iter=15, n=3) {
    v = v[!is.na(v)]
    km = kmeans(v, 3, nstart=n, iter.max=iter)
    bounds = rep(NA, k)
    for (i in 1:k) {
        bounds[i] = max(v[km$cluster==i])
    }
    bounds = sort(bounds)[1:(k-1)]
    return(bounds)
}

#THIS IS A REALLY SLOW FILE, STUFF HASNT BEEN OPTIMIZED WELL
#FOR NOW YOU JUST GOTTA SIT THROUGH IT. **RUN OVERNIGHT**

###############################################################
# Read in data
###############################################################

# Original smoke data at: https://www.ospo.noaa.gov/Products/land/hms.html
# here all the individual day files are stored in a list
smoke = read_rds("data/smoke/smoke_plumes.rds")
smoke = smoke[names(smoke) < paste0(max(as.numeric(years))+1,"0101")]

# Original fire data at: https://www.ospo.noaa.gov/Products/land/hms.html
# here all the individual day files are stored in a list
fire = read_rds("data/fire/hms_fires.RDS")
fire = fire[names(fire) < paste0(max(as.numeric(years))+1,"0101")]
data_ll = readRDS("data/clean/national_grid.RDS")

# width is the number fire points are buffered to merge into fire clusters
width = 2.9/111 #km/number of km at eq (to get in lat lon units)
#2.9 because it's a 4km grid cell, to reach diagonal need sqrt(2^2+2^2)


###############################################################
# Process fire
###############################################################

# get data for each grid cell, for each year. 
j = 1

# loop through years
for (i in 1:length(years)) {
    
    # get the fires that happened during the year of interest
    y = years[i]
    year_fire = grepl(paste0("^", y), names(fire))
    first = max(c(1, Position(function(x){x==T}, year_fire)-3))
    last = max(c(length(year_fire) - Position(function(x){x==T}, rev(year_fire)) + 1))
    year_fire = fire[first:last]
    year_fire_sp = as.list(rep(NA, length(year_fire)))
    start_loop = ifelse(first == 1, 1, 4)
    
    # loop through days in the year
    prog = txtProgressBar(min=0, max=length(year_fire), initial=0, char="-", style=3)
    for (k in start_loop:length(year_fire)) {
        
        #sf to sp
        date = names(year_fire)[[k]]
        
        # using fire data for the 3 days previous as well to capture long burning fires. 
        # figure out which days to combine to get the correct set.
        if (k>4) {
            f = year_fire[(k-3):k]
            f = rbind(f[[1]][, c("geometry")], 
                      f[[2]][, c("geometry")], 
                      f[[3]][, c("geometry")], 
                      f[[4]][, c("geometry")])
        } else {f = year_fire[[k]][,  c("geometry")]}
        if (nrow(f)==0) {next}
        
        # convert f to an SP object
        f = tryCatch({
            as_Spatial(f)
        }, error = function(e) {
            f = f[!is.na(st_is_valid(f)) & !st_is_empty(f),]
            f = st_cast(f, "POINT")
            as_Spatial(f)
        })
        
        #buffer by 'width' to merge adjacent pixels
        f_buf = gBuffer(f, byid=T,  width=width, capStyle="SQUARE", quadsegs=1)
        f_buf = st_cast(st_union(st_as_sf(f_buf)), "POLYGON")
        f_buf = SpatialPolygonsDataFrame(as_Spatial(f_buf), data.frame(id=1:length(f_buf)), match.ID=F)
        f_buf$area = gArea(f_buf, byid=T)*12321
        f_buf$date = date
        f = f_buf
        
        #create new ids for clustered fires so that they can all be merged at the end
        f$ID = j:(nrow(f)+j-1)
        row.names(f) = as.character(j:(length(f$ID)+j-1))
        j = length(f$ID)+j
        year_fire_sp[[k]] = f
        
        if (k %% 5 == 0) {setTxtProgressBar(prog, k)}
    }
    
    # rbind the list of results to get one massive data frame
    w = sapply(year_fire_sp, function(x){"SpatialPolygonsDataFrame" %in%  class(x)})
    year_fire = year_fire_sp[w] #remove the ones that were empty or broken
    year_fire_df = rbindlist(lapply(year_fire, function(x){x@data}), fill=T)
    year_fire = unlist(lapply(year_fire, function(x){x@polygons}))
    
    # convert to a SpatialPolygonsDataFrame
    row.names(year_fire_df) = year_fire_df$ID
    year_fire = SpatialPolygonsDataFrame(SpatialPolygons(year_fire), year_fire_df)
    crs(year_fire) = crs(data_ll)
    saveRDS(year_fire, paste0("data/fire/clusters_", y, ".RDS"))
  
    print(y)
}

###############################################################
# Process smoke
###############################################################

# read in data to start allocating plumes to fires
smoke_length = sapply(smoke, function(x){nrow(x) > 0})
smoke = smoke[smoke_length]
fire_comb = list()
fire_date = "2006"
fire_year = readRDS("data/fire/clusters_2006.RDS")
crs(fire_year) = crs(data_ll)


# runs slowly (1.5 hrs)
j = 1
k = 1
prog = txtProgressBar(min=0, max=length(smoke), initial=0, char="-", style=3)
for (i in 1:length(smoke)) {
    
    # get the date and relevant smoke
    date_str = names(smoke)[[i]]
    date = as.Date(date_str, format="%Y%m%d")
    x = smoke[[i]]
    
    # make sure you have the right fire file loaded. If you looped to a new year this 
    # will read in the fires for the new year
    if (substr(date_str, 1, 4) != fire_date) {
        fire_date = substr(date_str, 1, 4)
        fire_year = readRDS(paste0("data/fire/clusters_",  fire_date, ".RDS"))
        crs(fire_year) = crs(data_ll)
    }
    
    if(nrow(x) == 0) {next}
    
    # convert smoke data to SP
    x = tryCatch({
        as_Spatial(x)
    }, error = function(e) {
        x = x[!is.na(st_is_valid(x)),]
        x = st_cast(x, "POLYGON")
        as_Spatial(x)
    })
    x$lat = x$lon = x$f_area = NA
    
    # create a buffer around plumes so that it can capture fires slightly outside plume
    x_buf = gBuffer(x, byid=T, width=0.15) #buffer plumes a bit in case they don't overlap
    x_buf = spTransform(x_buf, crs(data_ll))
    f_buf = fire_year[fire_year$date == date_str, ]
    
    if(nrow(f_buf) != 0) {
        fires_over_all = over(x_buf, f_buf, returnList = T)
        
        #find the largest fire that each plume overlaps
        fires_over = unlist(lapply(fires_over_all, 
                                   FUN=function(y){if(nrow(y)>0){y[y$area == max(y$area),]$id[1]}else{NA}}))
        x[!is.na(fires_over), c("lat", "lon")]  = coordinates(f_buf[fires_over[!is.na(fires_over)],])
        x[!is.na(fires_over), "f_area"] = f_buf[fires_over[!is.na(fires_over)],]$area
    }
    
    #get smoke area
    x$area = sapply(slot(x, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area"))
    x = spTransform(x, crs(data_ll))
    
    #fix ids so that the smoke files can all be merged
    x$ID = j:(nrow(x)+j-1)
    row.names(x) = as.character(j:(length(x$ID)+j-1))
    j = length(x$ID)+j
    smoke[[i]] = x
    
    #print progress
    if (i %% 100 == 0) {setTxtProgressBar(prog, i)}
}


#rbind all the days  of smoke files and merge data back to the polygons
smoke_df = rbindlist(lapply(smoke, function(x){x@data}), fill=T)
smoke = unlist(lapply(smoke, function(x){x@polygons}))
smoke = SpatialPolygonsDataFrame(SpatialPolygons(smoke), smoke_df)
crs(smoke) = crs(data_ll)

#save the parsed smoke data
writeOGR(smoke, "data/clean/full_smoke_data", "full_smoke_data", 
         driver="ESRI Shapefile", overwrite_layer=T)

#get the plumes over each grid id for all time
num = over(data_ll, smoke, returnList=T) #~13 minutes
ids = unique(data_ll$id)
for (i in 1:length(ids)) {num[[i]]$id = ids[i]}

# loop through the list of smoke for each id and get distance to origin fire
# takes a while (±5 minutes)
prog = txtProgressBar(min=0, max=length(num), initial=0, char="-", style=3)
for (i in 1:length(num)) {
    cur = num[[i]]
    cur$fire_dist = NaN
    id_loc = coordinates(data_ll[data_ll$id==cur$id[1], ])
    
    #get distance between grid cell and origin fire of each plume
    nearest = FNN::get.knnx(cur[!is.na(cur$lat), c("lat", "lon")], id_loc, 
                            k=nrow(cur[!is.na(cur$lat),]))
    cur[!is.na(cur$lat), ]$fire_dist = nearest$nn.dist[order(nearest$nn.index)]
    
    num[[i]] = cur
    if (i %% 50 == 0) {setTxtProgressBar(prog, i)}
}

#combine the overlaps in to one big data table
num = rbindlist(num)
saveRDS(num, "data/clean/smoke_overlaps.RDS")

#define the bins for smoke size, fire size, and fire dist
#used to all divide by .33 and .66
med_smoke = quantile(num$area, 0.33, na.rm=T)  
large_smoke = quantile(num$area, 0.66, na.rm=T) 
med_dist = quantile(num$fire_dist, 0.33, na.rm=T) 
large_dist = quantile(num$fire_dist, 0.66, na.rm=T) 
med_fire = quantile(num$f_area, 0.33, na.rm=T) 
large_fire = quantile(num$f_area, 0.66, na.rm=T) 

#define the kmeans categories
num$scaled_f_area = scale(num$f_area)
num$scaled_area = scale(num$area)
num$scaled_fire_dist = scale(num$fire_dist)
km = kmeans(num[!is.na(num$f_area), c("scaled_f_area", "scaled_area", "scaled_fire_dist")], 
            8, nstart=5, iter.max=20)
num$category = 0
num[!is.na(num$f_area), ]$category = km$cluster
num$category[num$category==0] = NA

# good plot to investigate how plumes have been grouped 
#sub = sample.int(nrow(num), 1000)
#plot_ly(x=log(num$f_area[sub]+1), y=log(num$area[sub]+1), z=num$fire_dist[sub], 
#        type="scatter3d", mode="markers", color=as.factor(num$category[sub]), 
#        alpha=0.7, size=1) %>%
#    layout(scene = list(xaxis = list(title = "fire area"),
#                        yaxis = list(title = "smoke area"),
#                        zaxis = list(title = "distance to fire")))

#get linear weighted smoke vars
num$weightedlin_large = ifelse(num$area>large_smoke, 1, 0)/num$fire_dist
num$weightedlin_med = ifelse(num$area>med_smoke, 1, 0)/num$fire_dist - num$weightedlin_large

#get fire area weighted smoke vars
num$weightedlin_fire_large = (ifelse(num$area>large_smoke, 1, 0)/num$fire_dist) * (num$f_area/1000)
num$weightedlin_fire_med = ((ifelse(num$area>med_smoke, 1, 0)/num$fire_dist) * (num$f_area/1000)) - num$weightedlin_fire_large

#get capped linear weighted smoke vars
num$fire_dist_cap = num$fire_dist
num$fire_dist_cap[num$fire_dist_cap<1] = 1
num$weightedlin_cap_large = ifelse(num$area>large_smoke, 1, 0)/num$fire_dist_cap
num$weightedlin_cap_med = ifelse(num$area>med_smoke, 1, 0)/num$fire_dist_cap - num$weightedlin_cap_large

#get LOG weighted smoke vars
num$weightedlog_large = ifelse(num$area>large_smoke, 1, 0)/log(num$fire_dist_cap+1)
num$weightedlog_med = ifelse(num$area>med_smoke, 1, 0)/log(num$fire_dist_cap+1) - num$weightedlin_large

#summarize it by year
num$year = as.numeric(substr(num$date, 1, 4))
setkey(num, year, id) #set key by year to speed up 

d = num[,.(smoke=.N, 
           smoke_small=sum(area <= med_smoke, na.rm=T),
           smoke_med=sum(area > med_smoke & area <= large_smoke, na.rm=T),
           smoke_large=sum(area > large_smoke, na.rm=T),
           
           smoke_linear_large=sum(weightedlin_large, na.rm=T),
           smoke_linear_med=sum(weightedlin_med, na.rm=T),
           smoke_linearcap_large=sum(weightedlin_cap_large, na.rm=T),
           smoke_linearcap_med=sum(weightedlin_cap_med, na.rm=T),
           smoke_linearfire_large=sum(weightedlin_fire_large, na.rm=T),
           smoke_linearfire_med=sum(weightedlin_fire_med, na.rm=T),
           smoke_log_large=sum(weightedlog_large, na.rm=T),
           smoke_log_med=sum(weightedlog_med, na.rm=T),
           
           smoke_day=uniqueN(date), 
           smoke5=sum(Density%in%c("5.000",5)), 
           smoke16=sum(Density%in%c("16.000",16)), 
           smoke27=sum(Density%in%c("27.000",27)), 
           smoke_near=sum(fire_dist<=2, na.rm=T) + sum(is.na(fire_dist)),
           smoke_nearish=sum(fire_dist<=8, na.rm=T),
           smoke_farish=sum(fire_dist<=20, na.rm=T),
           smoke_far=sum(fire_dist>=20, na.rm=T),
           area=as.numeric(mean_na(area)),
           
           sfire_sdist = as.integer(sum_na(f_area<=med_fire & fire_dist<=med_dist)),
           sfire_mdist = as.integer(sum_na(f_area<=med_fire & fire_dist>med_dist & fire_dist<=large_dist)),
           sfire_ldist = as.integer(sum_na(f_area<=med_fire & fire_dist>large_dist)),
           mfire_sdist = as.integer(sum_na(f_area>med_fire & f_area<=large_fire & fire_dist<=med_dist)),
           mfire_mdist = as.integer(sum_na(f_area>med_fire & f_area<=large_fire & fire_dist>med_dist & fire_dist<=large_dist)),
           mfire_ldist = as.integer(sum_na(f_area>med_fire & f_area<=large_fire & fire_dist>large_dist)),
           lfire_sdist = as.integer(sum_na(f_area>large_fire & fire_dist<=med_dist)),
           lfire_mdist = as.integer(sum_na(f_area>large_fire & fire_dist>med_dist & fire_dist<=large_dist)),
           lfire_ldist = as.integer(sum_na(f_area>large_fire & fire_dist>large_dist)),
           
           cat1 = as.integer(sum_na(category==1)),
           cat2 = as.integer(sum_na(category==2)),
           cat3 = as.integer(sum_na(category==3)),
           cat4 = as.integer(sum_na(category==4)),
           cat5 = as.integer(sum_na(category==5)),
           cat6 = as.integer(sum_na(category==6)),
           cat7 = as.integer(sum_na(category==7)),
           cat8 = as.integer(sum_na(category==8)),
           
           fire_dist = mean(fire_dist, na.rm=T)),
        by=list(id, year)]
d$year = as.character(d$year)
d = merge(d, expand.grid(id=unique(data_ll$id), year=as.character(years)), by=c("id", "year"), all=T)
d[is.na(d)] = 0

saveRDS(d, "data/model_inputs/smoke_processed.RDS")
