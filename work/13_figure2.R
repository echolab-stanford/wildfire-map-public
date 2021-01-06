source("work/00_functions.R")


######################################################################################
# Read in the data
######################################################################################

#read in predictions
predictions = readRDS("data/clean/results_all.RDS")
predictions = predictions %>% group_by(id) %>% 
    dplyr::summarize(pop = mean(pop), perc_diff = mean(perc_diff), 
              diff = mean(diff), preds = mean(preds))

#attribute PM predictions to counties
data_ll = readRDS("data/clean/national_grid.RDS")
counties = readRDS("data/boundaries/counties.RDS")
data_ll_geo = fortify(data_ll, region="id") 
counties_geo = fortify(counties, region="GEOID") 

c = over(counties, data_ll, returnList=T)
for (i in  1:nrow(counties)) {
    c[[i]]$fips = counties$GEOID[i]
}
c = rbindlist(c)[, c("id", "fips")]
c = merge(c, predictions, all=T)
c = c %>% group_by(fips) %>% 
    dplyr::summarize(perc_diff = mean(perc_diff), 
                     diff = mean(diff), preds = mean(preds))

#merge it to the acs  data
acs = read_rds('data/census/us_county_average_pm25_income_race.rds')
acs$fips  =  str_pad(acs$fips, 5, "left", "0")
c = merge(c, acs, by="fips", all=T)
c = c %>% mutate(pop = n_white + n_black + n_asian + n_hisp, perc_white = n_white/pop)

#get mappable data_ll
data_ll = spTransform(data_ll, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
data_ll = data_ll[, c("id", "state")]
data_ll = st_as_sf(data_ll)


######################################################################################
# Make plots
######################################################################################

######################################################################################
# A

predictions = readRDS("data/clean/results_all.RDS")
p = predictions %>% dplyr::filter(year %in% 2006:2008) %>% 
    group_by(id) %>% 
    dplyr::summarise(perc_diff = mean(perc_diff), preds = mean(preds),  
                               diff=mean(diff))
p2 = predictions %>% dplyr::filter(year %in% 2016:2018) %>% 
    group_by(id) %>% dplyr::summarise(perc_diff = mean(perc_diff), preds = mean(preds),  
                               diff=mean(diff))

cols = get_col_regions()[seq(1,100,floor(100/8))]
diff = getJenksBreaks(c(p$diff, p2$diff), 9)
diff[1] =  0
lim  = c(-0.04, max(diff))

#this may not work every time due to known bug: https://github.com/tidyverse/ggplot2/issues/2252
data = merge(data_ll_geo, p, by="id")
pdf(file='images/Raw/Figure2a_percsmoke.pdf',width=10,height=10)
ggplot(data, aes(long,lat,group=group, fill=diff)) + # the data
    geom_polygon() + # make polygons
    scale_fill_gradientn(limits=lim, colors=cols, values=scales::rescale(diff)) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(data$lat)) + labs(fill="")
dev.off()

data = merge(data_ll_geo, p2, by="id")
pdf('images/Raw/Figure2b_percsmoke.pdf',width=10,height=10)
ggplot(data, aes(long,lat,group=group, fill=diff)) + # the data
    geom_polygon() + # make polygons
    scale_fill_gradientn(limits=lim, colors=cols, values=scales::rescale(diff)) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(data$lat)) + labs(fill="")
dev.off()

######################################################################################
# C

regions = c("NorthEast", "MidAtlantic", "SouthEast", "MidWest", "SouthernPlains", 
            "GreatPlains", "RockyMountains", "SouthWest", "NorthWest")
smokehr = c()

#loop through the data from Brey to find origin of the smoke in each region
for (r in regions) {
    load(paste0('data/Brey/',r,
                '_smokeHourSummary_2007-2014_season=6-9_gdas1_2day+6day_daily_plume=TRUE.RData'))  #gets loaded in as "dataList"
    smokehr = rbind(smokehr,apply(dataList[[1]],2,sum))
}

tot = apply(smokehr,1,sum) #sum up smoke hours by region of experience
smokehr = data.frame(region=regions,smokehr,total=tot)
us = names(smokehr)%in%regions
smokehr$us = apply(smokehr[,us],1,sum) # smoke from US regions
smokehr$pctoutside = 1- smokehr$us/smokehr$total  #smoke from outside the US

brey_region_mapping = read.csv('data/Brey/BreyStatesRegions.csv')
brey_region_mapping$fips = as.factor(str_pad(brey_region_mapping$fips, 2, "left", "0"))

shp = readRDS('data/boundaries/counties.RDS')
shp = sp::merge(shp, brey_region_mapping, by.x="STATEFP", by.y="fips", duplicateGeoms = T)
shp = clgeo_Clean(shp)
shp = gBuffer(shp, width=.03, byid=T)
groups = aggregate(shp, by = "region")
g = SpatialPolygonsDataFrame(gSimplify(groups, 0.01), data=groups@data)
x = fortify(g, region="region")
x = merge(x, smokehr[, c('region', 'total', 'pctoutside')], by.y="region", by.x="id")
x$pct_cut = cut(x$pctoutside, seq(0, 1, .1))

pdf('images/Raw/Figure2d_pctoutside.pdf',width=5,height=5)
ggplot(x, aes(long, lat, group=group,fill=pct_cut)) + # the data
    geom_polygon(aes(long, lat, group=group), fill=NA, color="black", lwd=.5, data=x) +
    geom_polygon() + # make polygons
    scale_fill_grey(start=1, end=0, drop=F) +
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(x$lat)) 
dev.off()


from_west = merge(x, smokehr[, c("region", "NorthWest", "SouthWest", "total")], 
                  by.y="region", by.x="id")
from_west$west_perc = (from_west$NorthWest + from_west$SouthWest)/from_west$total.x
from_west$west_perc = cut(from_west$west_perc, seq(0, 1, .1))

pdf('images/Raw/Figure2c_crossboundaries.pdf',width=5,height=5)
ggplot(from_west, aes(long, lat, group=group,fill=west_perc)) + # the data
    geom_polygon(aes(long, lat, group=group), fill=NA, color="black", lwd=.5, data=x) +
    geom_polygon() + # make polygons
    scale_fill_grey(start=1, end=0, drop=F) +
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(x$lat)) 
dev.off()



######################################################################################
# D

pdf('images/Raw/Figure2e_racepm.pdf',width=4,height=3.5)
ggplot(c, aes(perc_white, diff)) + 
    stat_bin_2d(bins=30) +
    scale_fill_gradientn(colors = brewer.pal(9,'Greys')[3:9], guide=F) +
    geom_smooth(method='lm', formula=y~x, color="dodgerblue4", size=1.5, se=F,fullrange=T) +
    ylab('predicted PM2.5 from smoke  (ug/m3)') + scale_y_continuous(limits = c(0,2.7)) + 
    xlab('Percent white') + scale_x_continuous(limits = c(0,1)) +
    theme_anne(font="sans", size=13)
dev.off()

######################################################################################
# E

pdf('images/Raw/Figure2f_racesmoke.pdf',width=4,height=3.5)
ggplot(c, aes(perc_white, preds)) + 
    stat_bin_2d(bins=30) +
    scale_fill_gradientn(colors = brewer.pal(9,'Greys')[3:9], guide=F) +
    geom_smooth(method='lm', formula=y~x, color="dodgerblue4", size=1.5, se=F,fullrange=T) +
    ylab('predicted PM2.5 (ug/m3)') + scale_y_continuous(limits = c(0,13.5)) + 
    xlab('Percent white') + scale_x_continuous(limits = c(0,1)) +
    theme_anne(font="sans", size=13)
dev.off()
