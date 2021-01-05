source("work/00_functions.R")


########  CROSSWALKS + READ IN

#modeling results
results_all = readRDS("data/clean/results_all.RDS")
data_ll = readRDS("data/clean/national_grid.RDS")
data_ll_geo = fortify(data_ll, region="id")

states = readOGR("data/boundaries/tl_2019_us_state", "tl_2019_us_state")
states = states[!states$STATEFP %in% c("02", "60", "78", "72", "69", "66", "60", "15"),]
states = gSimplify(states, 0.01)
states = fortify(states, region="STATEFP")

#epa region crosswalk
regions = data.frame(state=as.character(unique(as.character(results_all$state))), stringsAsFactors=F)
regions$epa_region = c(10, 10, 8, 8, 5, 5, 5, 1, 10, 8, 1, 8, 2, 1, 7, 7, 5, 1, 3, 9, 9,  
                       8, 5, 5, 1, 2, 1, 8, 7, 7, 3, 3, 3, 3, 4, 9, 6, 6, 4, 4, 6, 6, 4, 
                       4, 4, 4, 6, 4)
fips = read.csv("data/boundaries/state_fips_codes.csv")[, 1:2]
regions = merge(regions, fips, by="state")
regions$fips  = str_pad(regions$fips, 2, "left", "0")

#boundaries
area = readRDS("data/boundaries/counties.RDS")[, c("GEOID",  "ALAND", "AWATER")]
area$acres = (area$ALAND+area$AWATER)*0.000247105381 #convert to acres
counties = read.csv("data/boundaries/state_fips_codes.csv")[, c("state", "fips")]
counties$fips = str_pad(counties$fips, 2, "left", "0")

######### NEI READ  IN

nei2008 = read.csv("data/clean/2008_by_county_with_wildfire.csv")[, c(1, 6)]
nei2008$year = 2008
nei2011 = read.csv("data/clean/2011_by_county_with_wildfire.csv")[, c(1, 6)]
nei2011$year = 2011
nei2014 = read.csv("data/clean/2014_by_county_with_wildfire.csv")[, c(1, 6)]
nei2014$year = 2014
nei2017 = read.csv("data/clean/2017_by_county_with_wildfire.csv")[, c(2, 6)]
nei2017$year = 2017
names(nei2017) = c("state_and_county_fips_code", "PM25.PRI_wildfire", "year")
nei = rbindlist(list(nei2008, nei2011, nei2014, nei2017))
nei$state_and_county_fips_code = stringr::str_pad(as.character(nei$state_and_county_fips_code), 
                                                  5,  "left", "0")
nei$state = substr(nei$state_and_county_fips_code,  1, 2)
nei = merge_verbose(nei, area, by.x="state_and_county_fips_code",  by.y="GEOID")
nei = merge_verbose(nei, counties, by.x="state", by.y="fips")
nei = merge_verbose(nei, regions, by.x="state.y", by.y="state", all.x=T)
names(nei)[4] = "PM25_fire"

# get county level estimates of diff
area = clgeo_Clean(area)
x = over(area, data_ll, returnList=T)
FUN = function(x) {results_all %>% 
        dplyr::filter(id %in% x$id) %>% 
        group_by(year) %>% 
        dplyr::summarise(diff = mean(diff)) }
x = lapply(x, FUN)
for (i in 1:length(x)) {x[[i]]$GEOID = area[i,]$GEOID}
x = rbindlist(x)
x$year = as.integer(as.character(x$year))
x = merge(nei[, c("state_and_county_fips_code", "year", "PM25_fire", "acres", "epa_region", "state")], x, 
          by.x= c("state_and_county_fips_code", "year"), by.y=c("GEOID", "year"))
all = x %>% dplyr::group_by(state, year) %>% 
    dplyr::summarise(diff = mean(diff, na.rm=T), PM_scaled = mean(PM25_fire/acres, na.rm=T))
all$PM_scaled = all$PM_scaled*1000

names(all) = c("state", "year", "PM2.5 from smoke (model)", 
               "PM2.5 from smoke (NEI)")

pdf('images/FigureS5_modeltoNEI.pdf',width=5, height=5, useDingbats=F)
options(scipen=0)
panel(all[, 3:4], font="sans", size=26, a=0.8) 
dev.off()

######  ODELL

odell = data.frame(year=as.character(2006:2016), 
                   smoke= c(3.5, 3.7, 4.3, 2.5, 2.2, 2.1, 4.8, 3.7, 3.5, 5.7, 2.2),
                   full = c(8.1, 7.8, 7.8, 6.8, 5.8, 6.1, 9.1, 8.1, 7, 9, 5.5), 
                   type="odell", stringsAsFactors=F)
odell$perc = odell$smoke / ((odell$full-odell$smoke)*3 + odell$full) * 100
r  = results_all %>% 
    dplyr::filter(state %in% c("OR", "WA", "ID",  "MT")) %>% group_by(year) %>% 
    dplyr::summarise(perc=mean(perc_diff)*100) %>% mutate(type="model")
odell = rbind(odell[, c("year", "perc", "type")], r)
odell$year = as.numeric(as.character(odell$year))
odell$perc = as.numeric(as.character(odell$perc))

pdf('images/FigureS5_modeltoOdell.pdf',width=5,height=5, useDingbats=F)
ggplot(odell) + geom_line(aes(year, perc, group=type, color=type), size=3) + 
    theme_anne(font="sans", size=26) + xlab("Year") + 
    ylab("% PM2.5 from wildfire") + labs(color="") + ylim(0, 60) + 
    guides(color=F)
dev.off()


######  JAFFE

reg1 = Polygons(list(Polygon(cbind(c(-108.99, -116.8, -116.8, -108.99, -108.99), 
                                   c(41.01, 41.01, 49.13, 49.13, 41.01)))), ID=1)
reg2 = Polygons(list(Polygon(cbind(c(-102.2, -114.09, -114.09, -102.2, -102.2), 
                                   c(37.05, 37.05, 41.01, 41.01, 37.05)))), ID=2)
reg3 = Polygons(list(Polygon(cbind(c(-102.2, -114.09, -114.09, -102.2, -102.2), 
                                   c(30.57, 30.57, 37.05, 37.05, 30.57)))), ID=3)
reg4 = Polygons(list(Polygon(cbind(c(-118.0, -124.82, -124.82, -118.0, -118.0), 
                                   c(36.26, 36.26, 42.02, 42.02, 36.26)))), ID=4)
reg5 = Polygons(list(Polygon(cbind(c(-116.8, -125.24, -125.24, -116.8, -116.8), 
                                   c(42.02, 42.02, 49.13, 49.13, 42.02)))), ID=5)
jaffe_reg = SpatialPolygons( list(reg1, reg2, reg3, reg4, reg5) )
crs(jaffe_reg) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

x = over(jaffe_reg, data_ll, returnList = T)
for (i in 1:length(x)) {x[[i]] = as.data.frame(cbind(x[[i]]$id, i))}
x =  rbindlist(x)
names(x) = c("id", "jaffe_reg")
r = merge(results_all, x, by="id")
r = r %>% group_by(jaffe_reg) %>% dplyr::summarise(diff = mean(diff))
r$jaffe_diff = c(1.84, 1.09, 0.61, 0.81, 1.21)
#r$jaffe_diff = c(3.14, 2.26, 1.24, 1.28, 2.44)
r = melt(r,  "jaffe_reg")

pdf('images/FigureS5_modeltoJaffe.pdf',width=5,height=5, useDingbats=F)
ggplot(r, aes(x=jaffe_reg, y=value, fill=variable)) +
    geom_bar(stat="identity", width=.5, position = "dodge") + 
    xlab("Region") + ylab("PM2.5 from wildfires") + labs(fill="") + theme_anne(font="sans", size=26)  +  
    guides(fill=F)
dev.off()


######  FANN NUMBERS

r = results_all %>% dplyr::filter(year %in% 2008:2012) %>% group_by(year) %>%
    dplyr::summarise(mean = mean(diff))
r$type = "model"
r = rbind(data.frame(year=2008:2012, mean=c(1.1, 0.56, 0.7, 0.8, 0.92), type="Fann"), r)
r$type = factor(r$type, levels = c("model", "Fann"))

pdf('images/FigureS5_FannCompTime.pdf',width=5,height=5, useDingbats=F)
ggplot(r) + geom_line(aes(year, mean, color=type, group=type), size=3) + 
    theme_anne(font="sans", size=26) + guides(color=F) + 
    ylab("Mean PM2.5 from wildfire") + xlab("Year")
dev.off()


############################################################
# map - Fann
############################################################

breaks = c(-0.0001, 0.390, .86, 1.54, 2.44, 3.95, 6.55, 10.7, 16.8, 26.5, 49.3)
colors = c("#FEECD6", "#F7D8BC", "#F1C0A0", "#EBA987", "#E49171", "#DF775C", 
           "#D35C42", "#D24031", "#CA231B", "#C50A09")
results_all$diff_cat = cut(results_all$diff, breaks, 
                       labels=c("0 - 0.39", "0.39 - 0.86", "0.86 - 1.54", "1.54 - 2.44",
                                "2.44 - 3.95", "3.95 - 6.55", "6.55 - 10.7", "10.7 - 16.8", 
                                "16.8 - 26.5", "26.5 - 49.3"))

data = merge(data_ll_geo, results_all[results_all$year == 2008, ])
data = data[order(data$group, data$order), ]
eight = ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    geom_polygon(aes(long, lat, group=group), fill=NA, color="darkgrey", lwd=5, data=states) +
    coord_map("bonne", mean(data$lat)) + guides(fill=F) +
    xlim(-124.83,-58.62) + ylim(24.17, 49.38)

data = merge(data_ll_geo, results_all[results_all$year == 2009, ])
data = data[order(data$group, data$order), ]
nine = ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    geom_polygon(aes(long, lat, group=group), fill=NA, color="darkgrey", lwd=5, data=states) +
    coord_map("bonne", mean(data$lat)) + guides(fill=F) +
    xlim(-124.83,-58.62) + ylim(24.17, 49.38)

data = merge(data_ll_geo, results_all[results_all$year == 2010, ])
data = data[order(data$group, data$order), ]
ten = ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    geom_polygon(aes(long, lat, group=group), fill=NA, color="darkgrey", lwd=5, data=states) +
    coord_map("bonne", mean(data$lat)) + guides(fill=F) +
    xlim(-124.83,-58.62) + ylim(24.17, 49.38)

data = merge(data_ll_geo, results_all[results_all$year == 2011, ])
data = data[order(data$group, data$order), ]
eleven = ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    geom_polygon(aes(long, lat, group=group), fill=NA, color="darkgrey", lwd=5, data=states) +
    coord_map("bonne", mean(data$lat)) + guides(fill=F) +
    xlim(-124.83,-58.62) + ylim(24.17, 49.38)

data = merge(data_ll_geo, results_all[results_all$year == 2012, ])
data = data[order(data$group, data$order), ]
twelve = ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    geom_polygon(aes(long, lat, group=group), fill=NA, color="darkgrey", lwd=5, data=states) +
    coord_map("bonne", mean(data$lat)) + guides(fill=F) +
    xlim(-124.83,-58.62) + ylim(24.17, 49.38)

pdf("images/FigureS5-Fann2008.pdf",width=300,height=90, useDingbats=F)
eight
dev.off()

pdf("images/FigureS5-Fann2009.pdf",width=300,height=90, useDingbats=F)
nine
dev.off()

pdf("images/FigureS5-Fann2010.pdf",width=300,height=90, useDingbats=F)
ten
dev.off()

pdf("images/FigureS5-Fann2011.pdf",width=300,height=90, useDingbats=F)
eleven
dev.off()

pdf("images/FigureS5-Fann2012.pdf",width=300,height=90, useDingbats=F)
twelve
dev.off()

pdf(file="images/FigureS5-guide.pdf")
ggplot(data, aes(long, lat, group=group, fill=diff_cat)) + # the data
    geom_polygon(aes(fill=diff_cat)) + # make polygons
    scale_fill_manual(values=colors) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(data$lat))
dev.off()


############################################################
# map - Wilkins
############################################################

w = results_all %>%
    mutate(year = as.numeric(as.character(year))) %>%
    dplyr::filter(year >= 2008 & year <=  2012) %>%
    group_by(id) %>%
    summarise_at(vars(pm:preds0), mean, na.rm = T) %>%
    ungroup() %>%
    summarise(pm = mean(pm),
              preds = mean(preds),
              preds0 = mean(preds0))

perc = (w$preds - w$preds0)/w$preds
df = data.frame(value=c(perc*100, 10.5), variable=c("Our model", "Wilkins et al., 2018"))

pdf('images/FigureS5-Wilkenscheck.pdf',width=5,height=5, useDingbats=F)
ggplot(df, aes(y=value, x=variable, fill=variable)) +
    geom_bar(stat="identity", width=.5, position = "dodge") + 
    xlab("") + ylab("% of PM2.5 from wildfires") + labs(fill="") + 
    theme_anne(font="sans", size=26)  +  guides(fill=F) + ylim(0, 15)
dev.off()