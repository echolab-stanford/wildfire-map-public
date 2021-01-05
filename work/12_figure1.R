# written by MB, edited to add WUI and smoke pm by AD on 02/10/2020

source("work/00_functions.R")

######################################################################################

# plot burned area and suppression costs
  # data are from national interagency fire center
  #  https://www.nifc.gov/fireInfo/fireInfo_documents/SuppCosts.pdf
dt <- read_csv('data/emissions/supression_costs.csv')
dt$burned_area <- (dt$burned_area*0.404686)/1e6  #put in millions

#deflate cost time series
  # GDP deflator data are from SL FED:  https://fred.stlouisfed.org/series/USAGDPDEFAISMEI
def <- read_csv('data/emissions/USAGDPDEFAISMEI.csv')
def$year <- as.numeric(substr(def$DATE,1,4))
def$deflator <- def$USAGDPDEFAISMEI/def$USAGDPDEFAISMEI[def$year==2015]
dt <- left_join(dt,def,by="year")
dt$total_def <- dt$total/dt$deflator/1e9  #deflated billions of dollars

######################################################################################

# trends in PM2.5 from EPA: https://www.epa.gov/air-trends/particulate-matter-pm25-trends
f <- list.files('data/EPA_trend/')
epa <- data.frame(Year=2000:2018)
for (fl in f) {
  dta <- read_csv(paste0('data/EPA_trend/',fl)) 
  dta <- dta[,1:2]
  names(dta)[2] <- unlist(strsplit(unlist(strsplit(fl,".csv")),"PM25"))[2]
  epa <- left_join(epa,dta)
}
names(epa)[1] <- 'year'
dt <- left_join(dt,epa[,c('year',"National")])

######################################################################################

# abatzoglou fuel aridity data - we are using the mean across 8 aridity time series, as supplied in their SI data
fuel <- read_csv('data/emissions/abatzoglou_data.csv',skip=10)
fuel$mean_aridity <- fuel$`Z MEAN`
dt <- left_join(dt,fuel[,c("year","mean_aridity")])


######################################################################################

#wui data, read in years then combine 
w = read.csv("data/WUI/2006_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2006)
w1 = read.csv("data/WUI/2008_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2008)
w2 = read.csv("data/WUI/2011_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2011)
w3 = read.csv("data/WUI/2013_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2013)
w4 = read.csv("data/WUI/2016_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2016)
w5 = read.csv("data/WUI/2001_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2001)
w6 = read.csv("data/WUI/2004_combined_wui_hh_data.csv") %>% summarise(wui = sum(total_wui_households)/1000000, year=2004)
wui  = rbind(w, w1, w2, w3, w4, w5, w6)
dt = merge(dt,  wui, all.x=T)
dt$wui = na.approx(dt$wui,na.rm="FALSE")

######################################################################################

#read in the full modeling data with results
results = readRDS("data/clean/results_all.RDS")
regions = data.frame(state=as.character(unique(results$state)), stringsAsFactors=F)
regions$epa_region = c(10, 10, 8, 8, 5, 5, 5, 1, 10, 8, 1, 8, 2, 1, 7, 7, 5, 1, 3, 9, 9,  
                       8, 5, 5, 1, 2, 1, 8, 7, 7, 3, 3, 3, 3, 4, 9, 6, 6, 4, 4, 6, 6, 4, 
                       4, 4, 4, 6, 4)
full_data = merge(results, regions, by="state")
full_data$year = as.numeric(as.character(full_data$year))

#aggregate pm data to national level
full_nat = full_data %>% group_by(year) %>% 
    dplyr::summarise(pm = mean(preds), pm_nosmoke=mean(preds-diff), 
              pm_smoke=mean(diff), smoke_days = mean(smoke_day), 
              perc=mean(diff/preds, na.rm=T))
full_data = full_data %>% group_by(year, epa_region) %>% 
    dplyr::summarise(pm=mean(preds), pm_nosmoke=mean(preds-diff), 
              pm_smoke=mean(diff), smoke_days=mean(smoke_day), 
              perc=mean(diff/preds, na.rm=T))

######################################################################################

#load in fire perimeter data
prescribed = read.csv("data/fire/prescribed_burn_acres.csv") %>% 
    dplyr::filter(region != "AK") 
prescribed$region[prescribed$region=="SO"] = "NO"  #group two californias to ~match EPA
prescribed = prescribed %>%  group_by(year, region) %>%
    dplyr::summarise(acres = sum(acres))
prescribed$hectares = (prescribed$acres*0.404686)/100000 #acres to mil hectares

######################################################################################

#load in improve data for organic carbon plot
imp1 = read.csv("data/improve/IMPROVE_1988-2006.txt", na.strings="-999")[, c("Date", "State", "SiteCode", "OCf.Value")]
imp2 = read.csv("data/improve/IMPROVE_2007.txt", na.strings="-999")[, c("Date", "State", "SiteCode", "OCf.Value")]
imp3 = read.csv("data/improve/IMPROVE_2008-2016.txt", na.strings="-999")[, c("Date", "State", "SiteCode", "OCf.Value")]
imp4 = read.csv("data/improve/IMPROVE_2017.txt", na.strings="-999")[, c("Date", "State", "SiteCode", "OCf.Value")]
imp5 = read.csv("data/improve/IMPROVE_2018.txt", na.strings="-999")[, c("Date", "State", "SiteCode", "OCf.Value")]
imp = rbind(imp1, imp2, imp3, imp4, imp5)
imp = merge(imp, regions, by.x="State", by.y="state")
imp$month = substr(as.character(imp$Date), 1, 2)
imp$year = as.numeric(substr(as.character(imp$Date), 7, 10))

#get summary stats for improve
nat_imp = imp %>% 
    dplyr::filter(month %in% c("07", "08", "09")) %>% 
    group_by(Date, SiteCode, month, year) %>%
    dplyr::summarise(carbon=mean(OCf.Value, na.rm=T)) %>%
    group_by(month, SiteCode, year) %>%
    dplyr::summarise(carbon=mean(carbon, na.rm=T)) %>%
    group_by(year) %>%  
    dplyr::summarise(carbon=mean(carbon, na.rm=T))

imp = imp %>% 
    dplyr::filter(month %in% c("07", "08", "09")) %>% 
    group_by(Date, SiteCode, month, year, epa_region) %>%
    dplyr::summarise(carbon=mean(OCf.Value, na.rm=T)) %>%
    group_by(month, SiteCode, year, epa_region) %>%
    dplyr::summarise(carbon=mean(carbon, na.rm=T)) %>%
    group_by(year, epa_region) %>%  
    dplyr::summarise(carbon=mean(carbon, na.rm=T))


######################################################################################
# functions
######################################################################################

# add trend line to existing plot
addtrend <- function(x,custyr=F,years=2000:2018, data=dt,lty=1) {
  yrs <- data$year[is.na(data[,x])==F]
  if (custyr==T) {yrs=years}
  fmla <- as.formula(paste0(x,"~ year"))
  mod <- lm(fmla,data=data[data$year%in%yrs,])
  print(summary(mod)$coefficients['year',4])  #print p-value on trend
  cf <- coef(mod)
  yy=cf[1]+cf[2]*yrs
  lines(yrs,yy,col="red",lty=lty)
}

# plot both time series and estimated trends
toplot <- function(x,main,ylab,label,ylim,xaxis=NULL,col="red") {
  plot(dt[,c("year",x)],type="l",col="grey",ylab="",axes=F,ylim=ylim,xlab="")
  axis(2,las=1)
  if (!is.null(xaxis)) {axis(1,at=c(seq(xaxis,2015,5),2018))}
  yrs <- dt$year[is.na(dt[,x])==F]
  fmla <- as.formula(paste0(x,"~ year"))
  mod <- lm(fmla,data=dt)
  print(summary(mod)$coefficients['year',4])  #print p-value on trend
  cf <- coef(mod)
  yy=cf[1]+cf[2]*yrs
  lines(yrs,yy,col=col)
  mtext(main, side=3,line=-1,adj=0,cex=0.8)
  mtext(paste0("trend = ",round(cf[2],2),label), side=3,line=-2,adj=0,cex=0.7)
  title(ylab=ylab, line=2.2, cex.lab=1.2)
}



######################################################################################
# MAKE THE PLOT
######################################################################################

pdf(file="images/Final/Figure1.pdf",width=7,height=8,useDingbats = F)
#layout(matrix(1:8, ncol=2), width=scale(c(2018-1985, 2018-2000), center=F))
par(mar=c(3,4,0,1),mfcol=c(4,2))

######################################################################################
# burned area

toplot("burned_area","Burned area",ylim=c(0,4.5),ylab="million hectares",label=" million hectares/year")

######################################################################################
# fuel aridity

toplot("mean_aridity","Fuel aridity",ylim=c(-1,1.5),ylab="aridity (sd)",label=" sd/year")

######################################################################################
# wui

toplot("wui","Houses in WUI",ylim=c(42.5, 49.5),ylab="million houses",
       label=" million houses/year", col="blue")

######################################################################################
# suppression costs

toplot("total_def","Suppression costs",ylim=c(0,3),ylab=" costs ($ billion)",
       label="$ billion/year", xaxis=1985)

######################################################################################
# acres of prescribed burns
 
prescribed = prescribed[prescribed$year <= 2018, ]

plot(1,type="n",xlim=c(1985,2018),ylim=c(0, 25),las=1, ylab="",axes=F,xlab="")
title(ylab="hectares prescribed (100k's)", line=2.2, cex.lab=1.2)
axis(2,las=1)
for (i in c("EA", "GB", "NR", "RM", "SA", "SW", "NO", "NW")) {
    f = prescribed[prescribed$region==i, ]
    if (i %in% c("NW", "NO")) {
        if (i == "NW") {nm="Northwest"} else if (i == "NO") {nm = "West"} else {nm=""}
        col="grey20"
        text(2017,f[f$year==2018,"hectares"],nm,cex=0.4)
        lines(f$year,f$hectares,col=col)
    } else {
        col="grey"
        lines(f$year,f$hectares,col=col)
    }
}

yrs = 1998:2018
cf1 <- coef(lm("hectares ~ as.numeric(year)",data=prescribed[prescribed$region=="SA",]))
yy=cf1[1]+cf1[2]*yrs
lines(yrs,yy,col="red")

cf <- coef(lm("hectares ~ as.numeric(year)",data=prescribed[prescribed$region!="SA",]))
yy=cf[1]+cf[2]*yrs
lines(yrs,yy,col="red")

mtext("Hectares of prescribed burn", side=3,line=-1,adj=0,cex=0.8)
mtext(paste0("trend = ", round(cf1[2]*100000,0), " hectares/year in the South; ", 
             round(cf[2]*100000,0), " everywhere else"), 
      side=3,line=-2,adj=0,cex=0.7)


######################################################################################
# smoke days by year

plot(1,type="n",xlim=c(1985,2018),ylim=c(0, 100),las=1, ylab="",axes=F,xlab="")
title(ylab="smoke days", line=2.2, cex.lab=1.2)
axis(2,las=1)
for (i in 1:10) {
    f = full_data[full_data$epa_region == i, ]
    if (i %in% 9:10) {
        if (i == 10) {nm="Northwest"}
        if (i == 9) {nm = "West"}
        col="grey20"
        text(2017,f[f$year==2018,"smoke_days"],nm,cex=0.4)
        lines(f$year,f$smoke_days,col=col)
    } else {
        col="grey"
        lines(f$year,f$smoke_days,col=col)
    }
}

yrs = 2006:2018
cf <- coef(lm("smoke_days ~ as.numeric(year)",data=full_nat))
yy=cf[1]+cf[2]*yrs
lines(yrs,yy,col="blue")
mtext("Smoke days per year", side=3,line=-1,adj=0,cex=0.8)
mtext(paste0("trend = ",round(cf[2],2)," more days/year"), side=3,line=-2,adj=0,cex=0.7)

######################################################################################
# plot PM over time, by epa region - 2006-2018

plot(1,type="n",xlim=c(1985,2018),ylim=c(5, 18),las=1,ylab="",axes=F,xlab="")
title(ylab="PM2.5", line=2.2, cex.lab=1.2)
axis(2,las=1)
for (i in 2:dim(epa)[2]) {
    nm <- names(epa)[i]
    if (nm == "National") {next}
    if (nm %in% c('West','Northwest')) {
        col="grey20"
        text(2017,epa[epa$year==2018,i],nm,cex=0.4)
        lines(epa$year,epa[,i],col=col)
    } else {
        col="grey"
        lines(epa$year,epa[,i],col=col)
    }
}

#add the red line for average
x = "National"
yrs <- dt$year[is.na(dt[,x])==F]
fmla <- as.formula(paste0(x,"~ year"))
cf <- coef(lm(fmla,data=dt))
yy=cf[1]+cf[2]*yrs
lines(yrs,yy,col="red")
mtext("PM2.5", side=3,line=-1,adj=0,cex=0.8)
mtext(paste0("trend = ",round(cf[2],2)," ug/m3/year"), side=3,line=-2,adj=0,cex=0.7)


######################################################################################
# PM FROM SMOKE

plot(1,type="n",xlim=c(1985,2018),ylim=c(0,100),las=1, ylab="",axes=F,xlab="")
title(ylab="% PM2.5 from smoke", line=2.2, cex.lab=1.2)
axis(2,las=1)
axis(1,at=c(seq(1985,2015,5),2018))
for (i in 1:10) {
    f = full_data[full_data$epa_region == i, ]
    if (i %in% 9:10) {
        if (i == 10) {nm="Northwest"}
        if (i == 9) {nm = "West"}
        col="grey20"
        text(2017,f[f$year==2018,"perc"]*100,nm,cex=0.4)
        lines(f$year,f$perc*100,col=col)
    } else {
        col="grey"
        lines(f$year,f$perc*100,col=col)
    }
}

yrs = 2006:2018
cf <- coef(lm("perc*100 ~ as.numeric(year)",data=full_nat))
yy=cf[1]+cf[2]*yrs
lines(yrs,yy,col="blue")
mtext("Predicted % of PM2.5 from wildfire smoke", side=3,line=-1,adj=0,cex=0.8)
mtext(paste0("trend = ",round(cf[2],2)," percentage points/year"), side=3,line=-2,adj=0,cex=0.7)

dev.off() 