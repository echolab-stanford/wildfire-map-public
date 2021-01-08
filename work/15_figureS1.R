source("work/00_functions.R")

# set plot theme
theme_set(theme_minimal())

# USA Borders
data(countriesLow)
borders = countriesLow[countriesLow$ADMIN=="United States of America" & countriesLow$TYPE!="Dependency",]


# --------------------------------------------------------
# Load data
# --------------------------------------------------------
load("data/improve/Summarized_IMPROVE_Data.Rdata")
load("data/improve/Processed_IMPROVE.Rdata")


# Full Time Series - Monthly
p = ggplot(data=regdat,aes(x=date,color=Region,group=Region)) + geom_line(aes(y=PM2.5),lty=3) + geom_line(aes(y=OC))
p = p + theme(legend.position = "none") + facet_grid(Region ~ .) + ylab(expression(paste("Organic Carbon and Total ",PM[2.5]," [",mu,"g ",m^-3,"]",sep=""))) + xlab("")
p = p + scale_color_brewer(palette="Dark2")

pdf(width=8,height=6,file="images/Raw/FigureS1c_IMPROVE_PM_and_OC_TimeSeries.pdf")
par(mar=c(5,4,1,1))
print(p)
dev.off()



# Full Time Series - Season Trends
imdat$Season = "Winter"
imdat$Season[as.numeric(imdat$Month) %in% c(3:6)] = "Spring"
imdat$Season[as.numeric(imdat$Month) %in% c(7:10)] = "Summer"
trends = imdat %>% group_by(Region,Season) %>% do(model=lm(OCf.Value~as.numeric(Year),data=.))

p = ggplot(data=imdat,aes(x=date,color=Region,group=Region)) + geom_smooth(aes(y=OCf.Value),method=lm,fill="lightgrey")
p = p + theme(legend.position = "none",plot.margin=unit(c(1,1,1,1),"lines")) + facet_grid(~Season) + ylab(expression(paste("Organic Carbon [",mu,"g ",m^-3,"]",sep=""))) + xlab("")
p = p + scale_color_brewer(palette="Dark2")

pdf(width=9,height=4,file="images/Raw/FigureS1b_IMRPOVE_OC_Seasons.pdf")
print(p)
dev.off()


# Map
improve.data.annual$Region = "Alaska/Hawaii"
improve.data.annual$Region[improve.data.annual$State %in% westlist] = "West"
improve.data.annual$Region[improve.data.annual$State %in% swlist] = "Southwest"
improve.data.annual$Region[improve.data.annual$State %in% mwlist] = "Central"
improve.data.annual$Region[improve.data.annual$State %in% selist] = "Southeast"
improve.data.annual$Region[!(improve.data.annual$State %in% c(selist,mwlist,swlist,westlist,"AK","HI"))] = "East"

world = ne_countries(scale="medium",returnclass="sf")
usa = world[227,]
stations = st_as_sf(improve.data.annual)

p = ggplot() + geom_sf(data=usa) + geom_sf(data=stations,aes(color=Region),cex=0.5)
p = p + theme(legend.position = "none") + xlim(c(-179,-60)) + ylim(c(15,75)) 
p = p + scale_color_brewer(palette="Dark2") 

pdf(width=3,height=3,file="images/Raw/FigureS1a_IMPROVE_Stations_Map.pdf")
par(mar=c(5,4,1,10))
print(p)
dev.off()


