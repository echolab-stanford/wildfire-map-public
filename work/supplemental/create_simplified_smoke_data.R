rm(list = ls())
library(raster)
library(tidyverse)
library(sf)
library(rgdal)
library(data.table)

# Written by Sam Heft-Neal
###################################           
#write or load plume matches     
###################################  


#load all smoke plume shapefiles into list named by date
files <- list.files("data/smoke/")
files <- files[grep(x = files, pattern = ".shp")]

#loop over raw daily smoke plume data files and read into list

smoke <- list(); length(smoke) <- length(files); names(smoke)<-substr(files, 10,17);
for(i in 1:length(files)){
    try(
        smoke[[i]] <- st_read(paste("data/smoke/",files[i], sep = ""))
    )
    try(
        smoke[[i]]$date <- as.numeric(substr(files[i],10,17))
    )
}#end i loop


#fill out smoke plume data but use try() because variables don't exist for some days
for(i in 1:length(smoke)){
    try(smoke[[i]]$ID2 <- 1:nrow(smoke[[i]]))
    if("Density" %in% names(smoke[[i]]) ==F){
        try( smoke[[i]]$Density <- NA)
    }
}


write_rds("data/clean/smoke_plumes.rds")

