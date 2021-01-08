source("work/00_functions.R")

data = readRDS("data/clean/epa_full_long_data_spatial_lag.RDS")
data = unique(data[, c("id", "physio_region", "physio_section")])
data_ll = readRDS("data/clean/national_grid.RDS")
data_ll_geo = fortify(data_ll, region="id") #this only has the coordinates
data_ll_geo = merge(data_ll_geo, data, by="id")

x = unionSpatialPolygons(data_ll, data_ll$physio_section)
x_geo = fortify(x, region="id")

data_ll_geo$physio_section = factor(data_ll_geo$physio_section, 
            c( "PUGET TROUGH", "WEST GULF COASTAL PLAIN", "SONORAN DESERT", "HIGH PLAINS", "SOUTHERN ROCKY MOUNTAINS",
              "SUPERIOR UPLAND", "CALIFORNIA TROUGH", "NEW ENGLAND", "MEXICAN HIGHLAND", "FLORIDIAN", 
              "NORTHERN ROCKY MOUNTAINS", "PUGET SOUND", "GREAT BASIN", "PIEDMONT", "PLAINS BORDER", 
              "EASTERN LAKE", "COLORADO PLATEAUS", "COLUMBIA PLATEAU", "CALIFORNIA COAST RANGES", 
              "DISSECTED TILL PLAINS", "EAST GULF COASTAL PLAIN", "MISSOURI PLATEAU, UNGLACIATED",
              "LOWER CALIFORNIAN", "CASCADE-SIERRA MOUNTAINS", "ADIRONDACK", "INTERIOR LOW PLATEAUS", 
              "OREGON COAST RANGE", "VALLEY AND RIDGE", "COLORADO PIEDMONT", "TILL PLAINS",
              "MIDDLE ROCKY MOUNTAINS", "OZARK PLATEAUS", "ADIRONDACKS", "EMBAYED", 
              "KLAMATH MOUNTAINS", "WYOMING BASIN", "WESTERN LAKE",  "MISSOURI PLATEAU, GLACIATED",
              "MISSISSIPPI ALLUVIAL PLAIN", "OSAGE PLAINS", "SEA ISLAND", 
              "WISCONSIN DRIFTLESS", "PECOS VALLEY", "APPALACHIAN PLATEAUS"))
map = ggplot() + # the data
    geom_polygon(aes(long, lat, group=group, fill=physio_section), data=data_ll_geo) + # make polygons 
    geom_polygon(aes(long, lat, group=group), data=x_geo, color="black", fill=NA, lwd=0.5) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) +
    coord_map("bonne", mean(data_ll_geo$lat)) + guides(fill=F)

ggsave("images/Raw/FigureS2.pdf", map, width=15, height=10, units="cm")
