source("work/00_functions.R")

data_ll = readRDS("data/clean/national_grid.RDS")

#load in fire perimeter data
prescribed = read.csv("data/fire/prescribed_burn_acres.csv") %>% 
    dplyr::filter(region != "AK") 
prescribed = prescribed %>% group_by(year, region) %>%
    dplyr::summarise(acres = sum(acres), 
                     total_acres = sum_na(as.numeric(as.character(total_acres)), na.rm=T))
prescribed$hectares = ((prescribed$total_acres)*0.404686)/100000 #to mil hectares
prescribed = prescribed[, c("year", "region", "hectares")]

#regions of the fire data
gacc = readOGR("data/boundaries/GACC",  "National_GACC_Current_20200226")
g = over(gacc, data_ll, returnList = T)
gid = gacc$GACCAbbrev
for (i in 1:length(g)) {if (nrow(g[[i]])>0) {g[[i]]$gacc = gid[i]}}
g =  rbindlist(g[2:length(g)])[, c("id", "gacc")]
g$gacc = substr(g$gacc, 1, 2)
g$gacc[g$gacc == "ON"] =  "NO"
g$gacc[g$gacc == "OS"] =  "SO"

#results
results = readRDS("data/clean/results_all.RDS")[, c("id", "year", "diff")]
results = merge(results, g, by="id")
results =  merge(results, prescribed, by.x=c("year", "gacc"), by.y=c("year", "region"))

results = results %>% group_by(gacc, year) %>% 
    dplyr::summarise(diff=mean(diff), hectares=mean(hectares))

results$xhat = summary(lm(log(hectares) ~ gacc*as.numeric(year), 
                          data=results))$residuals
results$yhat = summary(lm(log(diff) ~ gacc*as.numeric(year), 
                          data=results))$residuals

names(results)[5:6] = c("hectares burned (log deviation from trend)", 
                        "PM from smoke (log deviation from trend)")
gg = panel(results[, 5:6], font="sans")
ggsave("images/Final/FigureS6.pdf", gg)