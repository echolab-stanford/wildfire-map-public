source("work/00_functions.R")
library(svMisc)
library(latticeExtra)
library(BAMMtools)
library(ggplot2)
library(caret)

###############################################################
# import data
###############################################################

full_data = readRDS("data/clean/epa_full_long_data_spatial_lag.RDS") 
data_ll = readRDS("data/clean/national_grid.RDS")

###############################################################
# model
###############################################################

full_data$year = as.factor(full_data$year)
full_data$state = as.factor(full_data$state)

full_data$pop_adjall = rowMeans(full_data[, c("pop_adj1", "pop_adj2")])
full_data$mean_pop_adj = rowMeans(full_data[, c("pop", "pop_adj1")])

df = full_data[complete.cases(full_data %>% dplyr::select(-pm, -obs)), ]
df0 = df 
df0 = df0 %>% group_by(id) %>% mutate(weighted_mean_dist = max(weighted_mean_dist))
df0[, unique(c(names(df0)[grepl("smoke", names(df0))], names(df0)[grepl("sfire", names(df0))], 
               names(df0)[grepl("mfire", names(df0))], names(df0)[grepl("lfire", names(df0))], 
       "f_count",  "mean_f_count"))] = 0



#not  inc  aircraft,  gas, industrial process
call = "pm ~ year + physio_section + coaldist2 + coalplant + plant + mining + traffic + 
             taxi + construction + cattle + gas +
             bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) +
             pop_adjall + coalplant_adj1 + mining_adj1 + 
           
             sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
             mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist"


###############################################################
# save predicted data
###############################################################

results_all = get_kfold_preds(df, df0, call, weights=T)
results_all = process_preds(results_all, df0, call)
results_all = results_all[, c("id", "state", "physio_region", "year", "pm", "pop", 
                              "perc_diff", "diff", "preds", "preds0", "preds_kfold", 
                              "preds_kfold0", "smoke", "smoke_day")]
saveRDS(results_all, "data/clean/results_all.RDS")


###############################################################
# plot map results
###############################################################

data_ll_geo = fortify(data_ll, region="id")
cols = get_col_regions()[seq(1,100,length.out=9)]

plot_year = 2016
results = results_all[results_all$year %in% plot_year, ] %>% dplyr::group_by(id) %>%
    dplyr::summarise(diff = mean(diff))
data = merge(data_ll_geo, results, by="id")
at_diff = getJenksBreaks(data$diff, 9)

png(filename=paste0("images/", plot_year, "_smoke_", Sys.Date(), ".png"), 
    width=1500, height=1000)
ggplot(data, aes(long,lat,group=group, fill=diff)) + # the data
    geom_polygon() + # make polygons
    scale_fill_gradientn(limits=c(min(at_diff), max(at_diff)), colors=cols, 
                         values=scales::rescale(at_diff)) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) + 
    coord_map("bonne", mean(data$lat)) + labs(fill="")
dev.off()

