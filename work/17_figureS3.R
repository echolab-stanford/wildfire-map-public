source("work/00_functions.R")


get_preds = function(x, df, df0) {
    preds = predict(x, df)
    preds_0 = predict(x, df0)
    preds[preds < 0] = 0
    preds_0[preds_0 < 0] = 0
    p = data.frame(preds_0=preds_0, preds=preds)
    p$diff = (p$preds - p$preds_0)
    p$year = df$year
    p$id = df$id
    p$physio_section  = df$physio_section
    p$physio_region = df$physio_region
    p$state =  df$state
    return(p)
}

export_map = function(call_in, type, weights=T) {
    
    p = get_kfold_preds(df, df0, call_in, weights)
    p$diff = p$preds_kfold - p$preds_kfold0
    
    results = merge(data_ll_geo, 
                    p[p$year %in% 2016:2018, ] %>% dplyr::group_by(id) %>% 
                        dplyr::summarize(diff=mean(diff)), by="id")
    
    map = ggplot(results, aes(long,lat,group=group, fill=diff)) + # the data
        geom_polygon() + # make polygons
        scale_fill_gradientn(limits=lim, colors=cols, values=scales::rescale(diff), guide=F) + 
        theme(line = element_blank(),  # remove the background, tickmarks, etc
              axis.text=element_blank(),
              axis.title=element_blank(),
              panel.background = element_blank()) +
        coord_map("bonne", mean(results$lat))
    ggsave(paste0("images/Raw/FigureS3_", type, ".pdf"), map, width=15, height=10, units="cm")
    
    
    print(cor(p$pm, p$preds_kfold, use="c")^2)
    
    p = p %>% dplyr::group_by(year) %>% dplyr::summarize(perc=mean(diff/preds_kfold))
    p$type = type
    
    
    #p = unique(results[, c("id", "diff")]) #use to get full set for jenks breaks
    #p$type=type#use to get full set for jenks breaks
    return(p) 
}

###############################################################
# import data
###############################################################

full_data = readRDS("data/clean/epa_full_long_data_spatial_lag.RDS")
data_ll = readRDS("data/clean/national_grid.RDS")
data_ll_geo = fortify(data_ll, region="id") #this only has the coordinates


###############################################################
# basic data prep, train test split
###############################################################

full_data$year = as.factor(full_data$year)
full_data$state = as.factor(full_data$state)

df = full_data[complete.cases(full_data %>% dplyr::select(-pm, -obs)), ]

df0 = df 
df0[, c(names(df0)[grepl("fire", names(df0)) | grepl("smoke", names(df0))], 
        "f_count",  "mean_f_count")] = 0
df0$area = max(df0$area)

call =  "pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
            sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
            bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1)"

percentages = as.list(rep(NA, 4))
diff =  c(0, 0.28, 0.51, 0.77, 1.11, 1.47, 1.92, 2.52, 5.1) #the set if include log and lin weighted
cols = get_col_regions()[seq(1,100,floor(100/8))]
lim = c(0, ceiling(max(diff)))


###############################################################
# LOG WEIGHTED DISTANCE
###############################################################

p = export_map(paste0(call, "+ smoke_small + smoke_log_med + smoke_log_large"), 
               "log_distance_weighted")

percentages[[1]] = p


###############################################################
# SIMPLE SMOKE COUNT - linear
###############################################################

p = export_map(paste0(call, "+ smoke"), "simple_count")

percentages[[2]] = p


###############################################################
# BINNED COUNT
###############################################################

p = export_map(paste0(call,  "+ smoke_small + smoke_med + smoke_large"), 
               "binned_count")

percentages[[3]] = p


###############################################################
# BINNED COUNT BY FIRE AND SMOKE AND DIST
###############################################################

p = export_map(paste0(call, "+ sfire_sdist + sfire_mdist + sfire_ldist + 
                    mfire_sdist + mfire_mdist + mfire_ldist + 
                    lfire_sdist + lfire_mdist + lfire_ldist"), 
               "fire_binned_count")

percentages[[4]] = p


p = get_kfold_preds(df, df0, paste0(call, "+ sfire_sdist + sfire_mdist + sfire_ldist + 
                    mfire_sdist + mfire_mdist + mfire_ldist + 
                    lfire_sdist + lfire_mdist + lfire_ldist"), T)
p$diff = p$preds_kfold - p$preds_kfold0

results = merge(data_ll_geo, 
                p[p$year %in% 2016:2018, ] %>% dplyr::group_by(id) %>% 
                    dplyr::summarize(diff=mean(diff)), by="id")

pdf(file=paste0("images/Raw/FigureS3_withlegend.pdf"), width=1500, height=1000)
ggplot(results, aes(long,lat,group=group, fill=diff)) + # the data
    geom_polygon() + # make polygons
    scale_fill_gradientn(limits=lim, colors=cols, values=scales::rescale(diff)) + 
    theme(line = element_blank(),  # remove the background, tickmarks, etc
          axis.text=element_blank(),
          axis.title=element_blank(),
          panel.background = element_blank()) +
    coord_map("bonne", mean(results$lat))  
dev.off()


###############################################################
# COMBINING MODELS AND DRAWING LINES
###############################################################

percentages = rbindlist(percentages)
percentages$year = as.numeric(as.character(percentages$year))

png(file=paste0("images/Raw/FigureS3a_overtime.png"), width=20, height=10, res=650, units="in")
ggplot(percentages) + 
    geom_hline(yintercept=c(0, .1, .2, .3, .4, .5), color="grey", lwd=1, alpha=.7) +
    geom_line(aes(year, perc, group=type, color=type), lwd=3.5) + 
     theme_anne(font="sans") + 
    theme(axis.text=element_text(size=45), 
          axis.line = element_line(color="black", size=4)) +
    scale_y_continuous(breaks=c(0, .1, .2, .3, .4), 
                       labels=c("0%", "10%", "20%", "30%", "40%"), limits=c(0, .4)) +
    xlab("") + ylab("") + labs(color="Smoke definition")
dev.off()
