source("work/00_functions.R")

################################################################################
################################################################################
# Numbers for table of R^2 by region and year
################################################################################
################################################################################

vanD = readRDS("data/clean/VanD_comparisons.Rds")

region = vanD %>% group_by(physio_region) %>% 
    summarise(r2_vanD=cor(vanD_pm, pm, use="c")^2, r2_model=cor(preds, pm, use="c")^2)
region$physio_region = c("Appalachia", "East Coast", "Plains", "Plateaus", "Pacific Coast", "Rockies")
region = pivot_longer(region, c("r2_vanD", "r2_model"))

years = vanD %>% group_by(year) %>% 
    summarise(r2_vanD=cor(vanD_pm, pm, use="c")^2, r2_model=cor(preds, pm, use="c")^2)
years = pivot_longer(years, c("r2_vanD", "r2_model"))

overall = vanD %>% 
    summarise(r2_vanD=cor(vanD_pm, pm, use="c")^2, r2_model=cor(preds, pm, use="c")^2) %>%
    mutate(physio_region="Lower  48")
overall = pivot_longer(overall, c("r2_vanD", "r2_model"))
region = rbind(region,  overall)

r = ggplot(region, aes(physio_region, value, fill=name))  + 
    geom_bar(stat="identity", width=.5, position="dodge")  + 
    #theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(fill="") + 
    ylim(0, .8) +  xlab("Regions") + ylab("R^2") + guides(fill=F)

y = ggplot(years, aes(year, value, fill=name))  + 
    geom_bar(stat="identity", width=.5, position="dodge")  + 
    #theme(axis.text.x = element_text(angle=45, hjust=1)) + labs(fill="") + 
    ylim(0, .8) + xlab("Year") + ylab("R^2")

x = ggarrange(r, y, labels = c("a", "b"), ncol = 2, nrow = 1) 


################################################################################
################################################################################
# Numbers for table on removing co-variates
################################################################################
################################################################################


get_preds = function(x, df, df0) {
    preds = predict(x, df)
    preds_0 = predict(x, df0)
    preds[preds < 0] = 0
    preds_0[preds_0 < 0] = 0
    p = data.frame(preds_0=preds_0, preds=preds)
    diff = (p$preds - p$preds_0)/p$preds
    return(diff)
}

get_preds_pm = function(x, df, df0) {
    preds = predict(x, df)
    preds_0 = predict(x, df0)
    preds[preds < 0] = 0
    preds_0[preds_0 < 0] = 0
    p = data.frame(preds_0=preds_0, preds=preds)
    diff = (p$preds - df$pm)/df$pm
    return(diff)
}

###############################################################
# read in data
###############################################################

full_data = readRDS("data/clean/epa_full_long_data_spatial_lag.RDS") 
full_data$year = as.factor(full_data$year)
full_data$state = as.factor(full_data$state)

df = full_data[complete.cases(full_data %>% dplyr::select(-pm, -obs)), ]
df0 = df 
df0[, c(names(df0)[grepl("fire", names(df0)) | grepl("smoke", names(df0))], 
        "f_count",  "mean_f_count")] = 0
df0$area = max(df0$area)

################################################################################
################################################################################
# get numbers for table on sensitivity to covariates
################################################################################
################################################################################

m = lm(pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
           sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
baseline = mean(get_preds(m, df, df0))

m = lm(pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
           sqrt(mining) + sqrt(tourism) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
traffic = mean(get_preds(m, df, df0))

m = lm(pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
           sqrt(tourism) + sqrt(traffic) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
mining = mean(get_preds(m, df, df0))

m = lm(pm ~ year + physio_section +  sqrt(cattle) + 
           sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
coal = mean(get_preds(m, df, df0))

m = lm(pm ~ year + physio_section + coaldist2 + coalplant + plant + 
           sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
cattle = mean(get_preds(m, df, df0))

m = lm(pm ~ year + physio_section +
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
all = mean(get_preds(m, df, df0))

m = lm(pm ~ sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)
all_fe = mean(get_preds(m, df, df0))


baseline
coal
cattle
mining
traffic
all
all_fe

################################################################################
################################################################################
# get numbers for table on sensitivity to divisor
################################################################################
################################################################################

df = df[!is.na(df$pm), ]
df0 = df0[!is.na(df0$pm), ]

m = lm(pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
           sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
           bs(pop, knots = quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1) + 
           sfire_sdist + sfire_mdist + sfire_ldist + mfire_sdist + mfire_mdist + 
           mfire_ldist + lfire_sdist + lfire_mdist + lfire_ldist, 
       df)

mean(get_preds(m, df, df0))
mean(get_preds_pm(m, df, df0))


################################################################################
################################################################################
# Numbers for table percent improvement of r2 from smoke vars
################################################################################
################################################################################

full_data = readRDS("data/clean/epa_full_long_data_spatial_lag.RDS") 

full_data$year = as.factor(full_data$year)
full_data$state = as.factor(full_data$state)

df = full_data[complete.cases(full_data %>% dplyr::select(-pm, -obs)), ]
df0 = df 
df0[, c(names(df0)[grepl("fire", names(df0)) | grepl("smoke", names(df0))], 
        "f_count",  "mean_f_count")] = 0
df0$area = max(df0$area)

call =  "pm ~ year + physio_section + coaldist2 + coalplant + plant + sqrt(cattle) + 
sqrt(mining) + sqrt(tourism) + sqrt(traffic) + 
bs(pop, knots =  quantile(pop, c(0.6, 0.9, 0.95, 0.99)), degree = 1)"


base = get_kfold_preds(df, df0, call, weights=T)
log = get_kfold_preds(df, df0, paste0(call, "+ smoke_small + smoke_log_med + smoke_log_large")) 
count = get_kfold_preds(df, df0, paste0(call, "+ smoke")) 
binned = get_kfold_preds(df, df0, paste0(call, "+ smoke_small + smoke_med + smoke_large")) 
binned_fire = get_kfold_preds(df, df0, paste0(call, "+ sfire_sdist + sfire_mdist + sfire_ldist + 
                        mfire_sdist + mfire_mdist + mfire_ldist + 
                                      lfire_sdist + lfire_mdist + lfire_ldist")) 


base = cor(base$pm, base$preds_kfold, use="c")^2
log = cor(log$pm, log$preds_kfold, use="c")^2
count = cor(count$pm, count$preds_kfold, use="c")^2
binned = cor(binned$pm, binned$preds_kfold, use="c")^2
binned_fire = cor(binned_fire$pm, binned_fire$preds_kfold, use="c")^2

log = (log - base)/base
count = (count - base)/base
binned = (binned - base)/base
binned_fire = (binned_fire - base)/base

d = data.frame(model = c("Plumes binned by fire size and distance", 
                         "Sum of plumes, weighted by log distance to fire", 
                         "Plumes binned by size of plume", "Sum of plumes"), 
               improve = c(binned_fire, log, binned, count)*100)
d$model = factor(d$model, levels=unique(as.character(d$model)) )

print(d)
