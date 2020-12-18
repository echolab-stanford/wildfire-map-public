results_all = readRDS("data/clean/results_all.RDS")

if (!"VanD_comparisons.Rds" %in% list.files("data/clean"))  {
    out = list() #this takes a wile
    for  (i in 1:length(years)) {
        y = years[i]
        r = raster(paste0("data/vanD/vanD-", y, ".asc"))
        r = velox(r)
        x = r$extract(data_ll, fun=mean_na)
        x = data.frame(id = rownames(x), vanD_pm = x, year = y)
        out[[i]] = x
    }
    
    x = rbindlist(out)
    saveRDS(x[, c("id", "year", "vanD_pm")], "data/clean/VanD_comparisons.Rds")
} else {
    x = readRDS("data/clean/VanD_comparisons.Rds")
}

x = merge(x, results_all, by=c("id", "year"))
names(x)[c(6, 10, 3)] = c("Measured PM2.5",  "Predicted PM2.5", "Predicted PM2.5 - VanD")

pdf('images/FigureS4_pmtomodel.pdf',width=5,height=5,useDingbats=F)
panel(x[, c("Predicted PM2.5", "Measured PM2.5")], square=T, font="sans", size=15) + 
    geom_abline(intercept=0, slope=1, lwd=0.5, lty=2, color="lightgrey")
dev.off()

pdf('images/FigureS4_pmtovanD.pdf',width=5,height=5,useDingbats=F)
panel(x[, c("Predicted PM2.5 - VanD", "Measured PM2.5")], square=T, font="sans", size=15) + 
    geom_abline(intercept=0, slope=1, lwd=0.5, lty=2, color="lightgrey")
dev.off()



### HISTOGRAMS
cross = data.frame(physio_region = unique(x$physio_region), 
                   region_name = c("Pacific Coast", "Rockies", "Plains", "Plateaus", 
                                   "Appalachia", "East Coast"))
x = readRDS("data/clean/VanD_comparisons.Rds")
x = merge(x, results_all, by=c("id", "year"))
x = merge(x, cross)
y = x %>% dplyr::summarise(r2_model = cor(preds, pm, use="c")^2, 
                    r2_vanD = cor(vanD_pm, pm, use="c")^2) %>% mutate(region_name = "Lower 48")
z =  x %>% group_by(region_name) %>%
    dplyr::summarise(r2_model = cor(preds, pm, use="c")^2, r2_vanD = cor(vanD_pm, pm, use="c")^2)
z = rbind(z, y)
z = melt(z, id="region_name")

z$variable = factor(z$variable, levels=c("r2_model", "r2_vanD"))
z$region_name = factor(z$region_name, levels=c("Plateaus", "Rockies", "Pacific Coast", 
                                               "East Coast", "Lower 48", "Appalachia", 
                                               "Plains"))
pdf('images/FigureS4_regionR2.pdf',width=6,height=5,useDingbats=F)
ggplot(z) + geom_col(aes(x=region_name, y=value, fill=variable), position="dodge", width=0.7) + 
    theme_anne(size=20, font="sans") + xlab("") + ylab("R2") + guides(fill=F) + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.line.x = element_line(color = "black", size = 0.2), 
          axis.line.y = element_line(color = "black", size = 0.2))
dev.off()


z =  x %>% group_by(year) %>%
    dplyr::summarise(r2_model = cor(preds, pm, use="c")^2, r2_vanD = cor(vanD_pm, pm, use="c")^2)
z = melt(z, id="year")
z$variable = factor(z$variable, levels=c("r2_model", "r2_vanD"))

pdf('images/FigureS4_yearR2.pdf',width=6,height=4.5,useDingbats=F)
ggplot(z) + geom_col(aes(x=year, y=value, fill=variable), position="dodge", width=0.7) + 
    theme_anne(size=20, font="sans") + xlab("") + ylab("R2") + guides(fill=F) + 
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.line.x = element_line(color = "black", size = 0.2), 
          axis.line.y = element_line(color = "black", size = 0.2))
dev.off()
