source("work/00_functions.R")
#mortality cost

###############################################################
# read in preds
###############################################################

predictions =  readRDS("data/clean/results_all.RDS")[, c("id", "state", "pop", "year", "pm", "diff", "preds", "smoke_day")]
predictions = predictions %>% 
    dplyr::group_by(id) %>%
    dplyr::mutate(mean_diff = mean(diff), no_smoke = preds-diff, pop=mean(pop))

census =  read.csv("data/census/ACS_17_1YR_S0101.csv", skip=1, stringsAsFactors=F)
census = census[, c("Id2",  "Total..Estimate..Total.population", 
                    names(census)[grepl("Total..Estimate..AGE...", names(census))])]
census$population = census$Total..Estimate..Total.population
census$over65 = rowSums(census[, 16:20])/census$population
census$a65_70 = census$Total..Estimate..AGE...65.to.69.years/census$population
census$a70_75 = census$Total..Estimate..AGE...70.to.74.years/census$population
census$a75_80 = census$Total..Estimate..AGE...75.to.79.years/census$population
census$a80_85 = census$Total..Estimate..AGE...80.to.84.years/census$population
census$a85 = census$Total..Estimate..AGE...85.years.and.over/census$population
census = census[, c(1, 22:27)]
census$state = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "DC", "FL", "GA", "HI", 
                 "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", 
                 "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH",
                 "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT",  "VT", "VA", "WA", 
                 "WV", "WI", "WY", "PR")

#get the % in each group and multiply by the cell level total pop to get age group pop
predictions = merge(predictions, census, by="state")
predictions[, 12:17] = predictions[, 12:17] * predictions$pop
ages = predictions[, 12:17]
predictions$over65 = predictions$over65/100000

#get the total baseline mortality, weighted by distribution of age groups
age_weighting = c(sum(census$a65_70)+sum(census$a70_75), 
                  sum(census$a75_80)+sum(census$a80_85), sum(census$a85))
all_bmr = wtd.mean(c(1790.9, 4472.6, 13573.6), weights=age_weighting)

#in usa ncd = 89% of all deaths (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6211719/)
#   all_cause*.89 + lower respiratory infections
#cdc National Vital Statistics Reports Volume 64, Number 2
# adding lower respiratory infection deaths to ncd deaths
ncd_bmr = all_bmr*0.89 + wtd.mean(c(29.5+141.2, 103.7+367, 441+699.3), age_weighting)

###############################################################
# get the cost functions
###############################################################

# Increases of 10 μg per cubic meter in PM2.5 and of 10 ppb in ozone were associated with increases in all-cause mortality of 7.3% 
# DI 2017
rate = 0.0073
di = function(X, r=rate) {Y=1+(r*(X-5)); Y[X<5]=1; return((Y-1)*all_bmr)}  #APPROXIMATES LINE FROM DI

#burnett mortality for ncd+lri
theta = 0.1430
burnett_form = function(X, theta, alpha, mu, v) {
    X = sapply(X, function(x){max(0, x-2.4)}) 
    one = log(1+ (X/alpha))
    two = 1/(1 + exp(-(X-mu)/v))
    Y = exp(theta*one*two)
    return(Y)
}
burnett = function(X, t=theta) {
    
    ncd_lri = burnett_form(X, t, 1.6, 15.5, 36.8)
    ncd_lri = (ncd_lri-1)*ncd_bmr
    return(ncd_lri)
}

## get the epa beta
## using the different parametric distributions in the EPA documentation
set.seed(5)
expa = rtruncnorm(1000, a=0, mean=1.42, sd=0.89)
expc = rtruncnorm(1000, a=0, mean=1.2, sd=0.49)
expd = rtriangle(1000, 0.1, 1.6, 0.95)
expe = rtruncnorm(1000, a=0, mean=2, sd=0.61)
expg = rtruncnorm(1000, a=0, mean=1, sd=0.19)
expi = rtruncnorm(1000, a=0, b=2.273, mean=1.25, sd=0.53)
expj = rweibull(1000, 2.21, 1.41)
epa = c(expa, expc, expd, expe, expg, expi, expj)
beta = mean(epa/100)

epa = function(X, b=beta) {(exp(b*X)-1)*all_bmr}


###############################################################
# Plotting
###############################################################

min = 0
max = 30
X = seq(min, max, 0.5)
rates = data.frame(x=X, di=di(X),  burnett=burnett(X), epa=epa(X))


#####  PLOT THE FUNCTIONS

pdf(file="images/Figure3-c.pdf",width=4,height=4,useDingbats = F)
par(mar=c(3.5, 4, 0, 0))

plot(rates[,c("x","di")],type="l",col="cyan",ylab="",axes=F,xlab="", 
     ylim=c(min(rates[, -1]), max(rates[, 2:4], na.rm=T)), lwd=1.5)

w=which(rates$x==max)-2
lines(rates$x,rates$di,col="cyan", lwd=1.5)
text(max*.98,rates[w,]$di, "Di", cex=1.2, col="cyan", pos=2)

lines(rates$x,rates$epa,col="purple", lwd=1.5)
text(max*.98,rates[w,]$epa, "EPA", cex=1.2, col="purple", pos=2)

lines(rates$x,rates$burnett,col="orange", lwd=1.5)
text(max*.98,rates[w,]$burnett, "Burnett", cex=1.2, col="orange", pos=2)

axis(2,las=1)
axis(1)
#title(ylab="Mortality from PM2.5 per 100k per year")
mtext("Mortality from PM2.5 per 100k per year", side=2, line=3, adj=-0.5)

dev.off()

###############################################################
#mortality calculations
###############################################################

get_mortality = function(c) {Y = c(sum(di(c)*predictions$over65),
                                   sum(burnett(c)*predictions$over65),
                                   sum(epa(c)*predictions$over65))
                             return(Y)}

c = predictions$preds - (predictions$smoke_day*2.259)/365
miller = get_mortality(c)
miller_data = c

c = predictions$preds - predictions$diff + predictions$mean_diff*1.5
increased = get_mortality(c)
increased_data = c

c = predictions$preds - predictions$diff + predictions$mean_diff*1.25
increased_slightly = get_mortality(c)
increased_slightly_data = c


c = predictions$no_smoke
no_smoke = get_mortality(c)
no_smoke_data = c

#COCHRANE 2012 treatments prevented 4.8ha of burning for every hectare of promoted burning
c = predictions$preds - predictions$diff + predictions$mean_diff*.5
reduced_from_burning = get_mortality(c)
reduced_from_burning_data = c

#COCHRANE 2012 treatments prevented 4.8ha of burning for every hectare of promoted burning
c = predictions$preds - predictions$diff + predictions$mean_diff*.85
somereduced_from_burning = get_mortality(c)
somereduced_from_burning_data = c


#just averaging all smoke pm over relevant years
c = predictions$preds - predictions$diff + predictions$mean_diff
equal_with_burning = get_mortality(c)
equal_with_burning_data = c

#observed
c = predictions$preds
observed = get_mortality(c)

costs =  data.frame(cost=c("di", "burnett", "epa"), 
                    obs=observed, burn=equal_with_burning, 
                    burn_some_better=somereduced_from_burning, 
                    burn_better=reduced_from_burning, no_smoke=no_smoke,
                    increased_slightly=increased_slightly, increased=increased)#, miller=miller)
costs =  as.data.frame(t(costs))
names(costs) = as.character(unlist(costs[1,]))
costs = costs[-1,]
for (i  in 1:ncol(costs))  {costs[,i] = as.numeric(as.character(costs[,i]))}

norm = function(X) {X-X[1]}
costs = apply(costs, 2, norm)
costs = as.data.frame(costs)
costs = costs/13


###############################################################
# Plotting part II
###############################################################

bw = 0.5
sc = 0.8

d  = data.frame(
           preds=predictions$preds, 
           burn=equal_with_burning_data,
           burn1.8=somereduced_from_burning_data,
           burn3.9=reduced_from_burning_data,  
           no_smoke=predictions$no_smoke,
           increased_slightly=increased_slightly_data, 
           increased=increased_data)
breaks = seq(min, max, bw)
mx = max(sapply(d, function(x) max(hist(x,breaks=breaks,plot=F)$counts)))
clz = apply(sapply(c("grey"), col2rgb)/255, 2, 
            function(x) rgb(x[1], x[2], x[3], alpha=0.85))


pdf(file="images/Figure3-ab.pdf",width=8,height=6,useDingbats = F)
layout(matrix(c(1,2), ncol=2), width=c(1, 2))
par(mar=c(4, 0, 0, 0))

plot(1,type="n",xlim=c(0,18),ylim=c(1,length(d)+1),axes=F,ylab="",xlab="PM2.5")
axis(1)
#segments(x0=12, y0=1, y1=length(d)+1,  col="#0000FF90", lwd=0.5)


for (i in 1:length(d)) {
    
    j = length(d) - i + 1 
    
    # plot background baseline dist
    ch = hist(d[[1]],breaks=breaks,plot=F)
    ch$counts = (ch$counts/mx*sc)[1:which(breaks==15)] #scale to max 1 by dividing by baseline max, then scale again
    ch$breaks = ch$breaks[1:which(breaks==15)]
    lines(c(rep(ch$breaks, each=2)[-1], 15), rep(ch$counts, each=2)+j, 
          lwd=1.5, col="black")
    
    hh = hist(d[[i]],breaks=breaks,plot=F)
    hh$counts = hh$counts/mx*sc
    rect(breaks,rep(j,length(breaks)),breaks+bw,hh$counts+j,col=clz,border=NA)
    
    if (i==1) {
        lines(c(rep(ch$breaks, each=2)[-1], 15), rep(ch$counts, each=2)+j, 
                         lwd=1.5, col="black")
    }  else {
        lines(c(rep(hh$breaks[1:30], each=2)[-1], 15), rep(hh$counts[1:30], each=2)+j, 
              lwd=1.5, col="grey50")
    }
}


#####  MAKE COLS

n = 5000
xvals = rev(1:nrow(costs))
at = seq(n*round(min(costs)/n),  n*round(max(costs)/n), n)
plot(1,1, col=NA, xlim=c(min(costs), max(costs)), ylim=c(1,length(d)+1), 
     axes=F, xlab="", ylab="")
segments(y0=1, y1=length(d)+1, x0=at, col="#7F7F7F80", lwd=0.5)

#make bars
rect(xleft=0, xright=costs$di, ybottom=xvals+0.5, ytop=xvals+0.7, col="cyan", lwd=0)
rect(xleft=0, xright=costs$burnett, ybottom=xvals+0.3, ytop=xvals+0.5, col="orange", lwd=0)
rect(xleft=0, xright=(costs$epa), ybottom=xvals+0.1, ytop=xvals+0.3, col="purple", lwd=0)

segments(y0=0.5, y1=8.5, x0=0,  col='black', lwd=0.5)

axis(1, at=at)
mtext("Additional premature deaths/year of adults 65+, rel to baseline", side=1, line=2.5)

dev.off()
