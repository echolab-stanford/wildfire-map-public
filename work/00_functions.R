library(BAMMtools)
library(caret)
library(cleangeo)
library(data.table)
library(dplyr)
library(gdata) #to read xls
library(geosphere)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(Hmisc)
library(hutils)
library(imputeTS) #for na_kalman
library(latticeExtra)
library(mapproj)
library(maptools)
library(ncdf4)
library(openxlsx) # to read xlsx
library(plyr)
library(raster)
library(RColorBrewer)
library(readr)
library(rgeos)
library(rgdal)
library(RSelenium)
library(sf)
library(signal) #for interpolation
library(sp)
library(splines)
library(stringr)
library(svMisc)
library(tidyr)
library(triangle)
library(truncnorm)
library(velox)
library(zoo)




##################################################################
# Settings
##################################################################

buffer = 10 #km
crs_using = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #"+init=epsg:4238" #
km = 100 #raster cell size in km across
years = 2006:2018
plant_radius = c(-20, -10, 20, 40)
max_dist = 10


##################################################################
# Functions
##################################################################

rmse = function(obs, pred) {
    sqrt(mean((pred-obs)^2))
}
sum_na = function(x, ...) if (all(is.na(x))) NaN else sum(x, na.rm = TRUE)
max_na = function(x, ...) if (all(is.na(x))) NaN else max(x, na.rm = TRUE)
min_na = function(x, ...) if (all(is.na(x))) NaN else min(x, na.rm = TRUE)

# function that takes input data with inconsistent years available and interpolates 
#   between them using a linear model
lm_interpolate = function(data, column, id="id", interp_years=2006:2018) {
    if (any(names(data) == "id") & id != "id") {
        names(data)[names(data) == "id"] = "old_id"
    }
    names(data)[names(data) == id] = "id"
    names(data)[names(data) == column] = "target"
    data$year = as.numeric(data$year)
    
    pad = expand.grid(id=unique(data$id), year=interp_years)
    pad = merge(pad, unique(data[, c("id", "year")]), by=c("id", "year"), all=T)
    data = merge(data, pad, by=c("id", "year"), all=T)
    data = data[order(data$year),]
    
    for (i in unique(data$id)) {
        cur = data[data$id == i, ]
        years = sort(cur[is.na(cur$target), ]$year)
        if (sum(!is.na(cur$target)) == 1) {
            pred = max(cur$target, na.rm=T)
        } else if (sum(!is.na(cur$target)) == 0) {
            next
        } else {
            m = lm(target ~ year, data=cur)
            pred = predict(m, data.frame(year=years))
            pred[pred < 0] = min(cur$target, na.rm=T)
        }
        data[data$id == i & data$year %in% years, "target"] = pred
    }
    
    names(data)[names(data) == "target"] = column
    names(data)[names(data) == "id"] = id
    
    return(data)
}

# function that takes input data with inconsistent years available and interpolates 
#   linearly between them
interpolate = function(data, column, id="id", interp_years=2006:2018) {
  
  if (any(names(data) == "id") & id != "id") {
    names(data)[names(data) == "id"] = "old_id"
  }
  names(data)[names(data) == id] = "id"
  names(data)[names(data) == column] = "target"
  
  pad = expand.grid(id=unique(data$id), year=interp_years, target=NA)
  data = data[order(data$year),]

  #loop over years you don't have
  for (i in unique(data$id)) {
      cur = data[data$id == i, ]
      cur = cur[complete.cases(cur), ]
      if (nrow(cur) == 1) {
          interp = cur$target
      }  else if (nrow(cur) == 0) {interp=NA} else {
          interp = interp1(cur$year, cur$target, interp_years, method="linear", extrap=T)
      }
      pad[pad$id == i, ]$target = interp
  }
  
  #combine all data
  names(pad) = c(id, "year", column)
  
  return(pad)
}

mean_interpolate = function(data, column, id="id", interp_years=2006:2018) {
    if (any(names(data) == "id") & id != "id") {
        names(data)[names(data) == "id"] = "old_id"
    }
    names(data)[names(data) == id] = "id"
    names(data)[names(data) == column] = "target"
    
    means = data %>% dplyr::group_by(id) %>% 
        dplyr::summarize(mean = mean(target, na.rm=T))
    
    pad = expand.grid(id=unique(data$id), year=interp_years)
    data = merge(data, pad, by=c("id", "year"), all=T)
    empty = data[is.na(data$target), ]
    empty = merge(empty, means, by="id", all.x=T)
    data[is.na(data$target), ] = empty[, c("id", "year", "mean")]
    
    names(data) = c("id", "year", column)
    
    return(data)
}

kalman_mean_interpolate =  function(data, column, id="id", interp_years=2006:2018) {
    if (any(names(data) == "id") & id != "id") {
        names(data)[names(data) == "id"] = "old_id"
    }
    names(data)[names(data) == id] = "id"
    names(data)[names(data) == column] = "target"
    data$year = as.numeric(data$year)
    
    pad = expand.grid(id=unique(data$id), year=interp_years)
    pad = merge(pad, unique(data[, c("id", "year")]), by=c("id", "year"), all=T)
    data = merge(data, pad, by=c("id", "year"), all=T)
    data = data[order(data$year),]
    
    #get ids for  kalman interp
    ids_kalman = data %>% group_by(id) %>% 
        dplyr::summarize(n=sum(!is.na(target)), u=length(unique(target[!is.na(target)]))) %>% 
        dplyr::filter(n>=3 & u>1)
    ids_kalman = ids_kalman$id
    
    #fill in the ones that can with kalman
    data[data$id %in% ids_kalman, ] = data %>% 
        dplyr::filter(id %in% ids_kalman) %>% 
        group_by(id) %>% 
        mutate(target = na_kalman(target, type="level"))
    data$target[data$target<0] = 0
    
    #get ids for mean
    ids_mean = data %>% group_by(id) %>% 
        dplyr::summarize(n=sum(!is.na(target)), u=length(unique(target[!is.na(target)]))) %>% 
        dplyr::filter(n>0 & (n<3 | u<=1))
    ids_mean = ids_mean$id
    
    #fill in the ones with  mean
    data[data$id %in% ids_mean, ] = mean_interpolate(data[data$id %in% ids_mean, ], "target")
    
    names(data)[names(data) == "target"] = column
    names(data)[names(data) == "id"] = id
    
    return(data)
}

get_fire_count = function(fire, dates_fire=names(fire), geo, geo_id="id") {
  
  #create fire_dists to store the distances for each day in
  fire_dists = list()
  #rename to use id
  if (any(names(geo) == "id") & geo_id != "id") {
    names(geo)[names(geo) == "id"] = "old_id"
  }
  names(geo)[names(geo) == geo_id] = "id"
  
  years = unique(substr(dates_fire, 1, 4))
  months = unique(substr(dates_fire, 5, 6))
  length(fire_dists) = length(length(years)*length(months))
  
  
  #loop through days and find distances for each
  i=1
  prog = txtProgressBar(min=0, max=length(years)*length(months), initial=0, 
                        char="-", style=3)
  for (year in years) {
      for (month in months) {
          
          cur = paste0(year, month)
          cur = which(startsWith(dates_fire, cur))
          cur = lapply(fire[cur], function(x) {if (!is.null(x)) {as.data.table(st_coordinates(x))}} )
          cur = lapply(cur, function(x) { if (!is.null(x)) {x[complete.cases(x), ]}  })
          cur = SpatialPoints(rbindlist(cur[sapply(cur, nrow) > 0]))
          proj4string(cur) = crs(geo)
          
          #find the nearest fire for each id
          num = over(geo, cur, returnList=T)
          num = sapply(num, length)
          
          fire_dists[[i]] = data.frame(id=geo$id, year=year, month=month, fires=num)
          names(fire_dists)[[i]] = paste0(year, month)
          
          setTxtProgressBar(prog, i)
          i = i + 1
      }
  }
  
  fire_dists = rbindlist(fire_dists)

  return(fire_dists)
}

get_smoke_plumes = function(smoke, geo, geo_id) {
  
  start = Sys.time()
  
  #rename to use id
  if (any(names(geo) == "id") & geo_id != "id") {
    names(geo)[names(geo) == "id"] = "old_id"
  }
  names(geo)[names(geo) == geo_id] = "id"
  
  #create a dataframe with a row for each mzip-day
  overlaps = as.list(rep(NA, length(names(smoke))))
  geo_df = data.frame(id = as.character(unique(geo$id)), den0=0)
  
  if (length(unique(geo$id)) != nrow(geo)) {warning("ID is not unique")}
  
  #loop through all the smoke files, where name(smoke) is the day
  prog = txtProgressBar(min=0, max=length(smoke), initial=0, char="-", style=3)
  for (i in 1:length(names(smoke))){
    
    day = names(smoke)[i] #the string representing the day
    s = smoke[[day]] #the shapefile
    temp_df = geo_df #df with row for each mzip
    
    if (nrow(s) == 0) {
      den = 0 
    } else {
      
      s = tryCatch({
        as_Spatial(s)
      }, error = function(e) {
        s = s[!is.na(st_is_valid(s)),]
        s = st_cast(s, "POLYGON")
        as_Spatial(s)
      })
      crs(s) = crs(geo)
      den = colSums(gIntersects(geo, s, byid = TRUE))
    }
    
    temp_df[, "den0"] = den
    
    temp_df$date = day
    overlaps[[i]] = temp_df
    
    if (i%%10 == 0) {
      setTxtProgressBar(prog, i)
    }
    
  }
  
  #combine all the days into one df and munge
  data = rbindlist(overlaps)
  data$smoke_day = ifelse(data$den0 > 0, 1, 0)
  data$year = as.character(substr(data$date, 1, 4))
  data = data %>% dplyr::group_by(id, year) %>% 
    dplyr::summarise(
      smoke_dayv = var(smoke_day, na.rm = T),
      smoke_day = sum(smoke_day, na.rm = T),
      den0v = var(den0, na.rm = T),
      den0 = sum(den0, na.rm = T)
    ) 
  
  end = Sys.time() - start
  print(end)
  
  #return df with smoke vals for all the geos provided 
  # for subset of smoke provided
  return(data)
}

# SOURCED FROM https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
    
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(old_path, "/Users/annedriscoll/opt/miniconda3/lib/python3.7/site-packages", sep = ":"))
    
    if (isTRUE(readpoly)) require(rgdal)
    if (is.null(pypath)) {
        pypath <- Sys.which('gdal_polygonize.py')
    }
    if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(dirname(pypath))
    if (!is.null(outshape)) {
        outshape <- sub('\\.shp$', '', outshape)
        f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
        if (any(f.exists))
            stop(sprintf('File already exists: %s',
                         toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                        sep='.')[f.exists])), call.=FALSE)
    } else outshape <- tempfile()
    if (is(x, 'Raster')) {
        require(raster)
        writeRaster(x, {f <- tempfile(fileext='.tif')})
        rastpath <- normalizePath(f)
    } else if (is.character(x)) {
        rastpath <- normalizePath(x)
    } else stop('x must be a file path (character string), or a Raster object.')
    system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                    pypath, rastpath, gdalformat, outshape)))
    if (isTRUE(readpoly)) {
        shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
        return(shp)
    }
    return(NULL)
}


get_kfold_preds = function(df, df0, call, weights=T, k=10) {
    
    set.seed(198432980)
    folds = unique(df[, c("id", "physio_section")])
    folds$fold = createFolds(as.factor(folds$physio_section), k=k, list=F)
    df = merge(df, folds[, c('id', 'fold')], by='id')
    
    df$preds_kfold = NA
    df$preds_kfold0 = NA
    
    for (i in 1:k) {
        if (weights) {
            model = lm(call, df[df$fold != i, ], weights=obs)
        } else { model = lm(call, df[df$fold != i, ]) }
        
        df[df$fold == i, ]$preds_kfold = predict(model, df[df$fold == i, ])
        df[df$fold == i, ]$preds_kfold0 = predict(model, df0[df$fold == i, ])
    }
    
    return(df)
}

process_preds = function(df, df0, call, weights=T) {
    if (weights) {
        model = lm(call, df, weights=obs)
    } else {
        model = lm(call, df)
    }
    df$preds = predict(model, df)
    df$preds0 = predict(model, df0)
    df$preds[df$preds < 0] = 0
    df$preds0[df$preds0 < 0] = 0
    df$diff = df$preds - df$preds0
    df$perc_diff = df$diff/df$preds
    df$perc_diff[df$preds == 0] = 0
    summary(df[, c('pm', 'preds', 'preds0', 'diff', 'perc_diff')])
    return(df)
}

detach_all = function() {
    invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
}

theme_anne = function(font="Avenir", size=10) {
    theme_tufte(base_size=size, base_family=font) %+replace% 
        theme(
            panel.background  = element_blank(),
            plot.background = element_rect(fill="transparent", colour=NA), 
            axis.line.x = element_line(color="black", size = .2), 
            axis.line.y = element_line(color="black", size = .2), 
            plot.title = element_text(hjust = 0.5)
        )
}

theme_blank = function(font="Avenir", size=10) {
    theme_tufte(base_size=size, base_family=font) %+replace% 
        theme(
            panel.background  = element_blank(),
            plot.background = element_blank(), 
            plot.title = element_blank(), 
            axis.line=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.title=element_blank(),
            legend.position="none"
        )
}

find_common_bounds = function(x, y, square) {
    if (square) {
        all = c(x, y)
        ret = c(min(all, na.rm=T), max(all, na.rm=T), min(all, na.rm=T), max(all, na.rm=T))
    } else {
        ret = c(min(x, na.rm=T), max(x, na.rm=T), min(y, na.rm=T), max(y, na.rm=T))
    }
    return(ret)
}

panel = function(data, xbounds=NULL, ybounds=NULL, laby=NULL, labx=NULL, annotate=T, annotation=NULL, betal=F,
                 lm=T, a=0.2, w=NULL, square=F, font="Avenir", size=10, dot_size=T) {
    
    # save names for axes
    labels = names(data)
    names(data) = c("x", "y")
    
    # define bounds and annotation location
    if (is.null(xbounds)) {xbounds=data$x}
    if (is.null(ybounds)) {ybounds=data$y}
    bounds = find_common_bounds(xbounds, ybounds, square)
    minx = bounds[1]
    maxx = bounds[2]
    miny = bounds[3]
    maxy = bounds[4]
    if (is.null(labx)) {labx=minx} 
    if (is.null(laby)) {laby=maxy} 
    
    # plot!
    plot = ggplot(data, aes(x, y)) + xlab(labels[1]) + ylab(labels[2]) + 
        xlim(minx, maxx) + ylim(miny, maxy) + theme_anne(font=font, size=size)
    if (!is.null(w) & dot_size) {plot = plot+geom_point(aes(size=w), alpha = a)+guides(size=FALSE)
    } else {plot = plot+geom_point(alpha = a)}
    
    if (lm) {
        plot = plot + geom_smooth(method="lm", se=FALSE, alpha=0.4, size=.5)
        if (annotate) {
            # model data to get r2 and Beta for annotation
            if (!is.null(w)) {model = lm(y~x, weights=w, data=data)
            } else {model = lm(y~x, data=data)}
            r2 = format(round(summary(model)$r.squared, 2), nsmall = 2)
            beta = round(summary(model)$coefficients[2], 2)
            
            if (!is.null(annotation) & betal) {lab = paste0(annotation, "\nr^2 =='", r2, "; Beta = ", beta,"'")}
            else if (!is.null(annotation) & !betal) {lab = paste(annotation, "\nr^2 =='", r2, "'")}
            else if (betal) {lab = paste0("r^2 =='", r2, "; Beta = ", beta,"'")}
            else if (!betal) {lab = paste("r^2 =='", r2, "'")}
            
            #add annotation to plot
            plot = plot + annotate("text", x=labx, y=laby, label=lab, color="grey46", 
                                   hjust=0, vjust=1, size=size/2.5, family=font, parse=T)
        } 
    }
    
    if (!lm & annotate) {
        if(!is.null(w)) {
            r2 = round(wtd.cor(data$x, data$y, w=w)[1]^2, 2)
        } else {
            r2 = round(cor(data$x, data$y, use ="c")^2, 2)
        }
        lab = paste("r^2 == ", r2)
        if (!is.null(annotation)) {lab = paste(annotation, "\nr^2 == ", r2)}
        plot = plot + annotate("text", x=labx, y=laby, label=lab, color="grey46", 
                               hjust=0, vjust=1, size=size/2.5, family=font, parse=T)
    }
    
    print(lab)
    return(plot)
}

panel_set = function(data, name, subtitles=NULL, n=1, w=NULL, save_path=NULL, font="Avenir", 
                     lm=rep(T, length(data)), annotate=T, betal=F, a=0.2, square=F, size=10, unif_bounds=T, dot_size=T) {
    plots = list(rep(NA, length(data)))
    
    xbounds = unlist(sapply(data, function(x) x[, 1]))
    xbounds = c(min(xbounds, na.rm=T)-.1, max(xbounds, na.rm=T)+.1)
    ybounds = unlist(sapply(data, function(x) x[, 2]))
    ybounds = c(min(ybounds, na.rm=T)-.1, max(ybounds, na.rm=T)+.1)
    
    #for each set of x/y passed in 
    for (i in 1:length(data)) {
        #data[[i]] needs to be a df with two columns, named as they should be 
        cur = data[[i]]
        weights = w[[i]]
        if (!unif_bounds) {
            xbounds = c(min(cur[, 1], na.rm=T)-abs(min(cur[, 1], na.rm=T)*.1), 
                        max(cur[, 1], na.rm=T)+abs(min(cur[, 1], na.rm=T)*.1))
            ybounds = c(min(cur[, 2], na.rm=T)-abs(min(cur[, 2], na.rm=T)*.1),
                        max(cur[, 2], na.rm=T)+abs(min(cur[, 2], na.rm=T)*.1))
        }
        gg = panel(cur, xbounds, ybounds, w=weights, font=font, lm=lm[i], size=size,
                   annotate=annotate, betal=betal, a=a, square=square, dot_size=dot_size) + ggtitle(subtitles[i])
        plots[[i]] = gg
    }
    
    arrange_panels(plots, n, save_path, name, font)
}



# a function that moves the downloaded data from the Downloads folder and names
unzip_rename = function(name) {
  x = file.info(list.files("~/Downloads", full.names=T))
  x = rownames(x[x$mtime ==  max(x$mtime),])
  unzip(x, exdir="./data/airport")
  
  Sys.sleep(abs(rnorm(1, mean=.2, sd=0.1)))
  
  y = file.info(list.files("./data/airport", full.names=T))
  y = rownames(y[y$mtime ==  max(y$mtime),])
  file.rename(from=y, to=paste0("./data/airport/",  name, ".csv"))
  
  y = read.csv(paste0("./data/airport/",  name, ".csv"))
  saveRDS(y, paste0("./data/airport/",  name, ".RDS"))
  
  file.remove(paste0("./data/airport/",  name, ".csv"))
  file.remove(x)
  
  return(T)
}


# get the mean of the variables in adjacent cells for each grid

get_adj = function(x) {
  one = full_data[full_data$id %in% names(adj1[[x]]), ]
  one = one = one %>% group_by(year) %>% 
    dplyr::select(all_of(adjacency_vars)) %>% 
    dplyr::summarize_all(funs(mean_na))
  names(one)[-1] = paste0(names(one)[-1], "_adj1")
  one$id = names(adj1)[x]
  
  two = full_data[full_data$id %in% names(adj2[[x]]), ]
  two = two = two %>% group_by(year) %>% 
    dplyr::select(all_of(adjacency_vars)) %>% 
    dplyr::summarize_all(funs(mean_na))
  names(two)[-1] = paste0(names(two)[-1], "_adj2")
  two$id = names(adj2)[x]
  
  three = full_data[full_data$id %in% names(adj3[[x]]), ]
  three = three = three %>% group_by(year) %>% 
    dplyr::select(all_of(adjacency_vars)) %>% 
    dplyr::summarize_all(funs(mean_na))
  names(three)[-1] = paste0(names(three)[-1], "_adj3")
  three$id = names(adj3)[x]
  
  comb = cbind(one, two %>% dplyr::select(-id, -year), three %>% dplyr::select(-id, -year))
  return(comb)
}





interpolate_cases = function(data, column, dont_use=c("id", "pm", "physio_section", "obs")) {
  
  cols = apply(data[is.na(data[,column]), ], 2, anyNA) #any na in cols
  cols = unique(c(names(data)[!cols], column)) #cols with no na
  cols = cols[!cols %in% dont_use] #remove cols in dont use
  model_data = data[, cols] #select the cols not in dontuse and w no nas
  model_data = model_data[complete.cases(model_data), ] #get complete cases
  
  not_model = data[is.na(data[,column]), ]
  cols = cols[cols != column] #cols not equal to of interest
  types = sapply(data[, cols], FUN=class) #get numeric cols
  num_cols = cols[types != "character" & types  != "factor"]
  fac_cols = cols[types == "character" | types  == "factor"]
  
  # check to see if there are any missing levels in the factor columns
  for (i in fac_cols) {
    if( any(!unique(not_model[, i]) %in% unique(model_data[, i])) ) {
      
      # if yes, don't use them to interpolate the column of interest
      cols = cols[cols != i]
    }
  }
  
  # create a model to predict the col of interest, using all columns still passing 
  # checks and additionally using the log for numeric columns
  m = lm(paste0(column, "~", paste(cols, collapse="+"), "+", 
                paste(paste0("log(", num_cols, "+1)"), collapse="+")), 
         model_data)
  print(summary(m)$r.squared)
  
  # return the model
  return(m)
}



get_km_bounds = function(v, k, iter=15, n=3) {
  v = v[!is.na(v)]
  km = kmeans(v, 3, nstart=n, iter.max=iter)
  bounds = rep(NA, k)
  for (i in 1:k) {
    bounds[i] = max(v[km$cluster==i])
  }
  bounds = sort(bounds)[1:(k-1)]
  return(bounds)
}

