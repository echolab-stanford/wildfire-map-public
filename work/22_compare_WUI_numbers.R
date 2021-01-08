source("00_functions.R")

########################################################################################
# Based on larger pipeline by: Jenny Xue
# Last edited by: Anne, Oct 2020
########################################################################################

#choose the locations to sample based on the grid we use for predictions
set.seed(983451034)
grid = readRDS("data/clean/national_grid.RDS")
grid = grid[sample(nrow(grid), 50),  ]

# read in the land cover data and only keep the area sampled
nlcd = raster ('data/wui/nlcd_2001_Land_Cover_L48_20190424.img')
grid = spTransform(grid, crs(nlcd))
nlcd = crop(nlcd, grid)  # very  slow
nlcd = mask(nlcd, grid)
