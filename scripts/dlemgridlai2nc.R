## Title: dlemlaimonthly2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-24


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/submit/nc_result/")

## read in non-veg luc
## convert to binary format
# landmask <- new("GDALReadOnlyDataset", 
#                 "S:/hao_shi/gcp18/trendyv7/mask/mask_new/w001001.adf")
# landmask_r <- asSGDF_GROD(landmask)
# gridded(landmask_r) <- TRUE
# landmask_r <- raster(landmask_r)
# proj4string(landmask_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:1415
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

gridlaiReadBin <- function(itime, scen) { ## itime starts from 1 to 118
  library(raster)
  lai_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  itime <- itime - 1
  mon <- itime %% 12
  year <- itime %/% 12 + 1900
  
  m0_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m0",
                                  scen, "/m0", scen, "yy", year, "m", mon, ".dlem"))
  
  m1_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m1",
                                  scen, "/m1", scen, "yy", year, "m", mon, ".dlem"))
  
  m2_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m2",
                                  scen, "/m2", scen, "yy", year, "m", mon, ".dlem"))
  
  m3_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m3",
                                  scen, "/m3", scen, "yy", year, "m", mon, ".dlem"))
  
  lai_r_m0 <- m0_s[["LAIgrid_rec"]]; lai_r_m1 <- m1_s[["LAIgrid_rec"]]
  lai_r_m2 <- m2_s[["LAIgrid_rec"]]; lai_r_m3 <- m3_s[["LAIgrid_rec"]]
  
  ## combine all masks
  lai_r <- lai_r_m0 + lai_r_m1 + lai_r_m2 + lai_r_m3
 
  ## previous year 
  m0_s_prev <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m0",
                                  scen, "/m0", scen, "yy", year-1, "m", mon, ".dlem"))
  
  m1_s_prev <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m1",
                                  scen, "/m1", scen, "yy", year-1, "m", mon, ".dlem"))
  
  m2_s_prev <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m2",
                                  scen, "/m2", scen, "yy", year-1, "m", mon, ".dlem"))
  
  m3_s_prev <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m3",
                                  scen, "/m3", scen, "yy", year-1, "m", mon, ".dlem"))
  
  lai_r_m0_prev <- m0_s_prev[["LAIgrid_rec"]]; lai_r_m1_prev <- m1_s_prev[["LAIgrid_rec"]]
  lai_r_m2_prev <- m2_s_prev[["LAIgrid_rec"]]; lai_r_m3_prev <- m3_s_prev[["LAIgrid_rec"]]
  
  ## combine all masks
  lai_r_prev <- lai_r_m0_prev + lai_r_m1_prev + lai_r_m2_prev + lai_r_m3_prev
  
  ## previous year
  
  lai_r <- lai_r - lai_r_prev
  ## write to array
  lai_arr[,] <- lai_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  
  return(lai_arr)
}


##define variables
lai_def <- ncvar_def("lai", "m2 m-2 (projected leaf area per grid)",
                     list(lon_dim, lat_dim, time_dim), fill_value,
                      "Leaf Area Index", "float", compression=9)

## create netcdf
scen="s0"
nc_lai <- nc_create(paste0("DLEM_",toupper(scen),"_lai.nc"), list(lai_def), force_v4=TRUE)

## parallel generation
cl <- makeCluster(6) ## note: can increase the number of cores here, e.g., 16
clusterExport(cl, c("%>%", "lon", "lat", "readDlemGridpool","readDlemBin"))
gridlai_res <- parLapply(cl, 1:1416, gridlaiReadBin, scen=scen)
stopCluster(cl)


gridlaiWrite_res <- lapply(1:1416, FUN=function(itime) {
  ncvar_put(nc_lai, lai_def, gridlai_res[[itime]],
            start=c(1, 1, itime), count=c(-1, -1, 1))
  cat("itime: ", itime, "\n")
})


## put additional attributes
for(nc_var in c("nc_lai")) {
  nc <- eval(parse(text=nc_var))
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  
  ## add global atts
  ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
  ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")
  
  ## close netcdf
  nc_close(nc)
}

rm(gridlai_res)





