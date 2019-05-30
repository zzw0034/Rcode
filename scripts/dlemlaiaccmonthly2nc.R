## Title: dlemlaiaccmonthly2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-09-05


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/submit/nc_result/")

gridlaiReadBin <- function(itime, scen) { ## itime starts from 1 to 118
  library(raster)
  
  itime <- itime - 1
  mon <- itime %% 12
  year <- itime %/% 12 + 1899
  
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
  e <- extent(-180, 180, -90, 90)
  lai_r <- extend(lai_r, e)
  
  
  return(lai_r)
}

lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:1415
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

lai_def <- ncvar_def("lai", "m2 m-2 (projected leaf area per grid)",
                     list(lon_dim, lat_dim, time_dim), fill_value,
                     "Leaf Area Index", "float", compression=9)
lai_arr <- array(numeric(), dim=c(length(lon), length(lat)))

for(scen in paste0("s", 0:4)) {
  cl <- makeCluster(12) ## note: can increase the number of cores here, e.g., 16
  clusterExport(cl, c("%>%", "readDlemGridpool","readDlemBin"))
  lai_s <- parLapply(cl, 1:12, gridlaiReadBin, scen=scen)
  stopCluster(cl)
  
  lai_s <- do.call(stack, lai_s)
  lai_acc <- stack(paste0("S:/hao_shi/gcp18/submit/nc_result/DLEM_", toupper(scen), "_lai_acc.nc"))
  lai_acc <- stack(lai_s, lai_acc)
  
  lai_def <- ncvar_def("lai", "m2 m-2 (projected leaf area per grid)",
                       list(lon_dim, lat_dim, time_dim), fill_value,
                       "Leaf Area Index", "float", compression=9)
  nc_lai <- nc_create(paste0("DLEM_",toupper(scen),"_lai.nc"), list(lai_def), force_v4=TRUE)
  
  lai_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  for(i in 13:nlayers(lai_acc)) {
    lai_r <- lai_acc[[i]] - lai_acc[[i-12]]
    lai_arr <- lai_r %>% values
    
    ncvar_put(nc_lai, lai_def, lai_arr, start=c(1, 1, i-12), count=c(-1, -1, 1))
    cat(i, "\n")
  }
}

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







