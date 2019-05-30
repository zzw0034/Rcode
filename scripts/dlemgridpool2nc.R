## Title: dlemgridpool2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


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

time <- 0:117
time_dim <- ncdim_def("time", "years since 1900-01-01 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

gridpoolReadBin <- function(itime, scen) { ## itime starts from 1 to 118
  library(raster)
  vegc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  litc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  soilc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  cwdc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  year <- itime - 1 + 1900
  
  m0_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m0",
                                   scen, "/m0", scen, "yy", year, ".dlem"))
  
  m1_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m1",
                                   scen, "/m1", scen, "yy", year, ".dlem"))
  
  m2_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m2",
                                   scen, "/m2", scen, "yy", year, ".dlem"))
  
  m3_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridpool/m3",
                                   scen, "/m3", scen, "yy", year, ".dlem"))
  
  vegc_r_m0 <- m0_s[["VegC_rec"]]; litc_r_m0 <- m0_s[["LitC_rec"]]; 
  soilc_r_m0 <- m0_s[["SoilC_rec"]]; cwdc_r_m0 <- m0_s[["CWDC_rec"]];
  
  vegc_r_m1 <- m1_s[["VegC_rec"]]; litc_r_m1 <- m1_s[["LitC_rec"]]; 
  soilc_r_m1 <- m1_s[["SoilC_rec"]]; cwdc_r_m1 <- m1_s[["CWDC_rec"]];
  
  vegc_r_m2 <- m2_s[["VegC_rec"]]; litc_r_m2 <- m2_s[["LitC_rec"]]; 
  soilc_r_m2 <- m2_s[["SoilC_rec"]]; cwdc_r_m2 <- m2_s[["CWDC_rec"]];
  
  vegc_r_m3 <- m3_s[["VegC_rec"]]; litc_r_m3 <- m3_s[["LitC_rec"]];
  soilc_r_m3 <- m3_s[["SoilC_rec"]]; cwdc_r_m3 <- m3_s[["CWDC_rec"]];
  
  ## combine all masks
  vegc_r <- vegc_r_m0 + vegc_r_m1 + vegc_r_m2 + vegc_r_m3
  litc_r <- litc_r_m0 + litc_r_m1 + litc_r_m2 + litc_r_m3
  soilc_r <- soilc_r_m0 + soilc_r_m1 + soilc_r_m2 + soilc_r_m3
  cwdc_r <- cwdc_r_m0 + cwdc_r_m1 + cwdc_r_m2 + cwdc_r_m3
      
  ## write to array
  vegc_arr[,] <- vegc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  litc_arr[,] <- litc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  soilc_arr[,] <- soilc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  cwdc_arr[,] <- cwdc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))

  return(list(vegc_arr, litc_arr, soilc_arr, cwdc_arr))
}


##define variables
vegc_def <- ncvar_def("cVeg", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                     "Carbon in Vegetation", "float", compression=9)

litc_def <- ncvar_def("cLitter", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                     "Carbon in Litter Pool", "float", compression=9)

soilc_def <- ncvar_def("cSoil", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                     "Carbon in Soil Pool", "float", compression=9)

cwdc_def <- ncvar_def("CWDC", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                      "Carbon in Coarse Wood Debris", "float", compression=9)

## create netcdf
scen="s0"
nc_vegc <- nc_create(paste0("DLEM_", toupper(scen), "_cVeg.nc"), list(vegc_def), force_v4=TRUE)
nc_litc <- nc_create(paste0("DLEM_", toupper(scen), "_cLitter.nc"), list(litc_def), force_v4=TRUE)
nc_soilc <- nc_create(paste0("DLEM_", toupper(scen), "_cSoil.nc"), list(soilc_def), force_v4=TRUE)
nc_cwdc <- nc_create(paste0("DLEM_", toupper(scen), "_CWDC.nc"), list(cwdc_def), force_v4=TRUE)

## parallel generation
cl <- makeCluster(16)
clusterExport(cl, c("%>%", "lon", "lat", "readDlemGridpool","readDlemBin"))
gridpool_res <- parLapply(cl, 1:118, gridpoolReadBin, scen=scen)
stopCluster(cl)


gridpoolWrite_res <- lapply(1:118, FUN=function(itime) {
  ncvar_put(nc_vegc, vegc_def, gridpool_res[[itime]][[1]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  ncvar_put(nc_litc, litc_def, gridpool_res[[itime]][[2]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  ncvar_put(nc_soilc, soilc_def, gridpool_res[[itime]][[3]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  ncvar_put(nc_cwdc, cwdc_def, gridpool_res[[itime]][[4]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  cat("itime: ", itime, "\n")
})


## put additional attributes
for(nc_var in c("nc_vegc", "nc_litc", "nc_soilc", "nc_cwdc")) {
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

rm(gridpool_res)





