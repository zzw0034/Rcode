## Title: dlembiomestate2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-09-21


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/submit/nc_result/")

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:1415
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

biomestateReadBin <- function(itime, scen) { ## itime starts from 1 to 118
  library(raster)
  vegc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  itime <- itime - 1
  mon <- itime %% 12
  year <- itime %/% 12 + 1900
  
  m0_s <- readDlemBiomestate("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/biomestate/m0",
                                  scen, "/m0", scen, "yy", year, "m", mon, ".dlem"))
  
  m1_s <- readDlemBiomestate("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/biomestate/m1",
                                  scen, "/m1", scen, "yy", year, "m", mon, ".dlem"))
  
  m2_s <- readDlemBiomestate("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/biomestate/m2",
                                  scen, "/m2", scen, "yy", year, "m", mon, ".dlem"))
  
  m3_s <- readDlemBiomestate("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/biomestate/m3",
                                  scen, "/m3", scen, "yy", year, "m", mon, ".dlem"))
  
  leafc_r_m0 <- m0_s[["leaf_c_91"]]; sapwc_r_m0 <- m0_s[["sapw_c_91"]]; 
  heartwc_r_m0 <- m0_s[["heartw_c_91"]]; reprodc_r_m0 <- m0_s[["reprod_c_91"]];
  
  leafc_r_m1 <- m1_s[["leaf_c_91"]]; sapwc_r_m1 <- m1_s[["sapw_c_91"]]; 
  heartwc_r_m1 <- m1_s[["heartw_c_91"]]; reprodc_r_m1 <- m1_s[["reprod_c_91"]];
  
  leafc_r_m2 <- m2_s[["leaf_c_91"]]; sapwc_r_m2 <- m2_s[["sapw_c_91"]]; 
  heartwc_r_m2 <- m2_s[["heartw_c_91"]]; reprodc_r_m2 <- m2_s[["reprod_c_91"]];
  
  leafc_r_m3 <- m3_s[["leaf_c_91"]]; sapwc_r_m3 <- m3_s[["sapw_c_91"]]; 
  heartwc_r_m3 <- m3_s[["heartw_c_91"]]; reprodc_r_m3 <- m3_s[["reprod_c_91"]];
  
  ## combine all masks
  leafc_r <- leafc_r_m0 + leafc_r_m1 + leafc_r_m2 + leafc_r_m3
  sapwc_r <- sapwc_r_m0 + sapwc_r_m1 + sapwc_r_m2 + sapwc_r_m3
  heartwc_r <- heartwc_r_m0 + heartwc_r_m1 + heartwc_r_m2 + heartwc_r_m3
  reprodc_r <- reprodc_r_m0 + reprodc_r_m1 + reprodc_r_m2 + reprodc_r_m3
  
  vegc_r <- leafc_r + sapwc_r + heartwc_r + reprodc_r
  
  ## write to array
  vegc_arr[,] <- vegc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  
  return(vegc_arr)
}


##define variables
vegc_def <- ncvar_def("cVeg", "g C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                      "Carbon in Crops", "float", compression=9)

## create netcdf
scen="s3"
nc_vegc <- nc_create(paste0("DLEM_", toupper(scen), "_cropVeg.nc"), list(vegc_def), force_v4=TRUE)

## parallel generation
cl <- makeCluster(16)
clusterExport(cl, c("%>%", "lon", "lat", "readDlemBiomestate","readDlemBin"))
biomestate_res <- parLapply(cl, 1:1416, biomestateReadBin, scen=scen)
stopCluster(cl)


biomestateWrite_res <- lapply(1:1416, FUN=function(itime) {
  ncvar_put(nc_vegc, vegc_def, biomestate_res[[itime]],
            start=c(1, 1, itime), count=c(-1, -1, 1))
  cat("itime: ", itime, "\n")
})


## put additional attributes
for(nc_var in c("nc_vegc")) {
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

rm(biomestate_res)





