## Title: dlemgridnbp2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-24


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/submit/nc_result/")

## read nbp
gridcReadBin <- function(itime, scen) {
  library(raster)
  
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
  
  c_r_m0 <- m0_s[["VegC_rec"]] + m0_s[["LitC_rec"]] + m0_s[["SoilC_rec"]]
  c_r_m1 <- m1_s[["VegC_rec"]] + m1_s[["LitC_rec"]] + m1_s[["SoilC_rec"]]
  c_r_m2 <- m2_s[["VegC_rec"]] + m2_s[["LitC_rec"]] + m2_s[["SoilC_rec"]]
  c_r_m3 <- m3_s[["VegC_rec"]] + m3_s[["LitC_rec"]] + m3_s[["SoilC_rec"]]
  
  ## combine all masks
  c_r <- c_r_m0 + c_r_m1 + c_r_m2 + c_r_m3
  return(c_r)
}

## parallel generation
scen <- "s1"

cl <- makeCluster(16)
clusterExport(cl, c("%>%", "readDlemGridpool", "readDlemBin"))
gridc_res <- parLapply(cl, 1:1416, gridcReadBin, scen=scen)
stopCluster(cl) 

##define variables
var_lnames <- c(nbp="Net Biospheric Production")
var_units <- c(nbp="kg C m-2 s-1")
vars <- c("nbp")
fill_value <- -99999.0

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:1415
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

days_mon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

## function for writing ncdf
writeNC <- function(var) {
  # library(ncdf4)
  
  var_def <- ncvar_def(var, var_units[var], list(lon_dim, lat_dim, time_dim), fill_value,
                       var_lnames[var], "float", compression=9)  
  
  nc <- nc_create(paste0("DLEM_", toupper(scen), "_", var, ".nc"), list(var_def), force_v4=TRUE)
  
  for(itime in 1:length(gridc_res)) {
    if(itime==1) {
      nbp_r <- raster(gridc_res[[itime]])
      nbp_r[] <- NA
    } else {
      nbp_r <- gridc_res[[itime]] - gridc_res[[itime-1]]
    }
    
    mon <- (itime - 1) %% 12 + 1
    secs <- days_mon[mon] * 24 * 3600
    
    ncvar_put(nc, var_def, (nbp_r * 10^-3 / secs)  %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720)),
              start=c(1, 1, itime), count=c(-1, -1, 1))
    cat("itime: ", itime, "\n")
  }
  
  
  ## put additional attributes
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  
  ## add global atts
  ncatt_put(nc, 0, "title", "DLEM land area for TRENDYv7, 2018")
  ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")
  
  ## close netcdf
  nc_close(nc)
  
  return(0)
}

## nc writing
res <- writeNC("nbp")





