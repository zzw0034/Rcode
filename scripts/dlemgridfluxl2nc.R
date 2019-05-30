## Title: dlemgridpool2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("G:/hao_shi/")

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

gridfluxReadBin <- function(itime, scen) { ## itime starts from 1 to 118
  library(raster)
  
  gpp_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  npp_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  ra_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  rh_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  nep_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  et_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  mrro_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  itime <- itime - 1
  mon <- itime %% 12
  year <- itime %/% 12 + 1900
  
  m0_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m0",
                                   scen, "/m0", scen, "yy", year, "m", mon, ".dlem"))
  
  m1_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m1",
                                   scen, "/m1", scen, "yy", year, "m", mon, ".dlem"))
  
  m2_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m2",
                                   scen, "/m2", scen, "yy", year, "m", mon, ".dlem"))
  
  m3_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m3",
                                   scen, "/m3", scen, "yy", year, "m", mon, ".dlem"))
  
  gpp_r_m0 <- m0_s[["g_gpp"]]; npp_r_m0 <- m0_s[["g_npp"]]; ra_r_m0 <- gpp_r_m0 - npp_r_m0
  rh_r_m0  <- m0_s[["g_rh"]]; nep_r_m0 <- npp_r_m0 - rh_r_m0; 
  et_r_m0 <- m0_s[["g_evap"]] + m0_s[["g_trans"]]
  mrro_r_m0 <- m0_s[["g_runoff_drain"]] + m0_s[["g_runoff_surf"]]
  
  gpp_r_m1 <- m1_s[["g_gpp"]]; npp_r_m1 <- m1_s[["g_npp"]]; ra_r_m1 <- gpp_r_m1 - npp_r_m1
  rh_r_m1  <- m1_s[["g_rh"]]; nep_r_m1 <- npp_r_m1 - rh_r_m1; 
  et_r_m1 <- m1_s[["g_evap"]] + m1_s[["g_trans"]]
  mrro_r_m1 <- m1_s[["g_runoff_drain"]] + m1_s[["g_runoff_surf"]]
  
  gpp_r_m2 <- m2_s[["g_gpp"]]; npp_r_m2 <- m2_s[["g_npp"]]; ra_r_m2 <- gpp_r_m2 - npp_r_m2
  rh_r_m2  <- m2_s[["g_rh"]]; nep_r_m2 <- npp_r_m2 - rh_r_m2; 
  et_r_m2 <- m2_s[["g_evap"]] + m2_s[["g_trans"]]
  mrro_r_m2 <- m2_s[["g_runoff_drain"]] + m2_s[["g_runoff_surf"]]
  
  gpp_r_m3 <- m3_s[["g_gpp"]]; npp_r_m3 <- m3_s[["g_npp"]]; ra_r_m3 <- gpp_r_m3 - npp_r_m3
  rh_r_m3  <- m3_s[["g_rh"]]; nep_r_m3 <- npp_r_m3 - rh_r_m3; 
  et_r_m3 <- m3_s[["g_evap"]] + m3_s[["g_trans"]]
  mrro_r_m3 <- m3_s[["g_runoff_drain"]] + m3_s[["g_runoff_surf"]]
  
  
  ## combine all masks
  gpp_r <- gpp_r_m0 + gpp_r_m1 + gpp_r_m2 + gpp_r_m3
  npp_r <- npp_r_m0 + npp_r_m1 + npp_r_m2 + npp_r_m3
  ra_r <- ra_r_m0 + ra_r_m1 + ra_r_m2 + ra_r_m3
  rh_r <- rh_r_m0 + rh_r_m1 + rh_r_m2 + rh_r_m3
  nep_r <- nep_r_m0 + nep_r_m1 + nep_r_m2 + nep_r_m3
  et_r <- et_r_m0 + et_r_m1 + et_r_m2 + et_r_m3
  mrro_r <- mrro_r_m0 + mrro_r_m1 + mrro_r_m2 + mrro_r_m3
  
  ## write to array
  gpp_arr[,] <- gpp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  npp_arr[,] <- npp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  ra_arr[,] <- ra_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  rh_arr[,] <- rh_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  nep_arr[,] <- nep_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  et_arr[,] <- et_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  mrro_arr[,] <- mrro_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))

  return(list(gpp_arr, npp_arr, ra_arr, rh_arr, nep_arr, et_arr, mrro_arr))
}

## parallel generation
scen <- "s1"

cl <- makeCluster(16)
clusterExport(cl, c("%>%", "lon", "lat", "readDlemGridflux", "readDlemBin"))
gridflux_res <- parLapply(cl, 1:1416, gridfluxReadBin, scen=scen)
stopCluster(cl)

##define variables
var_lnames <- c(gpp="Gross Primary Production", npp="Net Primary Production",
                ra="Autotrophic (Plant) Respiration", rh="Heterotrophic Respiration",
                nep="Net Ecosystem Production", evapotrans="Total Evapo-Transpiration",
                mrro="Total Runoff")

var_units <- c(gpp="kg C m-2 s-1", npp="kg C m-2 s-1",
                ra="kg C m-2 s-1", rh="kg C m-2 s-1",
                nep="kg C m-2 s-1", evapotrans="kg m-2 s-1",
                mrro="kg m-2 s-1")

vars <- c("gpp", "npp", "ra", "rh", "nep", "evapotrans", "mrro")
var_orders <- c(gpp=1, npp=2, ra=3, rh=4, nep=5, evapotrans=6, mrro=7)
var_coeffs <- c(gpp=10^-3, npp=10^-3, ra=10^-3, rh=10^-3, nep=10^-3, evapotrans=1, mrro=1)
days_mon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

## function for writing ncdf
writeNC <- function(var) {
  library(ncdf4)
  
  var_def <- ncvar_def(var, var_units[var], list(lon_dim, lat_dim, time_dim), fill_value,
                     var_lnames[var], "float", compression=9)  
  
  nc <- nc_create(paste0("DLEM_", toupper(scen), "_", var, ".nc"), list(var_def), force_v4=TRUE)
  
  gridfluxwrite_res <- lapply(1:1416, FUN=function(itime) {
    mon <- (itime - 1) %% 12 + 1
    secs <- days_mon[mon] * 24 * 3600
 
    ncvar_put(nc, var_def, gridflux_res[[itime]][[var_orders[var]]] * var_coeffs[var] / secs,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  })
  
  ## put additional attributes
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  
  ## add global atts
  ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
  ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")
  
  ## close netcdf
  nc_close(nc)
  
  return(0)
}

## parallel writing
cl <- makeCluster(4)
clusterExport(cl, c("%>%", "lon", "lat", "var_units", "var_orders", "var_coeffs", "days_mon",
                    "var_lnames", "scen", "gridflux_res", "fill_value", "lon_dim", "lat_dim", "time_dim"))
nc_res <- parLapply(cl, vars, writeNC)
stopCluster(cl)





