## Title: dlemgridpool2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/isimip/Rcode/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/isimip/nc_result")


## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:93
time_dim <- ncdim_def("time", "years since 2006-01-01 00:00:00", time, calendar="noleap")

fill_value <- 1.0*10^20

gridfluxReadBin <- function(itime) { ## itime starts from 1 to 94
  library(raster)
  
  gpp_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  npp_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  ra_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  rh_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  # nep_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  et_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  mrro_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  
  itime <- itime - 1
  # mon <- itime %% 12
  year <- itime - 1 + 2006
  
  m0_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/isimip/output/miroc5/rcp26/gridflux/",
                                    "yy", year, ".dlem"))

  gpp_r <- m0_s[["g_gpp"]]; 
  npp_r <- m0_s[["g_npp"]];
  ra_r <- gpp_r - npp_r
  rh_r <- m0_s[["g_rh"]]; 
  # nep_r_m0 <- npp_r_m0 - rh_r_m0; 
  et_r <- m0_s[["g_evap"]] + m0_s[["g_trans"]]
  mrro_r <- m0_s[["g_runoff_drain"]] + m0_s[["g_runoff_surf"]]
  #?qtot?
  
  ## write to array
  gpp_arr[,] <- gpp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  npp_arr[,] <- npp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  ra_arr[,] <- ra_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  rh_arr[,] <- rh_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  # nep_arr[,] <- nep_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  et_arr[,] <- et_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  mrro_arr[,] <- mrro_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))

  return(list(gpp_arr, npp_arr, ra_arr, rh_arr, et_arr, mrro_arr))
}

## parallel generation
# scen <- "s1"

cl <- makeCluster(5)
clusterExport(cl, c("%>%", "lon", "lat", "readDlemGridflux", "readDlemBin"))
gridflux_res <- parLapply(cl, 1:94, gridfluxReadBin)
stopCluster(cl)

##define variables
var_lnames <- c(gpp="Gross Primary Production", npp="Net Primary Production",
                ra="Autotrophic (Plant) Respiration", rh="Heterotrophic Respiration",
                evapotrans="Total Evapo-Transpiration",
                mrro="Total Runoff")

var_units <- c(gpp="kg C m-2", npp="kg C m-2",
                ra="kg C m-2", rh="kg C m-2",
                evapotrans="kg m-2",
                mrro="kg m-2")

vars <- c("gpp", "npp", "ra", "rh", "evapotrans", "mrro")
var_orders <- c(gpp=1, npp=2, ra=3, rh=4, evapotrans=5, mrro=6)
var_coeffs <- c(gpp=10^-3, npp=10^-3, ra=10^-3, rh=10^-3, evapotrans=1, mrro=1)
# days_mon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

## function for writing ncdf
writeNC <- function(var) {
  library(ncdf4)
  
  var_def <- ncvar_def(var, var_units[var], list(lon_dim, lat_dim, time_dim), fill_value,
                     var_lnames[var], "float", compression=9)  
  
  nc <- nc_create(paste0("DLEM_miroc5_ewembi_rcp26_rcp26soc_co2_",var,"_global_annual_2006_2099",
                          ".nc"), list(var_def), force_v4=TRUE)
  
  gridfluxwrite_res <- lapply(1:94, FUN=function(itime) {
    # mon <- (itime - 1) %% 12 + 1
    # secs <- days_mon[mon] * 24 * 3600
 
    ncvar_put(nc, var_def, gridflux_res[[itime]][[var_orders[var]]] * var_coeffs[var],
            start=c(1, 1, itime), count=c(-1, -1, 1))
    
    cat("itime: ", itime, "\n")
  })
  
  ## put additional attributes
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  
  ## add global atts
  ncatt_put(nc, 0, "title", "DLEM output for isimip2B, 2019")
  ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")
  
  ## close netcdf
  nc_close(nc)
  
  return(0)
}

## parallel writing
cl <- makeCluster(5)
clusterExport(cl, c("%>%", "lon", "lat", "var_units", "var_orders", "var_coeffs", 
                    "var_lnames", "gridflux_res", "fill_value", "lon_dim", "lat_dim", "time_dim"))
nc_res <- parLapply(cl, vars, writeNC)
stopCluster(cl)





