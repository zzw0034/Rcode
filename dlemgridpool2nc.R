## Title: dlemgridpool2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/isimip/output/miroc5/nc_result/")


## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 0:93
time_dim <- ncdim_def("time", "years since 2006-01-01 00:00:00", time, calendar="noleap")

fill_value <- 1.0*10^20

gridpoolReadBin <- function(itime) { ## itime starts from 1 to 94
  library(raster)
 
  litc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
  soilc_arr <- array(numeric(), dim=c(length(lon), length(lat)))
 
  
  year <- itime - 1 + 2006
  
  m0_s <- readDlemGridpool("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("S:/hao_shi/isimip/output/miroc5/rcp26/gridpool/","yy", year, ".dlem"))
  

  
 
  litc_r_m0 <- m0_s[["LitC_rec"]]; 
  soilc_r_m0 <- m0_s[["SoilC_rec"]]; 

  
  litc_r <- litc_r_m0 
  soilc_r <- soilc_r_m0 
      
  ## write to array
  litc_arr[,] <- litc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
  soilc_arr[,] <- soilc_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))

  return(list(litc_arr, soilc_arr))
}


##define variables


litc_def <- ncvar_def("cLitter", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                     "Carbon in Litter Pool", "float", compression=9)

soilc_def <- ncvar_def("cSoil", "kg C m-2", list(lon_dim, lat_dim, time_dim), fill_value,
                     "Carbon in Soil Pool", "float", compression=9)


## create netcdf
scen="rcp26"
nc_litc <- nc_create(paste0("DLEM_miroc5_ewembi_rcp26_rcp26soc_co2_clitter_global_annual_2006_2099", ".nc")
                     , list(litc_def), force_v4=TRUE)
nc_soilc <- nc_create(paste0("DLEM_miroc5_ewembi_rcp26_rcp26soc_co2_csoil_global_annual_2006_2099", ".nc")
                     , list(soilc_def), force_v4=TRUE)


## parallel generation
cl <- makeCluster(4)
clusterExport(cl, c("%>%", "lon", "lat", "readDlemGridpool","readDlemBin"))
gridpool_res <- parLapply(cl, 1:94, gridpoolReadBin)
stopCluster(cl)


gridpoolWrite_res <- lapply(1:94, FUN=function(itime) {
  
  ncvar_put(nc_litc, litc_def, gridpool_res[[itime]][[1]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  ncvar_put(nc_soilc, soilc_def, gridpool_res[[itime]][[2]] * 10^-3,
            start=c(1, 1, itime), count=c(-1, -1, 1))
  
  cat("itime: ", itime, "\n")
})


## put additional attributes
for(nc_var in c( "nc_litc", "nc_soilc")) {
  nc <- eval(parse(text=nc_var))
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  
  ## add global atts
  ncatt_put(nc, 0, "title", "DLEM output for isimip2B, 2019")
  ncatt_put(nc, 0, "institution", "Auburn University, US")
  ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")
  
  ## close netcdf
  nc_close(nc)
}

rm(gridpool_res)





