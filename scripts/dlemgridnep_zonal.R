## Title: dlemgridnep_zonal
## Author: hzs0087@auburn.edu 
## Date: 2018-08-24


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/nc_result/")

## calculate land area 
landmask <- new("GDALReadOnlyDataset",
                "S:/hao_shi/gcp18/trendyv7/mask/mask_new/w001001.adf")
landmask_r <- asSGDF_GROD(landmask)
gridded(landmask_r) <- TRUE
landmask_r <- raster(landmask_r)
proj4string(landmask_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

flake_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_lake.bin")
focean_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_ocean.bin")
friver_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_river.bin")

fland_r <- landmask_r - flake_r - focean_r - friver_r
landarea_r <- area(fland_r) * fland_r

## read nep
gridnepReadBin <- function(itime, scen) {
  library(raster)
 
  year <- itime - 1 + 1700
  
  m0_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m0",
                                  scen, "/m0", scen, "yy", year, ".dlem"))
  
  m1_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m1",
                                  scen, "/m1", scen, "yy", year, ".dlem"))
  
  m2_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m2",
                                  scen, "/m2", scen, "yy", year, ".dlem"))
  
  m3_s <- readDlemGridflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                           paste0("S:/hao_shi/gcp18/trendyv7/output/gridfluxes/m3",
                                  scen, "/m3", scen, "yy", year, ".dlem"))
  
  nep_r_m0 <- m0_s[["g_npp"]] - m0_s[["g_rh"]]
  nep_r_m1 <- m1_s[["g_npp"]] - m1_s[["g_rh"]]
  nep_r_m2 <- m2_s[["g_npp"]] - m2_s[["g_rh"]]
  nep_r_m3 <- m3_s[["g_npp"]] - m3_s[["g_rh"]]
  
  ## combine all masks
  nep_r <- nep_r_m0 + nep_r_m1 + nep_r_m2 + nep_r_m3
  nep_df <- (nep_r * landarea_r * 10^6) %>% as.data.frame(., xy=TRUE) %>% subset(., !is.na(layer))
  
  nep_global_stat <- sum(nep_df$layer) * 10^-15
  nep_tropic_stat <- subset(nep_df, y <= 30 & y >= -30)$layer %>% sum(.) * 10^-15
  nep_north_stat <- subset(nep_df, y > 30)$layer %>% sum(.) * 10^-15
  nep_south_stat <- subset(nep_df, y < -30)$layer %>% sum(.) * 10^-15
  
  return(c(year, nep_global_stat, nep_north_stat, nep_tropic_stat, nep_south_stat))
}

## parallel generation
for (scen in paste0("s", 0:4)) {
  cl <- makeCluster(16)
  clusterExport(cl, c("%>%", "readDlemGridflux", "readDlemBin", "landarea_r"))
  gridnep_res <- parSapply(cl, 1:318, gridnepReadBin, scen=scen)
  stopCluster(cl) 
  
  gridnep_res <- t(gridnep_res)
  names(gridnep_res) <- c("Year", "NEP_global", "NEP_north", "NEP_tropic", "NEP_south")
  write.table(gridnep_res, paste0("DLEM_", toupper(scen), "_zonalNEP.txt"), quote=FALSE, sep=",",
              row.names=FALSE, header=TRUE)
  cat("Scen: ", scen, "\n")
}


##define variables
var_lnames <- c(landArea="land area per grid")
var_units <- c(landArea="km2")
vars <- c("landArea")
fill_value <- -99999.0

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

## function for writing ncdf
writeNC <- function(var) {
  # library(ncdf4)
  
  var_def <- ncvar_def(var, var_units[var], list(lon_dim, lat_dim), fill_value,
                       var_lnames[var], "float", compression=9)  
  
  nc <- nc_create(paste0("DLEM_", var, ".nc"), list(var_def), force_v4=TRUE)
  ncvar_put(nc, var_def, landarea_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720)))
  
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
res <- writeNC("landArea")





