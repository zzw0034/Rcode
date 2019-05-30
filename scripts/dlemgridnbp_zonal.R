## Title: dlemgridnbp_zonal
## Author: hzs0087@auburn.edu 
## Date: 2018-08-24


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("G:/hao_shi/gcp18/submit/nc_result/")

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

## read nbp
gridcReadBin <- function(itime, scen) {
  library(raster)
 
  year <- itime - 1 + 1700
  
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
  
  c_r_m0 <- m0_s[["VegC_rec"]] + m0_s[["LitC_rec"]] + m0_s[["CWDC_rec"]] + m0_s[["SoilC_rec"]]
  c_r_m1 <- m1_s[["VegC_rec"]] + m1_s[["LitC_rec"]] + m1_s[["CWDC_rec"]] + m1_s[["SoilC_rec"]]
  c_r_m2 <- m2_s[["VegC_rec"]] + m2_s[["LitC_rec"]] + m2_s[["CWDC_rec"]] + m2_s[["SoilC_rec"]]
  c_r_m3 <- m3_s[["VegC_rec"]] + m3_s[["LitC_rec"]] + m3_s[["CWDC_rec"]] + m3_s[["SoilC_rec"]]
  
  ## combine all masks
  c_r <- c_r_m0 + c_r_m1 + c_r_m2 + c_r_m3
  c_df <- (c_r * landarea_r * 10^6) %>% as.data.frame(., xy=TRUE) %>% subset(., !is.na(layer))
  
  c_global_stat <- sum(c_df$layer) * 10^-15
  c_tropic_stat <- subset(c_df, y <= 30 & y >= -30)$layer %>% sum(.) * 10^-15
  c_north_stat <- subset(c_df, y > 30)$layer %>% sum(.) * 10^-15
  c_south_stat <- subset(c_df, y < -30)$layer %>% sum(.) * 10^-15
  
  return(c(year, c_global_stat, c_north_stat, c_tropic_stat, c_south_stat))
}

## parallel generation
for (scen in paste0("s", 0:4)) {
  cl <- makeCluster(16)
  clusterExport(cl, c("%>%", "readDlemGridpool", "readDlemBin", "landarea_r"))
  gridc_res <- parSapply(cl, 1:318, gridcReadBin, scen=scen)
  stopCluster(cl) 
  
  gridc_res <- t(gridc_res) %>% data.frame
  names(gridc_res) <- c("Year", "C_global", "C_north", "C_tropic", "C_south")
  write.table(gridc_res, paste0("DLEM_", toupper(scen), "_zonalC.txt"), quote=FALSE, sep=",",
              row.names=FALSE)
  
  gridnbp <- gridc_res[-1, ] - gridc_res[-nrow(gridc_res), ]
  gridnbp$Year <- gridc_res$Year[-1]
  names(gridnbp) <- c("Year", "NBP_global", "NBP_north", "NBP_tropic", "NBP_south")
  write.table(gridnbp, paste0("DLEM_", toupper(scen), "_zonalNBP.txt"), quote=FALSE, sep=",",
              row.names=FALSE)
  
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





