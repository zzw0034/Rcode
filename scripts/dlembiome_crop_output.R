## Title: dlembiome2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/result")

## read in non-veg luc
## convert to binary format
landmask <- new("GDALReadOnlyDataset", 
                "S:/hao_shi/gcp18/trendyv7/mask/mask_new/w001001.adf")
landmask_r <- asSGDF_GROD(landmask)
gridded(landmask_r) <- TRUE
landmask_r <- raster(landmask_r)
proj4string(landmask_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

landvalue_r <- landmask_r - 1

####################################
biomeReadBin <- function(year, scen) {
  library(raster)
  
  m0_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("D:/biomefluxes/m0",
                                   scen, "/m0", scen, "yy", year, ".dlem"))
  
  m1_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("D:/biomefluxes/m1",
                                   scen, "/m1", scen, "yy", year, ".dlem"))
  
  m2_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("D:/biomefluxes/m2",
                                   scen, "/m2", scen, "yy", year, ".dlem"))
  
  m3_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("D:/biomefluxes/m3",
                                   scen, "/m3", scen, "yy", year, ".dlem"))
  
  pft_ipft <- 91
  npp_pft_ipft <- paste0("g_npp_", pft_ipft)
      
      ## m0 calc
  if(npp_pft_ipft %in% names(m0_s)) {
    npp_r_m0 <- m0_s[[npp_pft_ipft]]
  } else {
    npp_r_m0 <- landvalue_r
  }
      
      ## m1 calc
  if(npp_pft_ipft %in% names(m1_s)) {
    npp_r_m1 <- m1_s[[npp_pft_ipft]]
  } else {
    npp_r_m1 <- landvalue_r
  }
      
      ## m2 calc
  if(npp_pft_ipft %in% names(m2_s)) {
    npp_r_m2 <- m2_s[[npp_pft_ipft]]
  } else {
    npp_r_m2 <- landvalue_r
  }
      
      ## m3 calc
  if(npp_pft_ipft %in% names(m3_s)) {
    npp_r_m3 <- m3_s[[npp_pft_ipft]]
  } else {
    npp_r_m3 <- landvalue_r
  }
      
  npp_r <- npp_r_m0 + npp_r_m1 + npp_r_m2 + npp_r_m3

  ## seperate crop and pastures
  writeRaster(npp_r, paste0("npp_trendy_", year, ".tif"), format="GTiff", overwrite=TRUE)
  
  return(0)
}


## create netcdf
scen <- "s3"

cl <- makeCluster(16)
clusterExport(cl, c("landvalue_r", "%>%", "readDlemBin",
                      "readDlemBiomeflux", "landmask_r"))
biomepft_res <- parLapply(cl, 1900:2017, biomeReadBin, scen=scen)
stopCluster(cl)

 