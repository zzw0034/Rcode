## Title: dlem2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-21


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("D:/hao_shi/gcp18/trendyv7/result")

## biome types (build a new type of managed pasture)
pft_types <- c(0, 91, 92, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 21, 22, 23, 61, 65, 66, 67, 68, 
               69, 70, 71, 72, 73, 74)
pft_names <- c("pft1=Urban", 'pft2=Crops', "pft3=Managed pasture", 'pft4=Tundra',
               'pft5=Boreal needle evergreen forest', 'pft6=Boreal needle deciduous forest',
               'pft7=Temperate broadleaf deciduous forest', 'pft8=Temperate broadleaf evergreen forest',
               'pft9=Temperate needle evergreen forest', 'pft10=Temperate needle deciduous forest', 
               'pft11=Tropical broadleaf deciduous forest', 'pft12=Tropical broadleaf evergreen forest',
               'pft13=Deciduous shrub', 'pft14=Evergreen shrub', 'pft15=C3 grass', 'pft16=C4 grass',
               'pft17=Herbaceous wetland', 'pft18=Boreal woody wetland', 'Temporate woody wetland', 'Tropical woody wetland',
               'South tundra', 'South Temperate broadleaf deciduous forest', 
               'South temperate broad evergreen forest', 'South temperate needle evergreen forest',
               'South temperate needle deciduous forest', 'South tropical broadleaf deciduous forest',
               'South tropical broadleaf evergreen forest', 'South deciduous shrub', 
               'South evergreen shrub', 'South C3 grass', 'South C4 grass')

## read in non-veg luc
## convert to binary format
landmask <- new("GDALReadOnlyDataset", 
                "S:/hao_shi/gcp18/trendyv7/mask/mask_new/w001001.adf")
landmask_r <- asSGDF_GROD(landmask)
gridded(landmask_r) <- TRUE
landmask_r <- raster(landmask_r)
proj4string(landmask_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

## read in managed pastures
tmp_s <- stack("S:/hao_shi/gcp18/data/luh2/states.nc", varname="pastr", bands=851:1168)
cl <- makeCluster(5)
clusterExport(cl, c("landmask_r", "%>%", "tmp_s"))
tmp_s_proj <- parLapply(cl, 1:nlayers(tmp_s), fun=function(i) {
  library(raster)
  tmp_r <- tmp_s[[i]] %>% projectRaster(., landmask_r)
  return(tmp_r)
})
stopCluster(cl)
pastr_s <- do.call(stack, tmp_s_proj)
rm(list=c("tmp_s", "tmp_s_proj"))

## read in water body and bareland
fbare_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_bare.bin")
fglacier_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_glacier.bin")
flake_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_lake.bin")
focean_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_ocean.bin")
friver_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_river.bin")

fland_r <- 1 - (fbare_r + fglacier_r + flake_r + focean_r + friver_r)
fland_r[fland_r <= 0] <- 0
fland_r <- fland_r %>% mask(., landmask_r)

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

pft <- 1:length(pft_types)
pft_dim <- ncdim_def("pft", "PFT", as.integer(pft))

time <- 365 * (0:317)
time_dim <- ncdim_def("time", "days since 1700-01-01 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

lcReadBin <- function(itime) {
  library(raster)
  var_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  
  fcro_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fcro/cro_", itime + 1700 - 1, ".bin"))
  fimp_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fimp/imp_", itime + 1700 - 1, ".bin"))
  fpas_r <- pastr_s[[itime]]
  fvegland_r <- fland_r - fimp_r 
  
  fpft1_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft1/fpft1_", itime + 1700 - 1, ".bin"))
  fpft2_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft2/fpft2_", itime + 1700 - 1, ".bin"))
  fpft3_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft3/fpft3_", itime + 1700 - 1, ".bin"))
  fpft4_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft4/fpft4_", itime + 1700 - 1, ".bin"))
  
  tpft1_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft1/tpft1_", itime + 1700 - 1, ".bin"))
  tpft2_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft2/tpft2_", itime + 1700 - 1, ".bin"))
  tpft3_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft3/tpft3_", itime + 1700 - 1, ".bin"))
  tpft4_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft4/tpft4_", itime + 1700 - 1, ".bin"))
  
  for(ipft in pft) {
    pft_ipft <- pft_types[ipft]
    
    if(pft_ipft==0) {
      var_arr[, , 1] <- values(fimp_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
    } else if(pft_ipft==91) {
      var_arr[, , 2] <- ((fcro_r * fvegland_r - fpas_r)* (landmask_r / landmask_r)) %>% values  %>% 
        c(rep(NA, 3*720), ., rep(NA, 3*720))
    } else if(pft_ipft==92) {
      var_arr[, , 3] <- values(fpas_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
    } else {
      tpft1_r_pft <- tpft1_r
      tpft1_r_pft[tpft1_r_pft != pft_ipft] <- 0
      tpft1_r_pft[tpft1_r_pft != 0] <- 1
      
      tpft2_r_pft <- tpft2_r
      tpft2_r_pft[tpft2_r_pft != pft_ipft] <- 0
      tpft2_r_pft[tpft2_r_pft != 0] <- 1
      
      tpft3_r_pft <- tpft3_r
      tpft3_r_pft[tpft3_r_pft != pft_ipft] <- 0
      tpft3_r_pft[tpft3_r_pft != 0] <- 1
      
      tpft4_r_pft <- tpft4_r
      tpft4_r_pft[tpft4_r_pft != pft_ipft] <- 0
      tpft4_r_pft[tpft4_r_pft != 0] <- 1
      
      fpft_r <- (fpft1_r*tpft1_r_pft + fpft2_r*tpft2_r_pft + fpft3_r*tpft3_r_pft + fpft4_r*tpft4_r_pft) * fvegland_r
      var_arr[, , ipft] <- values(fpft_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
    }
  }
  return(var_arr)
}


##define variables
var_name <- c("landCoverFrac")
var_unit <- "% of each grid cell"
var_lname <- "land cover fraction"

var_def <- ncvar_def(var_name, var_unit, list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                     var_lname, "float", compression=9)

## create netcdf
# nc_file <- paste0("DLEM_S3_landCoverFrac_", sprintf("%03d", itime), ".nc")
nc_file <- "DLEM_S3_landCoverFrac.nc"
nc <- nc_create(nc_file, list(var_def), force_v4=TRUE)

chunks <- (1:318) %/% 17
for(i in unique(chunks)) {
  ind <- which(chunks==i)
  
  ## parallel generation
  cl <- makeCluster(16)
  clusterExport(cl, c("landmask_r", "%>%", "lon", "lat", "pft", "readDlemBin",
                      "fland_r", "pastr_s", "pft_types"))
  lcpft_res <- parLapply(cl, ind, lcReadBin)
  stopCluster(cl)
  
  cat("chunk: ", i, "\n")
  ## put the variables
  lcWrite_res <- lapply(ind, FUN=function(itime) {
    ncvar_put(nc, var_def, lcpft_res[[which(ind==itime)]], start=c(1, 1, 1, itime), 
              count=c(-1, -1, -1, 1))
    cat("itime: ", itime, "\n")
  })
}


## put additional attributes
ncatt_put(nc, "lon", "axis", "X")
ncatt_put(nc, "lon", "long_name", "longitude")

ncatt_put(nc, "lat", "axis", "Y")
ncatt_put(nc, "lat", "long_name", "latitude")

## add global atts
ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
ncatt_put(nc, 0, "pfts", paste(pft_names, collapse=", "))
ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")

## close netcdf
nc_close(nc)
rm(lcpft_res)


## create s0 land frac
##define variables
var_name <- c("landCoverFrac")
var_unit <- "% of each grid cell"
var_lname <- "land cover fraction"

var_def <- ncvar_def(var_name, var_unit, list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                     var_lname, "float", compression=9)

## create netcdf
# nc_file <- paste0("DLEM_S3_landCoverFrac_", sprintf("%03d", itime), ".nc")
nc_file <- "DLEM_S0_landCoverFrac.nc"
nc <- nc_create(nc_file, list(var_def), force_v4=TRUE)

var_arr_1700 <- lcReadBin(1)

for(itime in 1:318) {
  ncvar_put(nc, var_def, var_arr_1700, start=c(1, 1, 1, itime), 
              count=c(-1, -1, -1, 1))
  cat("itime: ", itime, "\n")
}


## put additional attributes
ncatt_put(nc, "lon", "axis", "X")
ncatt_put(nc, "lon", "long_name", "longitude")

ncatt_put(nc, "lat", "axis", "Y")
ncatt_put(nc, "lat", "long_name", "latitude")

## add global atts
ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
ncatt_put(nc, 0, "pfts", paste(pft_names, collapse=", "))
ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")

## close netcdf
nc_close(nc)
rm(var_arr_1700)
















#########################################
# landCoverPFT <- function(itime, var_name, var_unit, var_lname) {
#   # cat(itime, "\n")
#   library(raster)
#   library(ncdf4)
#   
#   var_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
#   
#   fcro_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fcro/cro_", itime + 1700 - 1, ".bin"))
#   fimp_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fimp/imp_", itime + 1700 - 1, ".bin"))
#   fpas_r <- pastr_s[[itime]]
#   fvegland_r <- fland_r - fimp_r 
#   
#   fpft1_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft1/fpft1_", itime + 1700 - 1, ".bin"))
#   fpft2_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft2/fpft2_", itime + 1700 - 1, ".bin"))
#   fpft3_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft3/fpft3_", itime + 1700 - 1, ".bin"))
#   fpft4_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fpft4/fpft4_", itime + 1700 - 1, ".bin"))
#   
#   tpft1_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft1/tpft1_", itime + 1700 - 1, ".bin"))
#   tpft2_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft2/tpft2_", itime + 1700 - 1, ".bin"))
#   tpft3_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft3/tpft3_", itime + 1700 - 1, ".bin"))
#   tpft4_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/tpft4/tpft4_", itime + 1700 - 1, ".bin"))
#   
#   for(ipft in pft) {
#     pft_ipft <- pft_types[ipft]
#     
#     if(pft_ipft==0) {
#       var_arr[, , 1] <- values(fimp_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
#     } else if(pft_ipft==91) {
#       var_arr[, , 2] <- ((fcro_r * fvegland_r - fpas_r)* (landmask_r / landmask_r)) %>% values  %>% 
#         c(rep(NA, 3*720), ., rep(NA, 3*720))
#     } else if(pft_ipft==92) {
#       var_arr[, , 3] <- values(fpas_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
#     } else {
#       tpft1_r_pft <- tpft1_r
#       tpft1_r_pft[tpft1_r_pft != pft_ipft] <- 0
#       tpft1_r_pft[tpft1_r_pft != 0] <- 1
#       
#       tpft2_r_pft <- tpft2_r
#       tpft2_r_pft[tpft2_r_pft != pft_ipft] <- 0
#       tpft2_r_pft[tpft2_r_pft != 0] <- 1
#       
#       tpft3_r_pft <- tpft3_r
#       tpft3_r_pft[tpft3_r_pft != pft_ipft] <- 0
#       tpft3_r_pft[tpft3_r_pft != 0] <- 1
#       
#       tpft4_r_pft <- tpft4_r
#       tpft4_r_pft[tpft4_r_pft != pft_ipft] <- 0
#       tpft4_r_pft[tpft4_r_pft != 0] <- 1
#       
#       fpft_r <- (fpft1_r*tpft1_r_pft + fpft2_r*tpft2_r_pft + fpft3_r*tpft3_r_pft + fpft4_r*tpft4_r_pft) * fvegland_r
#       var_arr[, , ipft] <- values(fpft_r * (landmask_r / landmask_r)) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
#     }
#     cat(ipft, "\n")
#   }
#   
#   var_def <- ncvar_def(var_name, var_unit, list(lon_dim, lat_dim, pft_dim), fill_value,
#                        var_lname, "float", compression=9)
#   
#   ## create netcdf
#   nc_file <- paste0("DLEM_S3_landCoverFrac_", sprintf("%03d", itime), ".nc")
#   nc <- nc_create(nc_file, list(var_def), force_v4=TRUE)
#   
#   ## put the variables
#   ncvar_put(nc, var_def, var_arr)
#   
#   ## put additional attributes
#   ncatt_put(nc, "lon", "axis", "X")
#   ncatt_put(nc, "lon", "long_name", "longitude")
#   
#   ncatt_put(nc, "lat", "axis", "Y")
#   ncatt_put(nc, "lat", "long_name", "latitude")
#   
#   ## add global atts
#   ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
#   ncatt_put(nc, 0, "pfts", paste(pft_names, collapse=", "))
#   ncatt_put(nc, 0, "contact", "Hanqin Tian, tianh@auburn.edu; Hao Shi, hzs0087@auburn.edu")
#   
#   ## close netcdf
#   nc_close(nc)
#   
#   return(1)
# }
# 
# ##define variables
# lcpft_name <- c("landCoverFrac")
# lcpft_unit <- "% of each grid cell"
# lcpft_lname <- "land cover fraction"
# 
# ## parallel generation
# cl <- makeCluster(6)
# clusterExport(cl, c("landmask_r", "%>%", "lon", "lat", "pft", "readDlemBin", "fill_value",
#                     "fland_r", "lon_dim", "lat_dim", "pft_dim", "pastr_s", "pft_types",
#                     "pft_names"))
# lcpft_res <- parLapply(cl, 1:318, landCoverPFT, var_name=lcpft_name, var_unit=lcpft_unit,
#                        var_lname=lcpft_lname)
# stopCluster(cl)



