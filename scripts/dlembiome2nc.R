## Title: dlembiome2nc
## Author: hzs0087@auburn.edu 
## Date: 2018-08-22


## import packages
pacs <- c("ncdf4", "raster", "rgdal", "magrittr", "parallel")
imp <- lapply(pacs, require, character.only=TRUE)
source("S:/hao_shi/gcp18/scripts/read_dlem.R")

## set workspace
ori_dir <- setwd("D:/hao_shi/gcp18/trendyv7/result")

## biome types (build a new type of managed pasture)
pft_types <- c(0, 91, 92, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 21, 22, 23, 61, 65, 66, 67, 68, 
               69, 70, 71, 72, 73, 74)
pft_names <- c("Urban", 'Crops', "Managed pasture", 'Tundra', 'Boreal needle evergreen forest', 
               'Boreal needle deciduous forest',
               'Temperate broadleaf deciduous forest', 'Temperate broadleaf evergreen forest',
               'Temperate needle evergreen forest', 'Temperate needle deciduous forest', 
               'Tropical broadleaf deciduous forest', 'Tropical broadleaf evergreen forest',
               'Deciduous shrub', 'Evergreen shrub', 'C3 grass', 'C4 grass', 'Herbaceous wetland',
               'Boreal woody wetland', 'Temporate woody wetland', 'Tropical woody wetland',
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

landvalue_r <- landmask_r - 1

## read in managed pastures
tmp_s <- stack("S:/hao_shi/gcp18/data/luh2/states.nc", varname="pastr", bands=1051:1168)
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

## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

pft <- 1:length(pft_types)
pft_dim <- ncdim_def("pft", "PFT", as.integer(pft))

time <- 1:(118*12) - 1
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

biomeReadBin <- function(itime, scen) { ## itime starts from 1 to 1416
  library(raster)
  gpp_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  npp_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  nep_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  trans_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  et_arr <- array(numeric(), dim=c(length(lon), length(lat), length(pft)))
  
  itime <- itime - 1
  mon <- itime %% 12
  year <- itime %/% 12 + 1900
  
  m0_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("D:/hao_shi/gcp18/trendyv7/output/biomefluxes/m0",
                                   scen, "/m0", scen, "yy", year, "m", mon, ".dlem"))
  
  m1_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("D:/hao_shi/gcp18/trendyv7/output/biomefluxes/m1",
                                   scen, "/m1", scen, "yy", year, "m", mon, ".dlem"))
  
  m2_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin", 
                            paste0("D:/hao_shi/gcp18/trendyv7/output/biomefluxes/m2",
                                   scen, "/m2", scen, "yy", year, "m", mon, ".dlem"))
  
  m3_s <- readDlemBiomeflux("S:/hao_shi/gcp18/trendyv7/mask/maskidx_new.bin",
                            paste0("D:/hao_shi/gcp18/trendyv7/output/biomefluxes/m3",
                                   scen, "/m3", scen, "yy", year, "m", mon, ".dlem"))
  
  fcro_r <- readDlemBin(paste0("S:/hao_shi/gcp18/trendyv7/cohort/fcro/cro_", year, ".bin"))
  fcro_r[fcro_r > 0] <- 1
  fcro_r <- fcro_r * (landmask_r / landmask_r)
  
  fpastr_r <- pastr_s[[year - 1900 + 1]]
  fpastr_r[fpastr_r > 0] <- 1
  
  for(ipft in pft) {
    pft_ipft <- pft_types[ipft]
    
    if(pft_ipft==0) {
      gpp_arr[, , 1] <- values(landvalue_r) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      npp_arr[, , 1] <- values(landvalue_r) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      nep_arr[, , 1] <- values(landvalue_r) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      trans_arr[, , 1] <- values(landvalue_r) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      et_arr[, , 1] <- values(landvalue_r) %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      
    } else if(pft_ipft==91) {
      gpp_pft_ipft <- paste0("g_gpp_", pft_ipft)
      npp_pft_ipft <- paste0("g_npp_", pft_ipft)
      g_mr_pft_ipft <- paste0("g_mr_", pft_ipft)
      g_gr_pft_ipft <- paste0("g_gr_", pft_ipft)
      trans_pft_ipft <- paste0("g_trans_", pft_ipft)
      evap_pft_ipft <- paste0("g_evap_", pft_ipft)
      
      ## m0 calc
      if(gpp_pft_ipft %in% names(m0_s)) {
        gpp_r_m0 <- m0_s[[gpp_pft_ipft]]
        npp_r_m0 <- m0_s[[npp_pft_ipft]]
        nep_r_m0 <- npp_r_m0 - m0_s[[g_mr_pft_ipft]] - m0_s[[g_gr_pft_ipft]]
        trans_r_m0 <- m0_s[[trans_pft_ipft]]
        et_r_m0 <- trans_r_m0 + m0_s[[evap_pft_ipft]]
      } else {
        gpp_r_m0 <- landvalue_r
        npp_r_m0 <- landvalue_r
        nep_r_m0 <- landvalue_r
        trans_r_m0 <- landvalue_r
        et_r_m0 <- landvalue_r
      }
      
      ## m1 calc
      if(gpp_pft_ipft %in% names(m1_s)) {
        gpp_r_m1 <- m1_s[[gpp_pft_ipft]]
        npp_r_m1 <- m1_s[[npp_pft_ipft]]
        nep_r_m1 <- npp_r_m1 - m1_s[[g_mr_pft_ipft]] - m1_s[[g_gr_pft_ipft]]
        trans_r_m1 <- m1_s[[trans_pft_ipft]]
        et_r_m1 <- trans_r_m1 + m1_s[[evap_pft_ipft]]
      } else {
        gpp_r_m1 <- landvalue_r
        npp_r_m1 <- landvalue_r
        nep_r_m1 <- landvalue_r
        trans_r_m1 <- landvalue_r
        et_r_m1 <- landvalue_r
      }
      
      ## m2 calc
      if(gpp_pft_ipft %in% names(m2_s)) {
        gpp_r_m2 <- m2_s[[gpp_pft_ipft]]
        npp_r_m2 <- m2_s[[npp_pft_ipft]]
        nep_r_m2 <- npp_r_m2 - m2_s[[g_mr_pft_ipft]] - m2_s[[g_gr_pft_ipft]]
        trans_r_m2 <- m2_s[[trans_pft_ipft]]
        et_r_m2 <- trans_r_m2 + m2_s[[evap_pft_ipft]]
      } else {
        gpp_r_m2 <- landvalue_r
        npp_r_m2 <- landvalue_r
        nep_r_m2 <- landvalue_r
        trans_r_m2 <- landvalue_r
        et_r_m2 <- landvalue_r
      }
      
      ## m3 calc
      if(gpp_pft_ipft %in% names(m3_s)) {
        gpp_r_m3 <- m3_s[[gpp_pft_ipft]]
        npp_r_m3 <- m3_s[[npp_pft_ipft]]
        nep_r_m3 <- npp_r_m3 - m3_s[[g_mr_pft_ipft]] - m3_s[[g_gr_pft_ipft]]
        trans_r_m3 <- m3_s[[trans_pft_ipft]]
        et_r_m3 <- trans_r_m3 + m3_s[[evap_pft_ipft]]
      } else {
        gpp_r_m3 <- landvalue_r
        npp_r_m3 <- landvalue_r
        nep_r_m3 <- landvalue_r
        trans_r_m3 <- landvalue_r
        et_r_m3 <- landvalue_r
      }
      
      ## combine all masks
      gpp_r <- gpp_r_m0 + gpp_r_m1 + gpp_r_m2 + gpp_r_m3
      npp_r <- npp_r_m0 + npp_r_m1 + npp_r_m2 + npp_r_m3
      nep_r <- nep_r_m0 + nep_r_m1 + nep_r_m2 + nep_r_m3
      trans_r <- trans_r_m0 + trans_r_m1 + trans_r_m2 + trans_r_m3
      et_r <- et_r_m0 + et_r_m1 + et_r_m2 + et_r_m3
      
      ## seperate crop and pastures
      gpp_r <- gpp_r * fcro_r; npp_r <- npp_r * fcro_r; nep_r <- nep_r * fcro_r;
      trans_r <- trans_r * fcro_r; et_r <- et_r * fcro_r;
      
      ## write to array
      gpp_arr[, , 2] <- gpp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      npp_arr[, , 2] <- npp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      nep_arr[, , 2] <- nep_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      trans_arr[, , 2] <- trans_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      et_arr[, , 2] <- et_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      
      
    } else if(pft_ipft==92) {
      
      gpp_pft_ipft <- paste0("g_gpp_", 91)
      npp_pft_ipft <- paste0("g_npp_", 91)
      g_mr_pft_ipft <- paste0("g_mr_", 91)
      g_gr_pft_ipft <- paste0("g_gr_", 91)
      trans_pft_ipft <- paste0("g_trans_", 91)
      evap_pft_ipft <- paste0("g_evap_", 91)
      
      ## m0 calc
      if(gpp_pft_ipft %in% names(m0_s)) {
        gpp_r_m0 <- m0_s[[gpp_pft_ipft]]
        npp_r_m0 <- m0_s[[npp_pft_ipft]]
        nep_r_m0 <- npp_r_m0 - m0_s[[g_mr_pft_ipft]] - m0_s[[g_gr_pft_ipft]]
        trans_r_m0 <- m0_s[[trans_pft_ipft]]
        et_r_m0 <- trans_r_m0 + m0_s[[evap_pft_ipft]]
      } else {
        gpp_r_m0 <- landvalue_r
        npp_r_m0 <- landvalue_r
        nep_r_m0 <- landvalue_r
        trans_r_m0 <- landvalue_r
        et_r_m0 <- landvalue_r
      }
      
      ## m1 calc
      if(gpp_pft_ipft %in% names(m1_s)) {
        gpp_r_m1 <- m1_s[[gpp_pft_ipft]]
        npp_r_m1 <- m1_s[[npp_pft_ipft]]
        nep_r_m1 <- npp_r_m1 - m1_s[[g_mr_pft_ipft]] - m1_s[[g_gr_pft_ipft]]
        trans_r_m1 <- m1_s[[trans_pft_ipft]]
        et_r_m1 <- trans_r_m1 + m1_s[[evap_pft_ipft]]
      } else {
        gpp_r_m1 <- landvalue_r
        npp_r_m1 <- landvalue_r
        nep_r_m1 <- landvalue_r
        trans_r_m1 <- landvalue_r
        et_r_m1 <- landvalue_r
      }
      
      ## m2 calc
      if(gpp_pft_ipft %in% names(m2_s)) {
        gpp_r_m2 <- m2_s[[gpp_pft_ipft]]
        npp_r_m2 <- m2_s[[npp_pft_ipft]]
        nep_r_m2 <- npp_r_m2 - m2_s[[g_mr_pft_ipft]] - m2_s[[g_gr_pft_ipft]]
        trans_r_m2 <- m2_s[[trans_pft_ipft]]
        et_r_m2 <- trans_r_m2 + m2_s[[evap_pft_ipft]]
      } else {
        gpp_r_m2 <- landvalue_r
        npp_r_m2 <- landvalue_r
        nep_r_m2 <- landvalue_r
        trans_r_m2 <- landvalue_r
        et_r_m2 <- landvalue_r
      }
      
      ## m3 calc
      if(gpp_pft_ipft %in% names(m3_s)) {
        gpp_r_m3 <- m3_s[[gpp_pft_ipft]]
        npp_r_m3 <- m3_s[[npp_pft_ipft]]
        nep_r_m3 <- npp_r_m3 - m3_s[[g_mr_pft_ipft]] - m3_s[[g_gr_pft_ipft]]
        trans_r_m3 <- m3_s[[trans_pft_ipft]]
        et_r_m3 <- trans_r_m3 + m3_s[[evap_pft_ipft]]
      } else {
        gpp_r_m3 <- landvalue_r
        npp_r_m3 <- landvalue_r
        nep_r_m3 <- landvalue_r
        trans_r_m3 <- landvalue_r
        et_r_m3 <- landvalue_r
      }
      
      ## combine all masks
      gpp_r <- gpp_r_m0 + gpp_r_m1 + gpp_r_m2 + gpp_r_m3
      npp_r <- npp_r_m0 + npp_r_m1 + npp_r_m2 + npp_r_m3
      nep_r <- nep_r_m0 + nep_r_m1 + nep_r_m2 + nep_r_m3
      trans_r <- trans_r_m0 + trans_r_m1 + trans_r_m2 + trans_r_m3
      et_r <- et_r_m0 + et_r_m1 + et_r_m2 + et_r_m3
      
      ## seperate crop and pastures
      gpp_r <- gpp_r * fpastr_r; npp_r <- npp_r * fpastr_r; nep_r <- nep_r * fpastr_r;
      trans_r <- trans_r * fpastr_r; et_r <- et_r * fpastr_r;
      
      ## write to array
      gpp_arr[, , 3] <- gpp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      npp_arr[, , 3] <- npp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      nep_arr[, , 3] <- nep_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      trans_arr[, , 3] <- trans_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      et_arr[, , 3] <- et_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      
    } else {
      gpp_pft_ipft <- paste0("g_gpp_", pft_ipft)
      npp_pft_ipft <- paste0("g_npp_", pft_ipft)
      g_mr_pft_ipft <- paste0("g_mr_", pft_ipft)
      g_gr_pft_ipft <- paste0("g_gr_", pft_ipft)
      trans_pft_ipft <- paste0("g_trans_", pft_ipft)
      evap_pft_ipft <- paste0("g_evap_", pft_ipft)
      
      ## m0 calc
      if(gpp_pft_ipft %in% names(m0_s)) {
        gpp_r_m0 <- m0_s[[gpp_pft_ipft]]
        npp_r_m0 <- m0_s[[npp_pft_ipft]]
        nep_r_m0 <- npp_r_m0 - m0_s[[g_mr_pft_ipft]] - m0_s[[g_gr_pft_ipft]]
        trans_r_m0 <- m0_s[[trans_pft_ipft]]
        et_r_m0 <- trans_r_m0 + m0_s[[evap_pft_ipft]]
      } else {
        gpp_r_m0 <- landvalue_r
        npp_r_m0 <- landvalue_r
        nep_r_m0 <- landvalue_r
        trans_r_m0 <- landvalue_r
        et_r_m0 <- landvalue_r
      }
      
      ## m1 calc
      if(gpp_pft_ipft %in% names(m1_s)) {
        gpp_r_m1 <- m1_s[[gpp_pft_ipft]]
        npp_r_m1 <- m1_s[[npp_pft_ipft]]
        nep_r_m1 <- npp_r_m1 - m1_s[[g_mr_pft_ipft]] - m1_s[[g_gr_pft_ipft]]
        trans_r_m1 <- m1_s[[trans_pft_ipft]]
        et_r_m1 <- trans_r_m1 + m1_s[[evap_pft_ipft]]
      } else {
        gpp_r_m1 <- landvalue_r
        npp_r_m1 <- landvalue_r
        nep_r_m1 <- landvalue_r
        trans_r_m1 <- landvalue_r
        et_r_m1 <- landvalue_r
      }
      
      ## m2 calc
      if(gpp_pft_ipft %in% names(m2_s)) {
        gpp_r_m2 <- m2_s[[gpp_pft_ipft]]
        npp_r_m2 <- m2_s[[npp_pft_ipft]]
        nep_r_m2 <- npp_r_m2 - m2_s[[g_mr_pft_ipft]] - m2_s[[g_gr_pft_ipft]]
        trans_r_m2 <- m2_s[[trans_pft_ipft]]
        et_r_m2 <- trans_r_m2 + m2_s[[evap_pft_ipft]]
      } else {
        gpp_r_m2 <- landvalue_r
        npp_r_m2 <- landvalue_r
        nep_r_m2 <- landvalue_r
        trans_r_m2 <- landvalue_r
        et_r_m2 <- landvalue_r
      }
      
      ## m3 calc
      if(gpp_pft_ipft %in% names(m3_s)) {
        gpp_r_m3 <- m3_s[[gpp_pft_ipft]]
        npp_r_m3 <- m3_s[[npp_pft_ipft]]
        nep_r_m3 <- npp_r_m3 - m3_s[[g_mr_pft_ipft]] - m3_s[[g_gr_pft_ipft]]
        trans_r_m3 <- m3_s[[trans_pft_ipft]]
        et_r_m3 <- trans_r_m3 + m3_s[[evap_pft_ipft]]
      } else {
        gpp_r_m3 <- landvalue_r
        npp_r_m3 <- landvalue_r
        nep_r_m3 <- landvalue_r
        trans_r_m3 <- landvalue_r
        et_r_m3 <- landvalue_r
      }
      
      ## combine all masks
      gpp_r <- gpp_r_m0 + gpp_r_m1 + gpp_r_m2 + gpp_r_m3
      npp_r <- npp_r_m0 + npp_r_m1 + npp_r_m2 + npp_r_m3
      nep_r <- nep_r_m0 + nep_r_m1 + nep_r_m2 + nep_r_m3
      trans_r <- trans_r_m0 + trans_r_m1 + trans_r_m2 + trans_r_m3
      et_r <- et_r_m0 + et_r_m1 + et_r_m2 + et_r_m3
      
      ## write to array
      gpp_arr[, , ipft] <- gpp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      npp_arr[, , ipft] <- npp_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      nep_arr[, , ipft] <- nep_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      trans_arr[, , ipft] <- trans_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
      et_arr[, , ipft] <- et_r %>% values %>% c(rep(NA, 3*720), ., rep(NA, 3*720))
    }
  }
  return(list(gpp_arr, npp_arr, nep_arr, trans_arr, et_arr))
}


## create netcdf
scen <- "s1"

# chunks <- (1:1416) %/% 16
# for (i in unique(chunks)) {
#     ind <- which(chunks==i)
# 
#     cl <- makeCluster(16)
#     clusterExport(cl, c("landvalue_r", "%>%", "lon", "lat", "pft", "readDlemBin",
#                         "readDlemBiomeflux", "pastr_s", "pft_types", "landmask_r"))
#     biomepft_res <- parLapply(cl, ind, biomeReadBin, scen=scen)
#     stopCluster(cl)
#     
#     save(biomepft_res, file=paste0("biomepft_",scen, "_", i, ".RData"))
# 
#     cat("chunk: ", i, "\n")
# }

##define variables
gpp_def <- ncvar_def("gpppft", "kg C m-2 s-1", list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                     "Vegtype level GPP", "float", compression=9)

npp_def <- ncvar_def("npppft", "kg C m-2 s-1", list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                     "Vegtype level NPP", "float", compression=9)

nep_def <- ncvar_def("neppft", "kg C m-2 s-1", list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                     "Vegtype level NEP", "float", compression=9)

trans_def <- ncvar_def("tanspft", "W m-2", list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                       "Vegtype level transpiration", "float", compression=9)

et_def <- ncvar_def("evapotranspft", "W m-2", list(lon_dim, lat_dim, pft_dim, time_dim), fill_value,
                    "Vegtype level evapotranspiration", "float", compression=9)

nc_gpp <- nc_create(paste0("DLEM_",toupper(scen), "_gpppft.nc"), list(gpp_def), force_v4=TRUE)
nc_npp <- nc_create(paste0("DLEM_",toupper(scen), "_npppft.nc"), list(npp_def), force_v4=TRUE)
nc_nep <- nc_create(paste0("DLEM_",toupper(scen), "_neppft.nc"), list(nep_def), force_v4=TRUE)
nc_trans <- nc_create(paste0("DLEM_",toupper(scen), "_transpft.nc"), list(trans_def), force_v4=TRUE)
nc_et <- nc_create(paste0("DLEM_",toupper(scen), "_evapotranspft.nc"), list(et_def), force_v4=TRUE)

## parallel generation
chunks <- (1:1416) %/% 16
days_mon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

for (i in unique(chunks)) {
  ind <- which(chunks==i)

  cl <- makeCluster(16)
  clusterExport(cl, c("landvalue_r", "%>%", "lon", "lat", "pft", "readDlemBin",
                      "readDlemBiomeflux", "pastr_s", "pft_types", "landmask_r"))
  biomepft_res <- parLapply(cl, ind, biomeReadBin, scen=scen)
  stopCluster(cl)

  cat("chunk: ", i, "\n")

  ## put the variables
  ## put the variables
  biomeWrite_res <- lapply(ind, FUN=function(itime) {
    mon <- (itime - 1) %% 12 + 1
    secs <- days_mon[mon] * 24 * 3600

    ncvar_put(nc_gpp, gpp_def, biomepft_res[[which(ind==itime)]][[1]] * 10^-3 / secs,
              start=c(1, 1, 1, itime), count=c(-1, -1, -1, 1))

    ncvar_put(nc_npp, npp_def, biomepft_res[[which(ind==itime)]][[2]] * 10^-3 / secs,
              start=c(1, 1, 1, itime), count=c(-1, -1, -1, 1))

    ncvar_put(nc_nep, nep_def, biomepft_res[[which(ind==itime)]][[3]] * 10^-3 / secs,
              start=c(1, 1, 1, itime), count=c(-1, -1, -1, 1))

    ncvar_put(nc_trans, trans_def, biomepft_res[[which(ind==itime)]][[4]] * 2.45 * 10^6 / secs,
              start=c(1, 1, 1, itime), count=c(-1, -1, -1, 1))

    ncvar_put(nc_et, et_def, biomepft_res[[which(ind==itime)]][[5]] * 2.45 * 10^6 / secs,
              start=c(1, 1, 1, itime), count=c(-1, -1, -1, 1))
    cat("itime: ", itime, "\n")
  })
}

## put additional attributes
for(nc_var in c("nc_gpp", "nc_npp", "nc_nep", "nc_trans", "nc_et")) {
  nc <- eval(parse(text=nc_var))
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
}

rm(biomepft_res)





