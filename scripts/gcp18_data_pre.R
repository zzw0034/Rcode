## Title: gcp18_data_prepare
## Author: hzs0087@auburn.edu 
## Date: 2018-07-02

library(raster)
library(rgdal)
library(magrittr)
library(parallel)
library(R.utils)
library(bigmemory)

ori_dir <- setwd("S:/hao_shi/gcp18/trendyv7/clim/")

## custom function
readGz <- function(in_file) {
  in_nc <- gunzip(in_file, temporary=TRUE, remove=FALSE)
  in_6h <- brick(in_nc, values=TRUE)
  
  ind <- (1:nlayers(in_6h) - 1) %/% 4
  ndays <- ind %>% unique %>% length
  
  var <- strsplit(in_file, '\\.') %>% sapply(., "[[", 5)
  
  env <- environment()
  cl <- makeCluster(16)
  clusterExport(cl, c("in_6h", "ind", "%>%", "var"), envir=env)
  
  in_daily <- parLapply(cl=cl, 1:ndays, fun=function(day) {
    library(raster)
    
    layers <- which(ind == day - 1)
    if(var == "dswrf") {
      in_6h_daily <- in_6h[[layers]] %>% sum
      in_6h_daily <- in_6h_daily / (3600 * 24)
    } else if(var == "pre") {
      in_6h_daily <- in_6h[[layers]] %>% sum
    } else if(var =="tmax") {
      in_6h_daily <- in_6h[[layers]] %>% max
      in_6h_daily <- in_6h_daily - 273.15
    } else if(var == "tmin") {
      in_6h_daily <- in_6h[[layers]] %>% min
      in_6h_daily <- in_6h_daily - 273.15
    } else if(var == "tmp") {
      in_6h_daily <- in_6h[[layers]] %>% mean
      in_6h_daily <- in_6h_daily - 273.15
    }
    
    return(in_6h_daily)
  })
  
  stopCluster(cl)
  file.remove(in_nc)
  
  in_daily <- do.call("brick", in_daily)
  return(in_daily)
}  

## list all gz files
in_files <- list.files(pattern="gz$", recursive=TRUE)
for(in_file in in_files) {
  in_daily <- readGz(in_file) 
  writeRaster(in_daily, sub("\\.gz", "", in_file), format="CDF", overwrite=TRUE)
}  

## convert to binary format
landmask <- new("GDALReadOnlyDataset", 
                "S:/hao_shi/gcp18/trendyv7/mask/mask_new/w001001.adf")
landmask_r <- asSGDF_GROD(landmask)
gridded(landmask_r) <- TRUE
landmask_r <- raster(landmask_r)
proj4string(landmask_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

nc_files <- list.files(pattern="nc$", recursive=TRUE)
cl <- makeCluster(16)
clusterExport(cl, c("%>%", "landmask_r"))
res <- parLapply(cl, nc_files, fun=function(nc_file) {
  library(raster)
  b <- brick(nc_file) %>% crop(., landmask_r)
  b_vec <- values(b) %>% as.numeric
  b_vec[is.na(b_vec)] <- -9999
  
  var <- strsplit(nc_file, '\\.') %>% sapply(., "[[", 5)
  year <- strsplit(nc_file, '\\.') %>% sapply(., "[[", 6)
  
  bin_con <- file(paste0(var, year, ".bin") %>% file.path(var, .), "wb")
  writeBin(b_vec, bin_con, size=4)
  close(bin_con)
})
stopCluster(cl)

## potential vegetation
for(var in paste0("fveg", 1:4)) {
  var_obj <- new("GDALReadOnlyDataset",
                 paste0("S:/ISIMIP/input/cohort_new/potveg/", var, "/w001001.adf"))
  # getDriver(var_obj)
  # getDriverLongName(getDriver(var_obj))
  var_r <- asSGDF_GROD(var_obj)
  gridded(var_r) <- TRUE
  var_r <- raster(var_r)
  proj4string(var_r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  assign(paste0(var, "_r"), var_r)
}
rm(list=c("var_r", "var_obj"))

## luh2 crop fraction
ori_dir <- setwd("S:/hao_shi/gcp18/trendyv7/")

for(var in c("c3ann", "c3nfx", "c3per", "c4ann", "c4per", "pastr", "urban")) {
  tmp_s <- stack("../data/luh2/states.nc", varname=var)
  cl <- makeCluster(32)
  clusterExport(cl, c("landmask_r", "%>%", "tmp_s"))
  tmp_s_proj <- parLapply(cl, 1:nlayers(tmp_s), fun=function(i) {
    library(raster)
    tmp_r <- tmp_s[[i]] %>% projectRaster(., landmask_r)
    return(tmp_r)
  })
  stopCluster(cl)
  tmp_s_proj <- do.call(stack, tmp_s_proj)
  assign(paste0(var, "_s"), tmp_s_proj)
}
rm(list=c("tmp_s_proj", "tmp_s"))


## land cover data
readDlemBin <- function(bin_file) {
  bin_con <- file(bin_file, "rb")
  data_nu <- readBin(bin_con, "double", n=354 * 720, size=4)
  close(bin_con)
  
  r <- matrix(data_nu, nrow=354, byrow=TRUE) %>% 
    raster(., xmn=-180, xmx=180, ymn=-88.5, ymx=88.5, crs="+proj=longlat +datum=WGS84")
  return(r)
}

fbare_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_bare.bin")
fglacier_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_glacier.bin")
flake_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_lake.bin")
focean_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_ocean.bin")
friver_r <- readDlemBin("S:/ISIMIP/input/cohort/DLEM_cohort/glo_river.bin")

fvegland_r <- 1 - (fbare_r + fglacier_r + flake_r + focean_r + friver_r)
fvegland_r[fvegland_r <= 0] <- 0
fvegland_r <- fvegland_r %>% mask(., landmask_r)

chunks <- 1:nlayers(c3ann_s) %/% 48
nrow <- nrow(c3ann_s)
proj <- proj4string(c3ann_s)
ext <- extent(c3ann_s)
cro_list <- list()

## parallel computation
for(i in unique(chunks)) {
  ind <- which(chunks == i)
  c3ann_mat <- as.big.matrix(values(c3ann_s[[ind]])) 
  c3ann_des <- describe(c3ann_mat)
  
  c3nfx_mat <- as.big.matrix(values(c3nfx_s[[ind]])) 
  c3nfx_des <- describe(c3nfx_mat)
  
  c3per_mat <- as.big.matrix(values(c3per_s[[ind]])) 
  c3per_des <- describe(c3per_mat)
  
  c4ann_mat <- as.big.matrix(values(c4ann_s[[ind]])) 
  c4ann_des <- describe(c4ann_mat)
  
  c4per_mat <- as.big.matrix(values(c4per_s[[ind]])) 
  c4per_des <- describe(c4per_mat)
  
  pastr_mat <- as.big.matrix(values(pastr_s[[ind]])) 
  pastr_des <- describe(pastr_mat)
  
  urban_mat <- as.big.matrix(values(urban_s[[ind]])) 
  urban_des <- describe(urban_mat)
  
  cl <- makeCluster(16)
  clusterExport(cl, c("c3ann_des", "c3nfx_des", "c3per_des", "c4ann_des",
                      "c4per_des", "pastr_des", "nrow", "proj", "ext", "%>%",
                      "fvegland_r", "urban_des"))
  
  cro_list_ind <- parLapply(cl, 1:length(ind), fun=function(j) {
    library(raster)
    library(bigmemory)
    
    c3ann_r <- attach.big.matrix(c3ann_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    c3nfx_r <- attach.big.matrix(c3nfx_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    c3per_r <- attach.big.matrix(c3per_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    c4ann_r <- attach.big.matrix(c4ann_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    c4per_r <- attach.big.matrix(c4per_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    pastr_r <- attach.big.matrix(pastr_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    urban_r <- attach.big.matrix(urban_des) %>% .[, j] %>% 
      matrix(., nrow=nrow, byrow=TRUE) %>% raster(., xmn=ext[1], xmx=ext[2], 
                                                  ymn=ext[3], ymx=ext[4], crs=proj)
    
    cro_r <- c3ann_r + c3nfx_r + c3per_r + c4ann_r + c4per_r + pastr_r
    cro_r <- cro_r / (fvegland_r - urban_r)
    cro_r[cro_r < 0] <- 0
    cro_r[is.na(cro_r)] <- 0
    cro_r[cro_r > 1] <- 1
    
    return(cro_r)
    
  })
  stopCluster(cl)
  cro_list_ind <- do.call(stack, cro_list_ind)
  cro_list <- list(cro_list, cro_list_ind)
  cat(i, "\n")
  
}

cro_s <- do.call(stack, unlist(cro_list))

##write each year into bin
for(i in 851:nlayers(cro_s)) {
  for(var in c("cro", "fveg1", "fveg2", "fveg3", "fveg4")) {
    cro_r <- cro_s[[i]]
    
    if(var=="cro") {
      var_r <- cro_r
    } else {
      var_r <- (1 - cro_r) * eval(parse(text=paste0(var, "_r")))
    }
    
    b_vec <- var_r %>% values %>% as.numeric
    b_vec[is.na(b_vec)] <- 0
    
    var <- sub("veg", "pft", var)
    if(!dir.exists(var)) {
      dir.create(var)
    }
    
    bin_con <- file(paste0(var, "_", 849+i, ".bin") %>% file.path(var, .), "wb") 
    writeBin(b_vec, bin_con, size=4)
    close(bin_con)
    cat(var, "\n")
  }
  cat(i, "\n")
}

## fertilizer
fertl_c3ann_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="fertl_c3ann")
fertl_c3nfx_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="fertl_c3nfx")
fertl_c3per_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="fertl_c3per")
fertl_c4ann_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="fertl_c4ann")
fertl_c4per_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="fertl_c4per")

## irrigation
irrig_c3ann_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="irrig_c3ann")
irrig_c3nfx_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="irrig_c3nfx")
irrig_c3per_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="irrig_c3per")
irrig_c4ann_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="irrig_c4ann")
irrig_c4per_s <- stack("S:/hao_shi/gcp18/data/luh2/management.nc", varname="irrig_c4per")

fertl_res <- lapply(851:nlayers(fertl_c3ann_s), FUN=function(i) {
  c3ann_r <- c3ann_s[[i]]; c3nfx_r <- c3nfx_s[[i]]
  c3per_r <- c3per_s[[i]]; c4ann_r <- c4ann_s[[i]]
  c4per_r <- c4per_s[[i]];
  
  fertl_c3ann_r <- fertl_c3ann_s[[i]] %>% projectRaster(., landmask_r)
  fertl_c3nfx_r <- fertl_c3nfx_s[[i]] %>% projectRaster(., landmask_r)
  fertl_c3per_r <- fertl_c3per_s[[i]] %>% projectRaster(., landmask_r)
  fertl_c4ann_r <- fertl_c4ann_s[[i]] %>% projectRaster(., landmask_r)
  fertl_c4per_r <- fertl_c4per_s[[i]] %>% projectRaster(., landmask_r)
  
  irrig_c3ann_r <- irrig_c3ann_s[[i]] %>% projectRaster(., landmask_r)
  irrig_c3nfx_r <- irrig_c3nfx_s[[i]] %>% projectRaster(., landmask_r)
  irrig_c3per_r <- irrig_c3per_s[[i]] %>% projectRaster(., landmask_r)
  irrig_c4ann_r <- irrig_c4ann_s[[i]] %>% projectRaster(., landmask_r)
  irrig_c4per_r <- irrig_c4per_s[[i]] %>% projectRaster(., landmask_r)
  
  fcro_r <- c3ann_r + c3nfx_r + c3per_r + c4ann_r + c4per_r
  fertl_r <- (fertl_c3ann_r * c3ann_r + fertl_c3nfx_r * c3nfx_r + fertl_c3per_r * c3per_r +
    fertl_c4ann_r * c4ann_r + fertl_c4per_r * c4per_r) / fcro_r / 10 * landmask_r ## kg ha-1 yr-1 --> g m-2 yr-1
  
  irrig_r <- (irrig_c3ann_r * c3ann_r + irrig_c3nfx_r * c3nfx_r + irrig_c3per_r * c3per_r +
                irrig_c4ann_r * c4ann_r + irrig_c4per_r * c4per_r) / fcro_r * landmask_r
  irrig_r[irrig_r > 0] <- 1
  
  ## write fertilizer
  b_vec <- fertl_r %>% values %>% as.numeric
  b_vec[is.na(b_vec)] <- 0
  
  bin_con <- file(paste0("nfer_", 849+i, ".bin") %>% file.path("nfer", .), "wb") 
  writeBin(b_vec, bin_con, size=4)
  close(bin_con)
  
  ## write irrigation
  b_vec <- irrig_r %>% values %>% as.numeric
  b_vec[is.na(b_vec)] <- 0
  # b_vec <- round(b_vec)
  
  bin_con <- file(paste0("irr_", 849+i, ".bin") %>% file.path("irr", .), "wb") 
  writeBin(b_vec, bin_con, size=4)
  close(bin_con)
  
  cat(i, "\n")
})

for(i in 851:nlayers(urban_s)) {
  urban_r <- urban_s[[i]] %>% projectRaster(., landmask_r)
  urban_r <- urban_r * landmask_r
  
  b_vec <- urban_r %>% values %>% as.numeric
  b_vec[is.na(b_vec)] <- 0
  
  bin_con <- file(paste0("imp_", 849+i, ".bin") %>% file.path("cohort/fimp", .), "wb") 
  writeBin(b_vec, bin_con, size=4)
  close(bin_con)
}

## write co2 data
co2 <- read.csv("data/global_co2_ann_1700_2017.csv")
for(i in 1:nrow(co2)) {
  b_vec <- rep(co2$CO2[i], 354*720*365)
  
  bin_con <- file(paste0("c", co2$Year[i], ".bin") %>% file.path("co2/", .), "wb") 
  writeBin(b_vec, bin_con, size=4)
  close(bin_con)
  
}

## generate new masks
mask_ori <- readDlemBin("mask/mask_new.bin")
# maskidx_ori <- readDlemBin("mask/maskidx_new.bin")
ind_eff <- which(mask_ori[]==1) 
ind_each <- 0:(length(ind_eff) - 1) %/% 15620
for(i in unique(ind_each)) {
  mask_i <- raster(mask_ori)
  mask_i[] <- -9999
  mask_i[ind_eff[ind_each == i]] <- 1
  mask_i_con <- file(paste0("mask/mask_", i, ".bin"), "wb")
  writeBin(mask_i[], mask_i_con, size=4)
  close(mask_i_con)
}








