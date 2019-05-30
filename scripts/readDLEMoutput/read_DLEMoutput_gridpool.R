####################
#read DLEM gridpool#
####################
require(magrittr)
require(raster)
#==result <- readDlemGridpool()======
#==readDlemGridflux=============================================
readDlemGridpool <- function(mask_file,dlem_file) {
  ## custom function readDlemBin:
  readDlemBin <- function(bin_file) {
    bin_con <- file(bin_file, "rb")
    data_nu <- readBin(bin_con, "double", n=354 * 720, size=4)
    close(bin_con)
    r <- matrix(data_nu, nrow=354, byrow=TRUE) %>% 
      raster(., xmn=-180, xmx=180, ymn=-88.5, ymx=88.5, crs="+proj=longlat +datum=WGS84")
    return(r)
  }
  
  ## read dlem gridflux using raw method
  fileSize <- file.info(dlem_file)$size
  raw <- readBin(dlem_file, what="raw", n=fileSize);
  
  nbrOfRecords <- length(raw) %/% 20;
  
  raw <-  matrix(raw, nrow=20*4, byrow=FALSE)
  
  data_int1 <- lapply(1:2, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbrOfRecords/4)
    return(ints)
  })
  
  data_float1 <- lapply(3, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbrOfRecords/4)
    return(floats)
  })
  gridflux_df1 <- data.frame(do.call(cbind, data_int1),do.call(cbind, data_float1))
  
  data_int2 <- lapply(4:6, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbrOfRecords/4)
    return(ints)
  })
  
  data_float2 <- lapply(7:20, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbrOfRecords/4)
    return(floats)
  })
  
  gridflux_df2 <- data.frame(do.call(cbind, data_int2),do.call(cbind, data_float2))
  gridflux_df <- data.frame(gridflux_df1,gridflux_df2)
  
  var_names <- c("Findex_auto", "Findex_input","Farea", "year", "month", "day", 
                 "VegC_rec", "LitC_rec", "SoilC_rec", "CWDC_rec","LAIgrid_rec", "WoodagC_rec",
                 "VegN_rec", "LitN_rec", "SoilN_rec", "CWDN_rec", "FL_fine_rec", "FL_coarse_rec",
                 "FL_livewood_rec", "FL_liveherb_rec"
  )
  names(gridflux_df) <- var_names
  rm(raw)
  
  ## read mask file
  mask_r <- readDlemBin(mask_file)
  
  ## generate raster for a var at a day/mon/year
  genVarRaster <- function(col) {
    var_r <- raster(mask_r)
    var_r[mask_r[] != -9999] <- gridflux_df[, col]
    return(var_r)
  }
  
  gridflux_df_s <- lapply(1:ncol(gridflux_df), genVarRaster) %>% do.call(stack, .)
  names(gridflux_df_s) <- var_names
  return(gridflux_df_s)
}

#-------------------------------------------------------------------------------------------------------------------
mask_file <- "S:\\hao_shi\\gcp18\\trendyv7\\mask/mask_3.bin"
dlem_file <- "E:\\R\\Test/m3s0yy1701.dlem"
