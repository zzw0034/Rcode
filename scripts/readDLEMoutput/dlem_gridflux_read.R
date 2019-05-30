
## import library
require(magrittr)
require(raster)

mask_file <- "S:\\hao_shi\\gcp18\\trendyv7\\mask/mask_3.bin"
bin_file <- "E:\\R\\Test/yy1701.dlem"
###########################################################################################
## custom function
readDlemBin <- function(bin_file) {
  bin_con <- file(bin_file, "rb")
  data_nu <- readBin(bin_con, "double", n=354 * 720, size=4)
  close(bin_con)
  
  r <- matrix(data_nu, nrow=354, byrow=TRUE) %>% 
    raster(., xmn=-180, xmx=180, ymn=-88.5, ymx=88.5, crs="+proj=longlat +datum=WGS84")
  return(r)
}
######################################################
## read dlem gridflux using raw method
fileSize <- file.info(bin_file)$size
raw <- readBin(bin_file, what="raw", n=fileSize);

nbrOfRecords <- length(raw) %/% 79;
raw <-  matrix(raw, nrow=79*4, byrow=FALSE)

data_int <- lapply(1:5, FUN=function(var) {
  ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                  n=nbrOfRecords/4)
  return(ints)
})

data_float <- lapply(6:79, FUN=function(var) {
  floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                  n=nbrOfRecords/4)
  return(floats)
})

gridflux_df <- data.frame(do.call(cbind, data_int), do.call(cbind, data_float))
var_names <- c("Findex_auto", "Findex_input", "year", "month", "day", 
                        "g_gpp", "g_npp", "g_mr", "g_gr", "g_rh",
                        "g_harvc", "g_harvn", "g_estabC", "g_estabN", "g_proddecC", "g_lucc_cvrtC",
                        "g_lucc_cvrtN", "g_crop_manage_C", "g_crop_manage_N", "g_manureC", "g_manureN", 
                        "C_burntflux", "CH4_burntflux", "CO_burntflux", "CO2_burntflux", "NMHC_burntflux", 
                        "C_burntmort", "g_litleachDOC", "g_leach_DOC",
                        "g_leach_POC", "g_leach_DIC","g_ch4", "g_VOCisop", "g_VOCmono", "g_VOCovoc",
                        "g_VOCorvoc", "g_VOCco","g_nuptake", "g_NMin", "g_leakN", "g_leachNH4", 
                        "g_leachNO3", "g_litleachDON", "g_leachDON", "g_leachPON", 
                        "g_denitrif", "g_nitrif", "g_nFix", "g_Ndep", "g_NFer", "g_n2o", "g_no",
                        "g_nh3Emisn", "N_burntflux", "N2O_burntflux", "NOy_burntflux", "NH3_burntflux",
                        "N_burntmort","g_irrigation", "g_irrigation_withd", "g_evap",
                        "g_infil", "g_trans","g_runoff_drain","g_runoff_surf",
                        "g_stream_in", "g_stream_out","g_stream_et","g_stream_irrig",
                        "g_stream_outDOC", "g_stream_outPOC","g_stream_outNH4",
                        "g_stream_outNO3","g_stream_outDIN", "g_stream_outDON", "g_stream_outPON",
                        "g_stream_outDIC", "g_stream_rh","g_stream_nremoval"
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
# writeRaster()

#==readDlemGridflux===============================================================================================
readDlemGridflux <- function(mask_file,dlem_file) {
  ## custom function
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
  
  nbrOfRecords <- length(raw) %/% 79;
  raw <-  matrix(raw, nrow=79*4, byrow=FALSE)
  
  data_int <- lapply(1:5, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbrOfRecords/4)
    return(ints)
  })
  
  data_float <- lapply(6:79, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbrOfRecords/4)
    return(floats)
  })
  
  gridflux_df <- data.frame(do.call(cbind, data_int), do.call(cbind, data_float))
  var_names <- c("Findex_auto", "Findex_input", "year", "month", "day", 
                 "g_gpp", "g_npp", "g_mr", "g_gr", "g_rh",
                 "g_harvc", "g_harvn", "g_estabC", "g_estabN", "g_proddecC", "g_lucc_cvrtC",
                 "g_lucc_cvrtN", "g_crop_manage_C", "g_crop_manage_N", "g_manureC", "g_manureN", 
                 "C_burntflux", "CH4_burntflux", "CO_burntflux", "CO2_burntflux", "NMHC_burntflux", 
                 "C_burntmort", "g_litleachDOC", "g_leach_DOC",
                 "g_leach_POC", "g_leach_DIC","g_ch4", "g_VOCisop", "g_VOCmono", "g_VOCovoc",
                 "g_VOCorvoc", "g_VOCco","g_nuptake", "g_NMin", "g_leakN", "g_leachNH4", 
                 "g_leachNO3", "g_litleachDON", "g_leachDON", "g_leachPON", 
                 "g_denitrif", "g_nitrif", "g_nFix", "g_Ndep", "g_NFer", "g_n2o", "g_no",
                 "g_nh3Emisn", "N_burntflux", "N2O_burntflux", "NOy_burntflux", "NH3_burntflux",
                 "N_burntmort","g_irrigation", "g_irrigation_withd", "g_evap",
                 "g_infil", "g_trans","g_runoff_drain","g_runoff_surf",
                 "g_stream_in", "g_stream_out","g_stream_et","g_stream_irrig",
                 "g_stream_outDOC", "g_stream_outPOC","g_stream_outNH4",
                 "g_stream_outNO3","g_stream_outDIN", "g_stream_outDON", "g_stream_outPON",
                 "g_stream_outDIC", "g_stream_rh","g_stream_nremoval"
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




# Another method, very slow
# bin_file <- "D:/hao/gcp2018_trendyv7/m3s0yy1701.dlem"
# bin_con <- file(bin_file, "rb")
# 
# i=1;
# flag = 1;
# 
# gridflux <- c()
# 
# while(flag > 0) {
#   if(i < 6) {
#     what_type <- "integer"
#   } else {
#     what_type <- "double"
#   }
#   
#   tmp <- readBin(bin_con, what_type, size=4, n=1)
#   if(length(tmp)==0) {
#     flag=0
#   }
#   
#   cat(i, "\n")
#   gridflux <- c(gridflux, tmp)
#   
#   i=i+1
#   if(i > 79) i=1
# }
# 
# close(bin_con)
# 
# ## convert vec to df
# gridflux_df <- matrix(gridflux, ncol=79, byrow=TRUE) %>% data.frame





