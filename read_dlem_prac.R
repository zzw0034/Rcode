## Title: read_dlem
## Author: hzs0087@auburn.edu

## import packages
library(raster)
library(magrittr)
library(rgdal)

## custom functions
## read dlem binary input files
if(!exists("readDlemBin")) {
  readDlemBin <- function(bin_file) {
    bin_con <- file(bin_file, "rb")
    data_nu <- readBin(bin_con, "double", n=354 * 720, size=4)
    close(bin_con)
    r <- matrix(data_nu, nrow=354, byrow=TRUE) %>% 
      raster(., xmn=-180, xmx=180, ymn=-88.5, ymx=88.5, crs="+proj=longlat +datum=WGS84")
    return(r)
  }  
}

## read dlem biomeflux output
readDlemBiomeflux <- function(maskidx_file, dlem_biomeflux_file) {
  
  ## read dlem biomeflux using raw method
  file_size <- file.info(dlem_biomeflux_file)$size
  raw <- readBin(dlem_biomeflux_file, what="raw", n=file_size);
  
  ## wrap vector into matrix to accerate processing 
  nbr_of_records <- length(raw) %/% 61;
  raw <-  matrix(raw, nrow=61*4, byrow=FALSE)
  
  ## read the first 7 variables of interger type
  data_int <- lapply(1:7, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbr_of_records/4)
    return(ints)
  })
  
  ## read the other variables of float type
  data_float <- lapply(8:61, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbr_of_records/4)
    return(floats)
  })
  
  biomeflux_df <- data.frame(do.call(cbind, data_int), do.call(cbind, data_float))
  
  ## var list
  var_names <- c("Findex_auto", "Findex_input", "year", "month", "day", 
                 "ecoid","pft","g_gpp", "g_npp", "g_mr", "g_gr", "g_rh",
                 "g_harvc", "g_harvn", "g_estabC", "g_estabN", 
                 "crope_manage_C","crope_manage_N","g_manureC","g_manureN",
                 "C_burntflux", "CH4_burntflux", "CO_burntflux", "CO2_burntflux", "NMHC_burntflux", 
                 "C_burntmort", "g_litleachDOC", "g_leach_DOC",
                 "g_leach_POC", "g_leach_DIC","g_ch4", "g_VOCisop", "g_VOCmono", "g_VOCovoc",
                 "g_VOCorvoc", "g_VOCco","g_nuptake", "g_NMin", "g_leakN", "g_leachNH4", 
                 "g_leachNO3", "g_litleachDON", "g_leachDON", "g_leachPON", 
                 "g_denitrif", "g_nitrif", "g_nFix", "g_Ndep", "g_NFer", "g_n2o", "g_no",
                 "g_nh3Emisn", "N_burntflux", "N2O_burntflux", "NOy_burntflux", "NH3_burntflux",
                 "N_burntmort","g_irrigation", "g_evap","g_intercp","g_trans"
  )
  names(biomeflux_df) <- var_names
  rm(raw)
  
  ## read maskidx file
  maskidx_r <- readDlemBin(maskidx_file)
  
  ## generate raster for a var at a day/mon/year
  genVarStack <- function(id) {
    biomeflux_df_ecoid <- subset(biomeflux_df, ecoid==id)
    ecoid_s <- lapply(1:ncol(biomeflux_df_ecoid), FUN=function(col) {
      var_r <- maskidx_r
      var_r[var_r >=0] <- 0
      var_r[var_r < 0] <- NA
      ind <- which((maskidx_r[] %in% biomeflux_df_ecoid$Findex_input)==TRUE)
      var_r[ind] <- biomeflux_df_ecoid[, col]
      
      return(var_r)
    })
    
    ecoid_s <- do.call(stack, ecoid_s)
    cat(id, "\n")
    return(ecoid_s)
  }
  
  biomeflux_df_s <- lapply(unique(biomeflux_df$ecoid), genVarStack) %>% do.call(stack, .)
  names(biomeflux_df_s) <- outer(var_names, unique(biomeflux_df$ecoid), paste, sep="_") %>% as.vector
  
  return(biomeflux_df_s)
}

## read dlem gridflux output
readDlemGridflux <- function(maskidx_file, dlem_gridflux_file) {
  file_size <- file.info(dlem_gridflux_file)$size
  raw <- readBin(dlem_gridflux_file, what="raw", n=file_size);
  
  nbr_of_records <- length(raw) %/% 80;
  raw <-  matrix(raw, nrow=80*4, byrow=FALSE)
  
  data_int <- lapply(1:5, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbr_of_records/4)
    return(ints)
  })
  
  data_float <- lapply(6:80, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbr_of_records/4)
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
                 "g_nh3Emisn", "g_n2", "N_burntflux", "N2O_burntflux", "NOy_burntflux", "NH3_burntflux",
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
  maskidx_r <- readDlemBin(maskidx_file)
  
  ## generate raster for a var at a day/mon/year
  genVarRaster <- function(col) {
    var_r <- maskidx_r
    var_r[var_r >= 0] <- 0
    var_r[var_r < 0] <- NA
    ind <- which((maskidx_r[] %in% gridflux_df$Findex_input)==TRUE)
    var_r[ind] <- gridflux_df[, col]
    
    return(var_r)
  }
  
  gridflux_df_s <- lapply(1:ncol(gridflux_df), genVarRaster) %>% do.call(stack, .)
  names(gridflux_df_s) <- var_names
  
  return(gridflux_df_s)
}

## read dlem gridpool output
readDlemGridpool <- function(maskidx_file, dlem_gridpool_file) {
  
  ## read dlem gridpool using raw method
  file_size <- file.info(dlem_gridpool_file)$size
  raw <- readBin(dlem_gridpool_file, what="raw", n=file_size);
  
  nbr_of_records <- length(raw) %/% 20;
  raw <- matrix(raw, nrow=20*4, byrow=FALSE)
  
  data_int1 <- lapply(1:2, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbr_of_records/4)
    return(ints)
  })
  
  data_float1 <- lapply(3, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbr_of_records/4)
    return(floats)
  })
  
  gridpool_df1 <- data.frame(do.call(cbind, data_int1),do.call(cbind, data_float1))
  
  data_int2 <- lapply(4:6, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbr_of_records/4)
    return(ints)
  })
  
  data_float2 <- lapply(7:20, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbr_of_records/4)
    return(floats)
  })
  
  gridpool_df2 <- data.frame(do.call(cbind, data_int2),do.call(cbind, data_float2))
  gridpool_df <- data.frame(gridpool_df1,gridpool_df2)
  
  var_names <- c("Findex_auto", "Findex_input","Farea", "year", "month", "day", 
                 "VegC_rec", "LitC_rec", "SoilC_rec", "CWDC_rec","LAIgrid_rec", "WoodagC_rec",
                 "VegN_rec", "LitN_rec", "SoilN_rec", "CWDN_rec", "FL_fine_rec", "FL_coarse_rec",
                 "FL_livewood_rec", "FL_liveherb_rec")
  
  names(gridpool_df) <- var_names
  rm(raw)
  
  ## read mask file
  maskidx_r <- readDlemBin(maskidx_file)
  
  ## generate raster for a var at a day/mon/year
  genVarRaster <- function(col) {
    var_r <- maskidx_r
    var_r[var_r >=0] <- 0
    var_r[var_r < 0] <- NA
    
    ind <- which((maskidx_r[] %in% gridpool_df$Findex_input)==TRUE)
    var_r[ind] <- gridpool_df[, col]
    
    return(var_r)
  }
  
  gridpool_df_s <- lapply(1:ncol(gridpool_df), genVarRaster) %>% do.call(stack, .)
  names(gridpool_df_s) <- var_names
  
  return(gridpool_df_s)
}

## read dlem biomestate output
readDlemBiomestate <- function(maskidx_file, dlem_biomestate_file) {
  
  ## read dlem biomeflux using raw method
  file_size <- file.info(dlem_biomestate_file)$size
  raw <- readBin(dlem_biomestate_file, what="raw", n=file_size);
  
  ## wrap vector into matrix to accerate processing 
  nbr_of_records <- length(raw) %/% 55;
  raw <-  matrix(raw, nrow=55*4, byrow=FALSE)
  
  ## read the first 7 variables of interger type
  data_int <- lapply(1:7, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbr_of_records/4)
    return(ints)
  })
  
  ## read the other variables of float type
  data_float <- lapply(8:55, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbr_of_records/4)
    return(floats)
  })
  
  biomestate_df <- data.frame(do.call(cbind, data_int), do.call(cbind, data_float))
  
  ## var list
  var_names <- c("Findex_auto", "Findex_input", "year", "month", "day", 
                 "ecoid","pft","fbiomec", "ecoage", "age_disturb", "cleafmax_pot", "laiphenmax",
                 "g_suno3_pday", "g_shao3_pday", "canopywater", "leaf_c", 
                 "sapw_c","heartw_c","froot_c","croot_c",
                 "reprod_c", "wdag_c", "wdbg_c", "aom1ag_c", "aom1bg_c", 
                 "aom2ag_c", "aom2bg_c", "smb1_c",
                 "smb2_c", "smr_c","nom_c", "psom_c", "soil_doc", "soil_dic",
                 "soil_ch4", "leaf_n","sapw_n", "heartw_n", "froot_n", "croot_n", 
                 "reprod_n", "wdag_n", "wdbg_n", "aom1ag_n", 
                 "aom1bg_n", "aom2ag_n", "aom2bg_n", "smb1_n", "smb2_n", "smr_n", "nom_n",
                 "psom_n", "soil_don", "av_NH4", "av_NO3"
  )
  names(biomestate_df) <- var_names
  rm(raw)
  
  ## read maskidx file
  maskidx_r <- readDlemBin(maskidx_file)
  
  ## generate raster for a var at a day/mon/year
  genVarStack <- function(id) {
    biomestate_df_ecoid <- subset(biomestate_df, ecoid==id)
    ecoid_s <- lapply(1:ncol(biomestate_df_ecoid), FUN=function(col) {
      var_r <- maskidx_r
      var_r[var_r >=0] <- 0
      var_r[var_r < 0] <- NA
      ind <- which((maskidx_r[] %in% biomestate_df_ecoid$Findex_input)==TRUE)
      var_r[ind] <- biomestate_df_ecoid[, col]
      
      return(var_r)
    })
    
    ecoid_s <- do.call(stack, ecoid_s)
    cat(id, "\n")
    return(ecoid_s)
  }
  
  biomestate_df_s <- lapply(unique(biomestate_df$ecoid), genVarStack) %>% do.call(stack, .) ##unique(biomestate_df$ecoid)
  names(biomestate_df_s) <- outer(var_names, unique(biomestate_df$ecoid), paste, sep="_") %>% as.vector ##unique(biomestate_df$ecoid)
  
  return(biomestate_df_s)
}
