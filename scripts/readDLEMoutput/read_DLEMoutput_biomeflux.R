require(magrittr)
require(raster)
#==result <- readDlemGridflux()======
#==readDlemGridflux=============================================
readDlemBiomeflux <- function(mask_file,dlem_file) {
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
  
  nbrOfRecords <- length(raw) %/% 61;
  raw <-  matrix(raw, nrow=61*4, byrow=FALSE)
  
  data_int <- lapply(1:7, FUN=function(var) {
    ints <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="integer", size=4,
                    n=nbrOfRecords/4)
    return(ints)
  })
  
  data_float <- lapply(8:61, FUN=function(var) {
    floats <- readBin(con=as.vector(raw[(4*(var - 1) + 1): (4*var), ]), what="double", size=4,
                      n=nbrOfRecords/4)
    return(floats)
  })
  
  gridflux_df <- data.frame(do.call(cbind, data_int), do.call(cbind, data_float))
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
  names(gridflux_df) <- var_names
  rm(raw)
  
  ## read mask file
  mask_r <- readDlemBin(mask_file)
  
  ## generate raster for a var at a day/mon/year
  genVarStack <- function(id) {
    gridflux_df_ecoid <- subset(gridflux_df, ecoid==id)
    ecoid_s <- lapply(1:ncol(gridflux_df_ecoid), FUN=function(col) {
      var_r <- raster(mask_r)
      ind <- which((mask_r[] %in% gridflux_df_ecoid$Findex_input)==TRUE)
      var_r[ind] <- gridflux_df_ecoid[, col]
      
      return(var_r)
    })
    
    ecoid_s <- do.call(stack, ecoid_s)
    cat(id, "\n")
    return(ecoid_s)
  }
  
  gridflux_df_s <- lapply(unique(gridflux_df$ecoid), genVarStack) %>% do.call(stack, .)
  names(gridflux_df_s) <- outer(var_names, unique(gridflux_df$ecoid), paste, sep="_") %>% as.vector
  return(gridflux_df_s)
}

#-------------------------------------------------------------------------------------------------------------------
mask_file <- "S:\\hao_shi\\gcp18\\trendyv7\\mask/maskidx_new.bin"
dlem_file <- "E:\\R\\Test/biomeflux_yy1702.dlem"
