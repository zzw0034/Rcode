
#average 20 yrs climate data .bin;
setwd("S:/hao_shi/gcp18/trendyv7/clim/")

for(var in c("dswrf", "pre", "tmax", "tmin", "tmp")) {
  avg <- 0
  for(year in 1901:1920){
    var_con <- file(paste0(var, "/", var, year, ".bin"), "rb")
    var_bin <- readBin(var_con, "double", n=354 * 720 * 365, size=4)
    avg <- var_bin + avg
    close(var_con)
    cat(year, "\n")
  }
  
  avg = avg / 20
  
  for(year in 1700:1900) {
    var_con_year <- file(paste0(var, "/", var, year, ".bin"), "wb")
    writeBin(avg, var_con_year, size = 4)
    close(var_con_year)
  }
}

