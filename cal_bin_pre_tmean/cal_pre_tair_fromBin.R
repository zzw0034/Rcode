
library(raster)
library(magrittr)

setwd("D:/Research_AU_EDGE/Rcode/cal_bin_pre_tmean/tair_mean/")

readDlemClim <- function(in_file) {
  s <- brick()
  data_nu <- readBin(file(in_file, "rb"), "double", n=640 * 608 * 365, size=4)
  for(day in 1:365) {
    mat_day <- data_nu[(640 * 608 * (day - 1)  + 1):(640 * 608 * day)]
    mat_day <- matrix(mat_day, nrow=640, byrow=TRUE)
    r_day <- raster(mat_day, xmn=-20, xmx=56, ymn=-40, ymx=40,crs="+proj=longlat +datum=WGS84")
    r_day[r_day < -1000] <- NA
    s <- addLayer(s, r_day)
  }
  
  return(s)
}

files <- list.files("//FW-2UA2020YRM/Yao_correct/climate/tair")[81:116]

for ( file in files) {
  tair <- readDlemClim(paste0("//FW-2UA2020YRM/Yao_correct/climate/tair/",file))
  tair_mean <- mean(tair)
  writeRaster(tair_mean,gsub('.bin','_mean.tif',file),options=c('TFW=YES'))
}


setwd("D:/Research_AU_EDGE/Rcode/cal_bin_pre_tmean/prec_year/")
files <- list.files("//FW-2UA2020YRM/Yao_correct/climate/prec")[81:116]

for ( file in files) {
  prec <- readDlemClim(paste0("//FW-2UA2020YRM/Yao_correct/climate/prec/",file))
  prec_annual <- sum(prec)
  writeRaster(prec_annual,gsub('.bin','_annual.tif',file),options=c('TFW=YES'))
}




# y1900 <- readDlemClim("tair1900.bin")
# y1900_mean <- mean(y1900)
# writeRaster(y1900_mean,paste0('y1900_mean','.tif'),options=c('TFW=YES'))



# for(var in c("dswrf", "pre", "tmax", "tmin", "tmp")) {
#   avg <- 0
#   for(year in 1901:1920){
#     var_con <- file(paste0(var, "/", var, year, ".bin"), "rb")
#     var_bin <- readBin(var_con, "double", n=354 * 720 * 365, size=4)
#     avg <- var_bin + avg
#     close(var_con)
#     cat(year, "\n")
#   }
#   
#   avg = avg / 20
#   
#   for(year in 1700:1900) {
#     var_con_year <- file(paste0(var, "/", var, year, ".bin"), "wb")
#     writeBin(avg, var_con_year, size = 4)
#     close(var_con_year)
#   }
# }
