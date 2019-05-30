
library(raster)
library(magrittr)

setwd("D:/Research_AU_EDGE/Rcode/cal_bin_pre_tmean/tair_month/")

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

monthdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
monthstartday <- c(0,31,59,90,120,151,181,212,243,273,304,334)+1
# int monthdays[12];//{31,28,31,30,31,30,31,31,30,31,30,31} days per month (days)
# int monthstartday[12];//{0,31,59,90,120,151,181,212,243,273,304,334} start day of each month (day)

files <- list.files("//FW-2UA2020YRM/Yao_correct/climate/tair")[81:116]

for ( file in files) {
  tair_mean <- brick()
  tair <- readDlemClim(paste0("//FW-2UA2020YRM/Yao_correct/climate/tair/",file))
  for (i in 1:12) {
    # tair_mean[i] <- mean(tair[[(monthstartday[i]):(monthstartday[i]+monthdays[i]-1)]])
    tair_mean <- addLayer(tair_mean, mean(tair[[(monthstartday[i]):(monthstartday[i]+monthdays[i]-1)]]))
    writeRaster(tair_mean[[i]],gsub('.bin',paste0(i,'_mean.tif'),file),options=c('TFW=YES'))
  }
 
}


setwd("D:/Research_AU_EDGE/Rcode/cal_bin_pre_tmean/prec_month/")
files <- list.files("//FW-2UA2020YRM/Yao_correct/climate/prec")[81:116]

for ( file in files) {
  prec_month <- brick()
  prec <- readDlemClim(paste0("//FW-2UA2020YRM/Yao_correct/climate/prec/",file))
  for (i in 1:12) {
    prec_month <- addLayer(prec_month, sum(prec[[(monthstartday[i]):(monthstartday[i]+monthdays[i]-1)]]))
    writeRaster(prec_month[[i]],gsub('.bin',paste0(i,'_mean.tif'),file),options=c('TFW=YES'))
    }

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
