#ues longtitude and latitude extract point climate data from global dataset, CRUNCEP
library(raster)
library(plyr)
library(reshape2)
library(rgdal)
library(magrittr)
library(stringr)
library(parallel)
# library(ncdf4)
## custom read daily binary file function
readDlemDailyBin <- function(bin_file) {
  s <- brick()
  bin_con <- file(bin_file, "rb")
  data_nu <- readBin(bin_con, "double", n=354 * 720 * 365, size=4)
  close(bin_con)
  
  for(day in 1:365) {
    # cat(day, "\n")
    mat_day <- data_nu[(354 * 720 * (day - 1)  + 1):(354 * 720 * day)]
    mat_day <- matrix(mat_day, nrow=354, byrow=TRUE)
    r_day <- raster(mat_day, xmn=-180, xmx=180, ymn=-88.5, ymx=88.5)
    proj4string(r_day) <- "+proj=longlat +datum=WGS84"
    r_day[r_day < -1000] <- NA
    s <- addLayer(s, r_day)
  }
  return(s)
}
### custom round 0.5 function--------------------------------------------------------------
round_0.5 <-  function(x){
t <- ifelse((x - floor(x)) <= 0.25,  floor(x),
            ifelse(((x - floor(x)) > 0.25 & (x - floor(x)) <= 0.75),  (floor(x)+0.5), 
                   ifelse((x - floor(x)) > 0.75 ,  (floor(x)+1) ,x )   
                   )  
              )

  return(t)
}
#end---------------------------------------------------------------------------------------


## set workspace
ori_dir <- setwd("S:/Trendy_2017/climate_CRUNCEP/prec")

pre_files <- list.files(pattern="\\.bin$") #157 years

#???????????????
#convert long lat to sp point
lat <- 45.07
long <- -123.93

site <- data.frame(x = round_0.5(long) ,y = round_0.5(lat))
coordinates(site) <- c("x", "y")
proj4string(site) <- "+proj=longlat +datum=WGS84"
#==============================================================
#read in csv then get all long lat
sites_table <-  read.csv("D:/Research_AU_EDGE/Rcode/site+ID.csv")
sites <- sites_table[,c(1,4,5)]
sites$Lat <- round_0.5(sites$Lat.)
sites$test <- round_0.5(sites$Long.)
#_________________________________________________________________________________________

#?????????????????????,????????????,??????????????????raster
#??????site ??????
for(file in pre_files){
  pre_s_daily <- readDlemDailyBin(file)
  pixel_day <- as.vector(extract( pre_s_daily, site)) 
  
  bin_con <- file(paste0("D:/Research_AU_EDGE/Rcode/prec/",file), "wb")
  writeBin(pixel_day, bin_con, size = 4)
  close(bin_con)
}




