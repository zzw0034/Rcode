## Title: trendy_climate
## Date:2018-08-23

library(raster)
library(plyr)
library(reshape2)
library(rgdal)
library(magrittr)
library(stringr)
library(parallel)
library(ncdf4)

## set workspace
ori_dir <- setwd("S:/hao_shi/gcp18/trendyv7/clim/pre/")

## custom read daily binary file
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

## write nc file
## define dimensions
lon <- seq(-179.75, 179.75, 0.5)
lon_dim <- ncdim_def("lon", "degrees_east", as.double(lon))

lat <- seq(89.75, -89.75, -0.5)
lat_dim <- ncdim_def("lat", "degrees_north", as.double(lat))

time <- 1:(118 * 12) - 1
time_dim <- ncdim_def("time", "months since 1900-01-15 00:00:00", time, calendar="noleap")

fill_value <- -99999.0

##define variables
var_name <- c("pr")
var_unit <- "kg m-2 s-1"
var_lname <- "precipitation"

var_def <- ncvar_def(var_name, var_unit, list(lon_dim, lat_dim, time_dim), fill_value,
                     var_lname, "float", compression=9)

## create netcdf
nc_file <- "E:\\Trendy_2018\\DLEM_S0_pr.nc" ## note: set up the path for nc_file
nc <- nc_create(nc_file, list(var_def), force_v4=TRUE)


var_arr <- array(numeric(), dim=c(length(lon), length(lat), length(time)))

## for s0, s1, s4
pre_s_daily_1700 <- readDlemDailyBin("pre1700.bin")
pre_s_m_1700 <- stackApply(pre_s_daily_1700, c(rep(1, 31), rep(2, 28), rep(3, 31), rep(4, 30),
                                               rep(5, 31), rep(6, 30), rep(7, 31), rep(8, 31),
                                               rep(9, 30), rep(10, 31), rep(11, 30), rep(12, 31)), mean)
pre_s_m_1700 <- pre_s_m_1700 / (24*3600)

for(iyear in 1:118) {#118
  for(imon in 1:12) {
    var_arr[,,(iyear-1)*12 + imon] <- pre_s_m_1700[[imon]] %>% values %>%
      c(rep(NA, 3*720), ., rep(NA, 3*720))
  } 
}

ncvar_put(nc, var_def, var_arr)

## put additional attributes
ncatt_put(nc, "lon", "axis", "X")
ncatt_put(nc, "lon", "long_name", "longitude")

ncatt_put(nc, "lat", "axis", "Y")
ncatt_put(nc, "lat", "long_name", "latitude")

## add global atts
ncatt_put(nc, 0, "title", "DLEM output for TRENDYv7, 2018")
ncatt_put(nc, 0, "contact", "Hanqin Tian, tianhan@auburn.edu; Hao Shi, hzs0087@auburn.edu")

## close netcdf
nc_close(nc)
