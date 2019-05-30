
library(ncdf4)
library(raster)
library(rasterVis)

setwd("D:/Research_AU_EDGE/R/N2O_nc_file")
list.files()
ORC_n2o <- brick("ORCHIDEE-CNP_S1_nSoil.nc4")


ORC_n2o_mean <- mean(ORC_n2o)

n2o_layer_sum <- NULL
n2o_layer_mean<- NULL

for(i in 1:1872) {
n2o_layer <- ORC_n2o[[i]]

writeRaster(n2o_layer,paste0('ORCHIDEE-CNP_S1_nSoil_',i,'.tif'),options=c('TFW=YES'))

n2o_layer_sum_temp <- cellStats(n2o_layer, "sum", na.rm=TRUE)
n2o_layer_mean_temp <- cellStats(n2o_layer, "mean",na.rm=TRUE)
n2o_layer_sum <- c(n2o_layer_sum, n2o_layer_sum_temp)
n2o_layer_mean <- c(n2o_layer_mean, n2o_layer_mean_temp)

}


levelplot(ORC_n2o_mean)

#length(n2o_layer_sum)
#length(n2o_layer_mean)
#i <- 1
#sum(!is.na(getValues(n2o_layer)))

