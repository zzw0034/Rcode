library(raster)
library(sp)
library(ggplot2)
data <- raster("D:/Research_AU_EDGE/Rcode/stat_lat/prep2015.tif")


zones <- init(data, v='y')
longzones <- init(data, v='x')

z <- zonal(data, zones, 'mean')
z_sd <- zonal(data, zones, 'sd')

longz <- zonal(data, longzones, 'mean')
longz_sd <- zonal(data, longzones, 'sd')
#------------------------------------------------------------------

plot(z)
plot(z_sd)

#make it as dataframe for ggplot
z_df <- as.data.frame(z) 
z_sd_df <- as.data.frame(z_sd)

#ggplot with 
lat_plot <- ggplot(z_df) +
  geom_line(aes(x=zone,y=mean),color= "red") + 
  geom_ribbon(aes(x=zone,ymin=z_df$mean-z_sd_df$sd, ymax=z_df$mean+z_sd_df$sd), 
              linetype=2, alpha=0.1,fill= "green") +
  coord_flip()

longz_df <- as.data.frame(longz) 
longz_sd_df <- as.data.frame(longz_sd)

#ggplot with 
long_plot <- ggplot(longz_df) +
  geom_line(aes(x=zone,y=mean),color= "red") + 
  geom_ribbon(aes(x=zone,ymin=longz_df$mean-longz_sd_df$sd, ymax=longz_df$mean+longz_sd_df$sd), 
              linetype=2, alpha=0.1,fill= "green") + 
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid"))

ggsave("D:/Research_AU_EDGE/Rcode/stat_lat/image.tiff", height=3, width=7, units='in', dpi=600)
