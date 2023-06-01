##### This script is supposed to take the mean of the 3 months of a season, and calculate the monthly SST, SBT and ICC at each region

library(rgdal)

# Load env layers

SST_winter <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Winter_SST.tiff")
SBT_winter <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Winter_SBT.tiff")
ICC_winter <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Winter_SIC.tiff")

SST_spring <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Spring_SST.tiff")
SBT_spring <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Spring_SBT.tiff")
ICC_spring <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Spring_SIC.tiff")

SST_summer <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Summer_SST.tiff")
SBT_summer <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Summer_SBT.tiff")
ICC_summer <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Summer_SIC.tiff")

SST_autumn <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Autumn_SST.tiff")
SBT_autumn <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Autumn_SBT.tiff")
ICC_autumn <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Envars/Autumn_SIC.tiff")

Bot_depth <- raster::raster("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Envars/Bio-Oracle/layers.grd")[[1]]

Bot_depth = -Bot_depth
Bot_depth <- raster::resample(Bot_depth, SST_winter, method='bilinear')

# Load shapefiles
bar_pol = readOGR("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Study_area_shapefiles/baarents_all.shp")
now_pol = readOGR("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Study_area_shapefiles/Norwegian Sea.shp")
nor_pol = readOGR("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/Study_area_shapefiles/North Sea.shp")

## SUMMER RASTERS #### 

lay_summer_bar = list()
for (i in 1:26){
  lay_summer_bar[[i]] = raster::mask(raster::stack(SST_summer[[i]], SBT_summer[[i]], Bot_depth, ICC_summer[[i]]),bar_pol)
  names(lay_summer_bar[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

lay_summer_nor = list()
for (i in 1:26){
  lay_summer_nor[[i]] = raster::mask(raster::stack(SST_summer[[i]], SBT_summer[[i]], Bot_depth, ICC_summer[[i]]),nor_pol)
  names(lay_summer_nor[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

## AUTUMN RASTERS ####

lay_autumn_now = list()
for (i in 1:26){
  lay_autumn_now[[i]] = raster::mask(raster::stack(SST_autumn[[i]], SBT_autumn[[i]], Bot_depth, ICC_autumn[[i]]), now_pol)
  names(lay_autumn_now[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

## Save all the raster files: 

save(lay_summer_bar,lay_autumn_now, lay_summer_nor,file="data/environmental_envelope_all.RData")


raster::plot(lay_summer_bar[[2]])


### Add SBT, SST, ICC and Depth