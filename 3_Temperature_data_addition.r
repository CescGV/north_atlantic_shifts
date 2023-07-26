##### This script is supposed to take the mean of the 3 months of a season, and calculate the monthly SST, SBT and ICC at each region

library(rgdal)

# Load env layers

SST_winter <- raster::stack("../Envars/SST_winter")
SBT_winter <- raster::stack("../Envars/SBT_winter")
ICC_winter <- raster::stack("../Envars/SIC_winter")
SST_summer <- raster::stack("../Envars/SST_summer")
SBT_summer <- raster::stack("../Envars/SBT_summer")
ICC_summer <- raster::stack("../Envars/SIC_summer")
SST_autumn <- raster::stack("../Envars/SST_autumn")
SBT_autumn <- raster::stack("../Envars/SBT_autumn")
ICC_autumn <- raster::stack("../Envars/SIC_autumn")

Bot_depth <- sdmpredictors::load_layers("BO_bathymean")
Bot_depth <- raster::resample(Bot_depth, SST_winter, method='bilinear')
Bot_depth <- -Bot_depth
Bot_depth[raster::values(Bot_depth) < 0] = NA
Bot_depth[raster::values(Bot_depth) > 600] = NA

# Load shapefiles
bar_pol = readOGR("../Envars/barents_poly.shp")
now_pol = readOGR("../Envars/norwegian_poly.shp")
nor_pol = readOGR("../Envars/north_poly.shp")

## SUMMER RASTERS #### 

lay_summer_bar = list()
for (i in 1:28){
  lay_summer_bar[[i]] = raster::mask(raster::stack(SST_summer[[i]], SBT_summer[[i]], Bot_depth, ICC_summer[[i]]),bar_pol)
  names(lay_summer_bar[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

lay_summer_nor = list()
for (i in 1:28){
  lay_summer_nor[[i]] = raster::mask(raster::stack(SST_summer[[i]], SBT_summer[[i]], Bot_depth, ICC_summer[[i]]),nor_pol)
  names(lay_summer_nor[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

## AUTUMN RASTERS ####

lay_autumn_now = list()
for (i in 1:28){
  lay_autumn_now[[i]] = raster::mask(raster::stack(SST_autumn[[i]], SBT_autumn[[i]], Bot_depth, ICC_autumn[[i]]), now_pol)
  names(lay_autumn_now[[i]]) = c("SST", "SBT", "Bot_depth", "ICC")
  print(i)
}

## Save all the raster files: 

save(lay_summer_bar,lay_autumn_now, lay_summer_nor,file="../data/environmental_envelope_all.RData")


raster::plot(lay_summer_nor[[2]])


### Add SBT, SST, ICC and Depth

