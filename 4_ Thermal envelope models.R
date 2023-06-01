
## 1. Add environmental data ####
library(tidyverse)
load("data/environmental_envelope_all.RData")
raw_data = read.csv("data/data_ready.csv") %>% filter(Year < 2020 & Year > 1993)
dat_yr_env = NULL

for (i in 1:length(unique(raw_data$Region))){
  dat = dplyr::filter(raw_data, raw_data$Region == unique(raw_data$Region)[i])
  if(unique(dat$Region) == "Barents Sea"){
    for (k in sort(unique(dat$Year))){
      dat_yr = dplyr::filter(dat, dat$Year == sort(unique(dat$Year))[k-min(dat$Year)+1]) 
      env = data.frame(raster::extract(lay_summer_bar[[k-1993]],data.frame(x = dat_yr$Longitude,
                                                                           y = dat_yr$Latitude)))
      dat_yr_env = rbind(dat_yr_env, cbind(dat_yr, env))
    }
  } else if(unique(dat$Region) == "Norwegian Sea"){
    for (k in sort(unique(dat$Year))){
      dat_yr = dplyr::filter(dat, dat$Year == sort(unique(dat$Year))[k-min(dat$Year)+1]) 
      env = data.frame(raster::extract(lay_autumn_now[[k-1993]],data.frame(x = dat_yr$Longitude,
                                                                           y = dat_yr$Latitude)))
      dat_yr_env = rbind(dat_yr_env, cbind(dat_yr, env))
    }
  } else if(unique(dat$Region) == "North Sea") {
    for (k in sort(unique(dat$Year))){
      dat_yr = dplyr::filter(dat, dat$Year == sort(unique(dat$Year))[k-min(dat$Year)+1]) 
      env = data.frame(raster::extract(lay_summer_nor[[k-1993]],data.frame(x = dat_yr$Longitude,
                                                                           y = dat_yr$Latitude)))
      dat_yr_env = rbind(dat_yr_env, cbind(dat_yr, env))
    }
  }
}
colSums(is.na(dat_yr_env))
#The NAs in sea ice are actually 0es

dat_yr_env$ICC[is.na(dat_yr_env$ICC)] = 0
write.csv(dat_yr_env, "data/shifts_and_env_data.csv")

## 2. Prepare PA dataset ####
library(mgcv)
raw_data = dat_yr_env[,-46]
PA_data = raw_data
AB_data = raw_data
AB_data = AB_data %>% group_by(Year, Species) %>% mutate(annual_mean_biomass = mean(Biomass,na.rm = T))

# All sites
PA_sites = PA_data %>% select(Longitude, Latitude) %>% distinct() %>% 
  mutate(Site = paste(Longitude, Latitude))


# Species
PA_sp = PA_data %>% select(Species) %>% distinct()

#Species per site
PA_species = PA_data %>% select(Species, Longitude, Latitude, Abundance) %>%
  distinct() %>% mutate(Site = paste(Longitude, Latitude), Abundance = 1) %>% select(-Longitude, -Latitude)

#Env data
env_data = PA_data %>% select(Longitude, Latitude, Bot_depth, SST, SBT, ICC) %>% distinct()

## Create the data frame
combinations = expand.grid(Species = PA_sp$Species, Site = PA_sites$Site)
data_PA = left_join(combinations, PA_species)
data_PA$Abundance[is.na(data_PA$Abundance)] = 0

data_PA = separate(data_PA, Site, into = c("Longitude","Latitude"), sep = " ")
data_PA$Longitude = as.numeric(data_PA$Longitude); data_PA$Latitude = as.numeric(data_PA$Latitude)
data_PA = left_join(data_PA, env_data)
data_PA$ICC = round(data_PA$ICC, 3)

## 3. Presence Absence models calibrations & evaluation ####

models_species_PA = list()
for (k in 1:length(unique(data_PA$Species))) {
  
  data_spe_s = filter(data_PA, Species == unique(data_PA$Species)[k])
  
  if(length(unique(data_spe_s$ICC)) > 5 & max(unique(data_spe_s$ICC), na.rm =T ) > 0.05){
    
  models_species_PA[[k]] = gam(Abundance ~ s(SST, k = 4) + s(Bot_depth, k = 4) + ICC +
                              s(SBT,k=4), data = data_spe_s, family = "binomial")
  }else{
    models_species_PA[[k]] = gam(Abundance ~ s(SST, k = 4) + s(Bot_depth, k = 4) +
                                   s(SBT,k=4), data = data_spe_s, family = "binomial")
  }
  print(paste("This is species", k, "of", length(unique(data_PA$Species))))  
}

## 4. Biomass models calibrations & evaluation ####
AB_data$ICC = round(AB_data$ICC,3)

models_species_AB = list()
for (k in 1:length(unique(AB_data$Species))) {
  tryCatch({
  
  data_spe_s = filter(AB_data, Species == unique(AB_data$Species)[k])
  
  if(length(unique(data_spe_s$ICC)) > 5 & max(unique(data_spe_s$ICC), na.rm =T ) >= 0.05){
    model_AB = gam(log(Biomass) ~ s(SST, k = 4) + s(Bot_depth, k = 4) + ICC +
                     s(SBT,k=4), data = data_spe_s)
  }else{
    model_AB = gam(log(Biomass) ~ s(SST, k = 4) + s(Bot_depth, k = 4) +
                     s(SBT,k= 4), data = data_spe_s)
  }}, error=function(e){})

  models_species_AB[[k]] = model_AB
  
  print(paste("This is species", k, "Barents Sea"))  
}

models = list(models_species_PA, models_species_AB)
saveRDS(models, "models_PA_AB_ICC_AB_Gaussian.rds")

## 5. Diagnostics ####
# I only want to know the diagnostic of the models including 
# the species selected in the analysis.

selected_sp = unique(read.csv("data/data_offshore_regional_realized_shift.csv")$Species)
extract_dev_gam = function(x) {summary(x)$dev.exp}

# PA models
auc_list = list()
for (i in 1:length(models_species_PA)){
  tryCatch({
  model = models_species_PA[[i]]
  data_spe_s = na.omit(filter(data_PA, Species == unique(PA_species$Species)[i]))
  predict = predict.gam(model, newdata = data_spe_s)
  auc_list[[i]] = pROC::auc(data_spe_s$Abundance, predict)
  }, error=function(e){})
}
auc = do.call(rbind, auc_list)

auc_sp = auc[which(PA_sp$Species %in% selected_sp)]
range(auc_sp)
mean(auc_sp)
hist(auc_sp,main = "AUC of presence-absence models", xlab = "AUC") # 26

# Abundance models
dv = unlist(lapply(models_species_AB, FUN = extract_dev_gam))
dv_sp = dv[which(unique(AB_data$Species) %in% selected_sp)]
range(dv_sp)
mean(dv_sp)

hist(dv_sp,main = "% Deviance explained of Abundance models", xlab = "% DV explained") # 26

# Plots

par(mfrow = c(1,2))
hist(auc_sp, main = "Histogram of AUC", xlab = "AUC") # 0.89
hist(dv_sp,
     main = "% Deviance explained of Abundance models", xlab = "% DV explained") # 26

data_models_diag = data.frame(AUC = auc_sp,
                         DV = dv_sp,
                         Species = unique(AB_data$Species)[which(unique(AB_data$Species) %in% selected_sp)])

data_models_diag %>% arrange(desc(DV))


## 6. Regional predictions: Annual thermal envelope projections ####

load("data/environmental_envelope_all.RData")
models_species_PA = readRDS("models_PA_AB_ICC_AB_Gaussian.rds")[[1]]
models_species_AB = readRDS("models_PA_AB_ICC_AB_Gaussian.rds")[[2]]

change_allbar = list()
for (j in 1:length(unique(AB_data$Species))){
  Region_r = "Barents all"
  
  Species = unique(AB_data$Species)[j]
  model_gam_AB = models_species_AB[[which(unique(AB_data$Species) == Species)]]
  model_gam_PA = models_species_PA[[which(unique(data_PA$Species) == Species)]]
  #model_rf = models_fitted[[2]][[which(unique(PA_species$Species) == Species)]]
  species_r = unique(PA_species$Species)[which(unique(PA_species$Species) == Species)]
  
  pred_summer_allbar = list()
  centr_summer_allbar = list()
  for (r in 1:26){ # And create annual predictions ~ 4 min x species
    # I create a data frame with the env data for predicting.
    pred_data = data.frame(raster::coordinates( lay_summer_bar[[r]]),
                           SBT = raster::values(lay_summer_bar[[r]]$SBT),
                           SST = raster::values(lay_summer_bar[[r]]$SST),
                           ICC = round(raster::values(lay_summer_bar[[r]]$ICC),3),
                           Bot_depth = raster::values(lay_summer_bar[[r]]$Bot_depth))
    
    pred_data$ICC[pred_data$ICC == 0] = NA
    pred_data$Bot_depth[pred_data$Bot_depth > 0] = 0
    pred_data$Bot_depth = abs(pred_data$Bot_depth)
    pred_data = na.omit(pred_data)
    pred_data = filter(pred_data, pred_data$Bot_depth > 10)
    #Latitude = coordinates(lay_summer_allbar[[r]])[,1],
    #Longitude = coordinates(lay_summer_allbar[[r]])[,2])
    
    pred = data.frame(pred_data,
                      pred_PA = predict(model_gam_PA, newdata = pred_data, type = "response"),
                      pred_AB = predict.gam(model_gam_AB, newdata = pred_data,type = "response")) #,
    
    
    correction = (sum(exp(residuals(model_gam_AB))))/length(residuals(model_gam_AB))    
    
    #pred_rf = predict(model_rf, newdata = pred_data))
    pred = pred %>% mutate(pred_fin = pred_PA * exp(pred_AB) * correction)
    
    centr_summer_allbar[[r]] = data.frame(Latitude  = weighted.mean(x = pred$y,w = pred$pred_fin, na.rm = T),
                                           Longitude = weighted.mean(x = pred$x,w = pred$pred_fin, na.rm = T), 
                                           Year = r + 1993)
    
    print(r)
  }
  centroids = do.call(rbind, centr_summer_allbar)
  change_allbar[[j]] = data.frame(Species = species_r, Region = Region_r, centroids, Season = "Summer")
  print(paste(species_r,j))
}
allbar_chan = do.call(rbind,change_allbar)
save(allbar_chan, file = "allbar_chan.RData")

change_now = list()
for (j in 1:length(unique(AB_data$Species))){
  Region_r = "Norwegian Sea"
  
  Species = unique(AB_data$Species)[j]
  model_gam_AB = models_species_AB[[which(unique(AB_data$Species) == Species)]]
  model_gam_PA = models_species_PA[[which(unique(data_PA$Species) == Species)]]
  #model_rf = models_fitted[[2]][[which(unique(PA_species$Species) == Species)]]
  species_r = unique(PA_species$Species)[which(unique(PA_species$Species) == Species)]
  
  pred_summer_norwegian = list()
  centr_summer_norwegian = list()
  for (r in 1:26){ # And create annual predictions ~ 4 min x species
    # I create a data frame with the env data for predicting.
    pred_data = data.frame(raster::coordinates(lay_autumn_now[[r]]),
                                   SBT = raster::values(lay_autumn_now[[r]]$SBT),
                                   SST = raster::values(lay_autumn_now[[r]]$SST),
                                   ICC = round(raster::values(lay_autumn_now[[r]]$ICC),3),
                                   Bot_depth = raster::values(lay_autumn_now[[r]]$Bot_depth))
    
    pred_data$ICC[is.na(pred_data$ICC)] = 0
    pred_data$Bot_depth[pred_data$Bot_depth > 0] = 0
    pred_data$Bot_depth = abs(pred_data$Bot_depth)
    pred_data = na.omit(pred_data)
    #Latitude = coordinates(lay_summer_now[[r]])[,1],
    #Longitude = coordinates(lay_summer_now[[r]])[,2])
    
    pred = data.frame(pred_data,
                      pred_PA = predict(model_gam_PA, newdata = pred_data, "response"),
                      pred_AB = predict(model_gam_AB, newdata = pred_data, "response")) #,
    
    correction = (sum(exp(residuals(model_gam_AB)))/length(residuals(model_gam_AB)))    
    
    #pred_rf = predict(model_rf, newdata = pred_data))
    pred = pred %>% mutate(pred_fin = pred_PA * exp(pred_AB)* correction)
    
    centr_summer_norwegian[[r]] = data.frame(Latitude  = weighted.mean(x = pred$y,w = pred$pred_fin),
                                           Longitude = weighted.mean(x = pred$x,w = pred$pred_fin), 
                                           Year = r + 1993)
    print(r)
  }
  centroids = do.call(rbind, centr_summer_norwegian)
  change_now[[j]] = data.frame(Species = species_r, Region = Region_r, centroids, Season = "Autumn")
  print(paste(species_r,j))
}
now_chan = do.call(rbind,change_now)

change_nor = list()
for (j in 1:length(unique(AB_data$Species))){
  Region_r = "North Sea"
  
  Species = unique(AB_data$Species)[j]
  model_gam_AB = models_species_AB[[which(unique(AB_data$Species) == Species)]]
  model_gam_PA = models_species_PA[[which(unique(data_PA$Species) == Species)]]
  #model_rf = models_fitted[[2]][[which(unique(PA_species$Species) == Species)]]
  species_r = unique(PA_species$Species)[which(unique(PA_species$Species) == Species)]
  
  pred_summer_north = list()
  centr_summer_north = list()
  for (r in 1:26){ # And create annual predictions ~ 4 min x species
    # I create a data frame with the env data for predicting.
    pred_data = data.frame(raster::coordinates(lay_summer_nor[[r]]),
                           SBT = raster::values(lay_summer_nor[[r]]$SBT),
                           SST = raster::values(lay_summer_nor[[r]]$SST),
                           ICC = round(raster::values(lay_summer_nor[[r]]$ICC),3),
                           Bot_depth = raster::values(lay_summer_nor[[r]]$Bot_depth))
    
    pred_data$ICC[is.na(pred_data$ICC)] = 0
    pred_data$Bot_depth[pred_data$Bot_depth > 0] = 0
    pred_data$Bot_depth = abs(pred_data$Bot_depth)
    pred_data = na.omit(pred_data)
    
    #Latitude = coordinates(lay_summer_nor[[r]])[,1],
    #Longitude = coordinates(lay_summer_nor[[r]])[,2])
    
    pred = data.frame(pred_data,
                      pred_PA = predict(model_gam_PA, newdata = pred_data, "response"),
                      pred_AB = predict(model_gam_AB, newdata = pred_data, "response")) #,
    
    correction = (sum(exp(residuals(model_gam_AB)))/length(residuals(model_gam_AB)))    
    
    #pred_rf = predict(model_rf, newdata = pred_data))
    pred = pred %>% mutate(pred_fin = pred_PA * exp(pred_AB)* correction)
    
    centr_summer_north[[r]] = data.frame(Latitude  = weighted.mean(x = pred$y,w = pred$pred_fin),
                                           Longitude = weighted.mean(x = pred$x,w = pred$pred_fin), 
                                           Year = r + 1993)
    
    
    print(r)
  }
  centroids = do.call(rbind, centr_summer_north)
  change_nor[[j]] = data.frame(Species = species_r, Region = Region_r, centroids, Season = "Summer")
  print(paste(species_r,j))
}
nor_chan = do.call(rbind,change_nor)

#load("allbar_chan.RData")
thermal_centroid_annual_only_gam = list(allbar_chan,now_chan, nor_chan)
saveRDS(thermal_centroid_annual_only_gam, file = "data/annual_chentroid_thermal_changes_only_gam.rds")

## 7. Calculate the slope of changes ##### 
thermal_centroid_annual_only_gam = readRDS("data/annual_chentroid_thermal_changes_only_gam.rds")

bar_chan = thermal_centroid_annual_only_gam[[1]]
now_chan = thermal_centroid_annual_only_gam[[2]]
nor_chan = thermal_centroid_annual_only_gam[[3]]

# Select species included in the analysis in each region

selected_sp = unique(read.csv("data/data_offshore_regional_realized_shift.csv")[,c("Species", "Region")])

bar_chan = bar_chan %>% filter(Species %in% selected_sp$Species[selected_sp$Region == "Barents Sea"])
now_chan = now_chan %>% filter(Species %in% selected_sp$Species[selected_sp$Region == "Norwegian Sea"])
nor_chan = nor_chan %>% filter(Species %in% selected_sp$Species[selected_sp$Region == "North Sea"])

# Calculate the slope of change in each of them
slopes_bar = list()
for (i in 1:length(unique(bar_chan$Species))){
  dat_s = filter(bar_chan, Species == unique(bar_chan$Species)[i])
  slopes_bar[[i]] = data.frame(lat_s = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,1],
                               lat_d = summary(lm(Latitude ~ Year, data = dat_s))$r.squared,
                               lat_p = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,4],
                               Species = unique(bar_chan$Species)[i]) 
}
slopes_bar = do.call(rbind,slopes_bar)

slopes_now = list()
for (i in 1:length(unique(now_chan$Species))){
  dat_s = filter(now_chan, Species == unique(now_chan$Species)[i])
  slopes_now[[i]] = data.frame(lat_s = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,1],
                               lat_d = summary(lm(Latitude ~ Year, data = dat_s))$r.squared,
                               lat_p = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,4],
                               Species = unique(now_chan$Species)[i]) 
}
slopes_now = do.call(rbind,slopes_now)

slopes_nor = list()
for (i in 1:length(unique(nor_chan$Species))){
  dat_s = filter(nor_chan, Species == unique(nor_chan$Species)[i])
  slopes_nor[[i]] = data.frame(lat_s = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,1],
                               lat_d = summary(lm(Latitude ~ Year, data = dat_s))$r.squared,
                               lat_p = summary(lm(Latitude ~ Year, data = dat_s))$coefficients[2,4],
                               Species = unique(nor_chan$Species)[i]) 
}
slopes_nor = do.call(rbind,slopes_nor)

slopes_bar$Region = "Barents Sea"; slopes_now$Region = "Norwegian Sea"; slopes_nor$Region = "North Sea"

slopes = rbind(slopes_bar, slopes_now, slopes_nor)
slopes$Region = factor(slopes$Region, levels = c("Barents Sea", "Norwegian Sea", "North Sea"))

# Simple visualisation of the overall regional thermal envelope shift rates
boxplot(slopes$lat_s ~ slopes$Region, ylim = c(-0.07,0.07))
abline(h = 0)

# Slopes moving in which direction

a = slopes %>% 
  mutate(agreement = case_when(lat_s > 0 ~ 1,
                               lat_s < 0 ~ 0))

table(a$Region, a$agreement)
#Save the results
save(slopes, file = "data/slopes_thermal_envelope.RData")



## 8. Save the table with both realized and thermal envelope shifts ####


realized_slope = read.csv("data/species_lat_shifts_slopes.csv")
load("data/slopes_thermal_envelope.RData")

colnames(slopes)
colnames(realized_slope)

realized_slope = realized_slope %>% select(Species, n, Region, lat_s, lat_p, lat_dev) %>% 
  rename(lat_r_s = lat_s, lat_r_p = lat_p, lat_r_dev = lat_dev)


slopes_table = left_join(realized_slope, slopes) %>% select(Species, Region, n, lat_r_s, lat_r_dev, lat_r_p, lat_s, lat_d, lat_p)
write.csv(slopes_table, "tables/slopes_table.csv")




