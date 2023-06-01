### 
### Paper: Species geographical latitudinal shifts in the North - Norwegian - Barents Seas
### Author: Cesc Gordó-Vilaseca
### Date: 25.05.2023
setwd("../Final code")
### 1. Load library and data ############################################

# Library 
library(tidyverse)
library(mgcv)
detach("package:raster", unload=TRUE) # In case it's loaded from previous script


# Data
data_full = read.csv("data/data_ready.csv") # Dataset after 1st filtration
fish_list = read.csv("data/fish_list_johannesen_2021.csv") # Species list allowed in Barents Sea Surveys
traits = read.csv("data/Species_traits.csv")


### 2. Second data filtration: sites ############################################

# I filter the data to only use the species in Johannesen et al. 2019 for the IMR data
# (Norwegian and Barents Sea). Also, I discard winter data for standardisation
# For the Barents Sea, I only keep data after 2004, because it's what the 
# data collectors recommended. I also fitler out coastal data. 

data_off = data_full %>% 
  filter(((Species %in% fish_list$Species & Region != "North Sea") | (Region == "North Sea")),
         Season != "Winter",coast_off_20 == "offshore", 
         ((Region == "Barents Sea" & Year > 2003) | Region != "Barents Sea")) 

ggplot(unique(data_off[,c("Longitude", "Latitude", "Year", "Region")])) +
  geom_point(aes(x = Longitude, y = Latitude, color = Region)) + facet_wrap(.~ Year)

## Data from 2014 and 2016 have no sampling in north Barents Sea, which 
## could bias the results. I discard those years for the Barents Sea, and data on top of Svalbard
## May alter the latitudinal shift distributions because of the continent
## Data east of 10 degrees Longitude in the North Sea is not there in the first years, which could 
## also have an effect
data_off = data_off %>% filter((Region == "Barents Sea"  & (Year!=2014) & (Year != 2016) & 
                                  !(Longitude < 25 & Latitude > 79) &
                                 !(Latitude > 80))|
                                  (Region == "North Sea" & 
                                     !(Longitude > 9 & Latitude > 56)) | 
                                  Region == "Norwegian Sea")

## Spatio temporal distribution of the data
ggplot(unique(data_off[,c("Longitude", "Latitude", "Year", "Region")])) + 
  geom_point(aes(x = Longitude, y = Latitude, col = Region)) + facet_wrap(.~Year)
ggplot(unique(data_off[,c("Longitude", "Latitude", "Year", "Region")])) + geom_boxplot(aes(x = Year, y = Longitude, group = Year)) + 
  facet_wrap(.~ Region)

ggplot(unique(data_off[,c("Longitude", "Latitude", "Year", "Region")])) + geom_boxplot(aes(x = Year, y = Latitude, group = Year)) + 
  facet_wrap(.~ Region)

# To take into account this variability across years, I add one 
# variable of anual mean latitude of the surveys in each point,
# which will be included in the regression

data_all_lat = data_off %>% dplyr::select(Longitude, Latitude, Year, Region) %>% 
  distinct() %>% group_by(Year, Region) %>% 
  summarise(annual_latitude = mean(Latitude),
            annual_longitude = mean(Longitude),
            annual_latitude_median = median(Latitude),
            annual_longitude_median  = median(Longitude)) %>% ungroup()

data_off = left_join(data_off, data_all_lat, by = c("Region", "Year")) 


### 3. Third data filtration: Rare species ###########
## I filter out very rare species, because I do not trust the centroid calculated with very few points

data_off %>% select(Year, Region, Longitude, Latitude) %>% distinct() %>%
  group_by(Region, Year) %>% count() %>% print(n = 200)
# Very little data in the Norwegian Sea in 2001... keep in mind

## I eliminate pecies with less than ry sites per year, and less than yr years 
## of data in each region

ry = 5
yr = 10

sp_no_sel_1 = data_off %>% dplyr::select(Species, Longitude, Latitude, Region, Bot_depth, Year) %>% group_by(Species, Region, Year) %>%
  mutate(sites_in_year = n()) %>% filter(sites_in_year < ry) %>%
  mutate(Reg_sp = paste(Region, Species, Year)) %>% ungroup() %>% dplyr::select(Reg_sp, Species, Region, Year) %>% distinct()

sp_no_sel_2 = data_off %>% mutate(Reg_sp = paste(Region, Species, Year)) %>% filter(!(Reg_sp %in% unique(sp_no_sel_1$Reg_sp))) %>% ungroup() %>% 
  dplyr::select(Region,Species, Year) %>% distinct() %>% group_by(Species, Region) %>%  mutate(years_found = n()) %>% 
  filter( years_found < yr) %>%
  mutate(Reg_sp2 = paste(Region, Species)) %>% ungroup() %>% dplyr::select(Reg_sp2) %>% distinct()

data_off = data_off %>% mutate(Reg_sp = paste(Region, Species, Year),
                               Reg_sp2 = paste(Region, Species)) %>% 
  filter(!(Reg_sp %in% unique(sp_no_sel_1$Reg_sp)) &
           !(Reg_sp2 %in% unique(sp_no_sel_2$Reg_sp2))) %>% ungroup()

# How many species left? 

unique(data_off$Species[data_off$Region == "Barents Sea"]) # 35 sp  
unique(data_off$Species[data_off$Region == "Norwegian Sea"]) # 25 sp
unique(data_off$Species[data_off$Region == "North Sea"]) # 57 sp
unique(data_off$Species) # 82

# Save this data, for example to count the number of records in each species

write.csv(data_off, "data/data_offshore_regional_realized_shift.csv")


### 4. Loop to calculate each species geographic shift
pin_f = list()
data_original = list()
models_lat_regions = list()

for (k in 1:length(unique(data_off$Region))) { 
  Region = unique(data_off$Region)[k]
  data = data_off %>% filter(Region == unique(data_off$Region)[k]) %>%
    select(Bot_depth, Longitude, Latitude, Biomass,
           annual_latitude_median,annual_longitude_median, Year,
           Species, Region, Season, annual_latitude, annual_longitude)
  pin = list()
  data_orig = list()
  
  models_lat = list()
  for (i in 1:length(unique(data$Species))){
    
    Region = unique(data$Region)
    x = unique(data$Species)[i]
    data_sp = data %>% filter(Species == x[1]) %>%
      group_by(Year, Season, annual_latitude, annual_longitude, annual_latitude_median,
               annual_longitude_median,Species, Region) %>%
      mutate(Total_biomass = sum(Biomass)) %>% 
      summarise(mean_lat = weighted.mean(weight = Biomass/Total_biomass, Latitude),
                mean_lon = weighted.mean(weight = Biomass/Total_biomass, Longitude),
                mean_depth = mean(Bot_depth))
    
    model_lat = gam(mean_lat ~ Year + annual_latitude +
                      s(mean_depth, k = 3) + annual_latitude_median, data = data_sp)
    
    
    lm_model_lat = summary(model_lat)
    intercept_lat = lm_model_lat$p.coeff[1]
    lat_s = lm_model_lat$p.coeff[2]
    lat_std_error = lm_model_lat$se[2]
    an_lat_s = lm_model_lat$p.coeff[3]
    an_lat_p = lm_model_lat$p.pv[3]
    lat_p = lm_model_lat$p.pv[2]
    lat_dev = lm_model_lat$dev.expl        
 
    
    pin[[i]] = data.frame(Species = x, n = length(unique(data_sp$Year)), Region = Region, 
                          intercept_lat = intercept_lat, lat_s = lat_s, seas_lat = 0,
                          an_lat_s = an_lat_s, an_lat_p, lat_p = lat_p, 
                          lat_dev = lat_dev,  lat_std_error = lat_std_error)
    data_orig[[i]] = data_sp
    models_lat[[i]] = model_lat
  }
  pin_f[[k]]= do.call(rbind, pin)
  data_original[[k]] = do.call(rbind, data_orig)
  models_lat_regions[[k]] = models_lat
  
}

realized_slope = do.call(rbind, pin_f) # slopes of regional shifts
realized_acent = do.call(rbind, data_original) # regional annual centroids 

realized_slope$lat_s = realized_slope$lat_s * 111 # Degrees of latitude to Km
realized_slope$Region = factor(realized_slope$Region, levels = c("Barents Sea","Norwegian Sea","North Sea"))

## 5. Plotting results ####

# Changes in general 
realized_slope %>% select(Species,Region, lat_s, lat_std_error) %>%
  mutate(change_lat = case_when(lat_s < 0 ~ 0,
                                lat_s > 0 ~ 1)) %>% 
  group_by(change_lat, Region) %>% count()

# Clear species changes
realized_slope %>% filter(lat_p < 0.05) %>% select(Species,Region, lat_s, lat_std_error) %>%
  mutate(change_lat = case_when(lat_s < 0 ~ 0,
                                lat_s > 0 ~ 1)) %>% 
  group_by(change_lat, Region) %>% count()

# Latitudinal shifts different from 0 if species as replicates? 
pdf("figures/Figure_community_shifts.pdf", width = 5, height = 5)
boxplot(lat_s ~ Region, data = realized_slope, outline = F,
        col = c("#FF7F50", "#16E0BD", "#FFF4E9"), xlab = "",
        main = "", ylab = "Latitudinal shift (km/yr)",
        ylim = c(-15,15))

means <- tapply(realized_slope$lat_s, realized_slope$Region, mean)

#add means as circles to each boxplot
points(means, pch=20, cex=1.5, col = "white")


abline(h =0, lty = 5)
text(1, 14, expression("*"), cex = 2) 
text(3, 14, expression("*"), cex = 2) 

text(1, -14, expression("n = 35"), cex = 1.5) 
text(2, -14, expression("n = 25"), cex = 1.5) 
text(3, -14, expression("n = 57"), cex = 1.5) 


dev.off()
# Statistical differences
# Non-parametric test
wilcox.test(realized_slope$lat_s[realized_slope$Region == "Barents Sea"], conf.int = T, correction = F) # Significant
wilcox.test(realized_slope$lat_s[realized_slope$Region == "Norwegian Sea"], conf.int = T,correction = F) # Not significant
wilcox.test(realized_slope$lat_s[realized_slope$Region == "North Sea"], conf.int = T) # Significant

#Parametric tests for comparison with other studies
t.test(realized_slope$lat_s[realized_slope$Region == "Barents Sea"], conf.int = T, correction = F) # Significant
t.test(realized_slope$lat_s[realized_slope$Region == "Norwegian Sea"], conf.int = T,correction = F) # Not significant
t.test(realized_slope$lat_s[realized_slope$Region == "North Sea"], conf.int = T) # Significant


shapiro.test(realized_slope$lat_s[realized_slope$Region == "Barents Sea"]) # Significant
shapiro.test(realized_slope$lat_s[realized_slope$Region == "Norwegian Sea"]) # Not significant
shapiro.test(realized_slope$lat_s[realized_slope$Region == "North Sea"]) # Significant



write.csv(realized_slope, "data/species_lat_shifts_slopes.csv")



## 6. Figures for Supplementary material of each species ####
realized_shifts_sig_lat = left_join(realized_slope, realized_acent) %>% filter(lat_p < 0.05)

pdf("figures/Figure_S2_barents_sea_species_shifts.pdf", width = 6, height = 6)
ggplot(filter(realized_shifts_sig_lat,Region == "Barents Sea"),
       aes(x = Year, y = mean_lat)) + 
  geom_point() +  geom_smooth(method = "lm") + ylab("Mean latitude (°)") +
  facet_wrap(.~ Species, scale = "free") + 
  theme(strip.text = element_text(face = "italic"))
dev.off()


pdf("figures/Figure_S3_norwegian_sea_species_shifts.pdf", width = 4.2, height = 4)
ggplot(filter(realized_shifts_sig_lat,Region == "Norwegian Sea"),
       aes(x = Year, y = mean_lat)) + 
  geom_point() +  geom_smooth(method = "lm") + ylab("Mean latitude (°)") +
  facet_wrap(.~ Species, scale = "free") + 
  theme(strip.text = element_text(face = "italic"))
dev.off()

pdf("figures/Figure_S4_north_sea_species_shifts.pdf", width = 8)
ggplot(filter(realized_shifts_sig_lat,Region == "North Sea"),
       aes(x = Year, y = mean_lat)) + 
  geom_point() +  geom_smooth(method = "lm") + ylab("Mean latitude (°)") +
  facet_wrap(.~ Species, scale = "free") + 
  theme(strip.text = element_text(face = "italic"))
dev.off()

## 7. MAP FIGURE 1 AND FIGURE S1 ####
library(rgdal)
library(tidyverse)
no_balt = rgdal::readOGR("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Envars/Dist_coast/world_countries.shp")
# Correct the extent of the Shapefile
e <- raster::extent(-10,45,50,83)
no_balt = raster::crop(no_balt, e)

base_map = ggplot() + 
  geom_polygon(data = no_balt, aes(x = long, y = lat, group = group), colour = "black", fill = NA) + 
  scale_fill_viridis_c() + theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text =element_text(size=10)) +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0))


surveys_map = distinct(data_off[,c("Longitude", "Latitude", "Year", "Region")])

pdf("Figure_2.pdf", width = 7, height = 6)

base_map + geom_point(data = surveys_map, aes(x = Longitude, y = Latitude, col = Region), size = 0.9) + 
  scale_color_manual(values=c("#ff7f50ff", "#16e0bdff", "khaki3"))+
  ylab("Latitude (°)") + xlab("Longitude (°)")

dev.off()

## Figure S1
pdf("figures/Figure_S1.pdf")


ggplot(unique(data_off[,c("Longitude", "Latitude", "Year", "Region")])) + 
  geom_polygon(data = no_balt, aes(x = long, y = lat, group = group), colour = "black", fill = NA) + 
  geom_point(aes(x = Longitude, y = Latitude, col = Region)) + 
  theme_bw() + 
  facet_wrap(.~Year, ncol = 5) +  
  scale_color_manual(values = c("#FF7F50", "#16E0BD", "khaki3"))

dev.off()

