# At the end of this script, I should have a dataset fully clean and complete with: 
# 180.851 Species records, with biomass, abundance, traits and env. variables
# Distributed across 15626 trawls (Long, Lat, Year, Month)


## Cleaning the original dataset #####
library(tidyverse)
detach("package:raster")
## Hauls data from IMR
data = read.csv("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Data/total_completed_data.csv") # 645648

data_quality <- data %>% filter(Quality_gear < 3) # 641688
data_quality <- data_quality %>% filter(Quality_haul < 3)#629443
data_quality <- data_quality %>% filter(Total_abund  > 0)#600459
data_quality <- data_quality %>% filter(Total_weight > 0)#561335
data_quality <- data_quality %>% filter(Opening      > 0)#490228
data_quality <- data_quality %>% filter(Bot_depth    > 0)#489572
data_quality <- data_quality %>% filter(Distance     > 0)#488793
#data_quality <- data_quality %>% filter(Year         <1987 | Year > 1988) # 467469
data_quality <- data_quality %>% filter(Year         > 1988)# 472464

## There are several Gears available, and its use has a clear temporal pattern
## After looking at the ammount of data per Gear, and the IMR document about Gear
## I only keep the number 3271, which is the Campelen 1800 amb 20 mm mesh size. 
table(data_quality$Gear)

data_quality <- data_quality %>% filter(Gear == 3271)# 277 991
table(data_quality$Gear, data_quality$Year)
# I filter the Measurement code to keep only rounded weight in kg, gr or tonnes
data_quality <- data_quality %>% 
  filter(Measurement == -1 |Measurement == 1 | Measurement == 5 | Measurement == 6 | Measurement == 9) 

# 277 982
table(data_quality$Measurement_raw) # The -1 means in Kg too

## What about opening, swept area, etc? 
boxplot(data_quality$Opening  ~ data_quality$Year)
boxplot(data_quality$Opening  ~ data_quality$Year, ylim = c(0,100))
boxplot(data_quality$Distance  ~ data_quality$Year) # One clear outlier
boxplot(data_quality$Distance  ~ data_quality$Year, ylim = c(0,5)) 

# Both measures slightly decline with time. 

shrimp_trawl <- data_quality %>% 
  mutate(Distance_in_m = data_quality$Distance *1852,
         swept_area = (Opening * Distance_in_m)/1000000)  %>% 
  filter(Opening < 100, Distance_in_m < 5000) # This eliminates some outliers

# 276 342
## Negative trend of swept_area with time, which needs to be considered
boxplot(shrimp_trawl$swept_area ~ shrimp_trawl$Year)

### Keep only fish species - Select tax groups -----
data_clean <- shrimp_trawl %>% filter(grepl(' ', valid_name), class %in% c("Actinopterygii",
                                                                           "Elasmobranchii", 
                                                                           "Holocephalii",
                                                                           "Myxini",
                                                                           "Petromyzonti")) %>% dplyr::rename(Species = valid_name) %>% # 214 207
  dplyr::select(Species, Year, Month, Longitude, Latitude, Bot_depth, Total_weight, Total_abund, Distance_in_m, swept_area)

## Remaining records --> 216 690
### Correcting detected mistakes 

# Surprising but possible
data_clean$Total_abund[data_clean$Species == "Hippoglossoides platessoides" & data_clean$Total_abund > 2500]

# Impossible, mistake in the data recording, I assume two extra zeroes in abundance
data_clean$Total_abund[data_clean$Species == "Cottunculus microps" & data_clean$Total_abund > 100] =  11

# Surprising but possible
data_clean[data_clean$Species == "Etmopterus spinax" & data_clean$Total_abund > 1000,]

# Impossible, mistake in the data recording, I assume two extra zeroes in abundance
data_clean[data_clean$Species == "Lycodes rossi" & data_clean$Total_abund > 70,] 
data_clean[121611,]$Total_abund = 75; data_clean[121611,]$Total_weight = 3.75 

# Surprising but possible
data_clean[data_clean$Species == "Pleuronectes platessa" & data_clean$Total_abund > 200,]

# Surprising but possible
data_clean[data_clean$Species == "Argentina silus" & data_clean$Total_abund > 3000,] 


# Correcting a name synonym
worms_names = worms::wormsbynames(unique(data_clean$Species))
worms_names %>% filter(status != "accepted")

data_clean$Species[data_clean$Species == "Gobiusculus flavescens"] = "Pomatoschistus flavescens"
data_clean$Species[data_clean$Species == "Phrynorhombus norvegicus"] = "Zeugopterus norvegicus"
data_clean$Species[data_clean$Species == "Cottunculus sadko"] = "Cottunculus microps"
## Some species are divided in length classes, but I am not looking at this, so I want to aggregate per haul

data_clean = data_clean %>% group_by(Bot_depth, Year, swept_area,Species, Longitude, Latitude) %>% 
  mutate(Total_abund = sum(Total_abund),
         Total_weight = sum(Total_weight)) %>% 
  mutate( Abundance = Total_abund/swept_area,
          Biomass = Total_weight/swept_area) %>% dplyr::select(-Total_abund, -Total_weight) 

hist(data_clean$Bot_depth)
# I restrict the data to 600 m because deeper than that the coverage is very poor. 
data_clean = filter(data_clean, Bot_depth < 601); hist(data_clean$Bot_depth)
unique(data_clean$Species) # 181

## 214 380

# Before adding the North Sea data, I want to make sure that I 
# can distinguish both datasets

data_clean$Source = "IMR"

### ADD NORTH SEA ICES DATA
north_sea = read.csv("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Final analysis/Traits analysis/ICES_data/north_sea_data.csv")

colnames(data_clean)
colnames(north_sea)

data_clean = data_clean %>% dplyr::select(-Distance_in_m) %>% mutate()
north_sea = north_sea %>% dplyr::select(colnames(data_clean))

# I add north sea data from ICES --> 214 380 + 187 810 = 402 190 records
data_clean = rbind(data_clean, north_sea) 
sum(is.na(data_clean))
#ALL DATA

orig_data = data_clean %>%
  mutate(Site = paste(Longitude, Latitude, Year, Month, Bot_depth, swept_area))

# The cleaned data without considering Months or homogeneity of trawl consists of 402 190 species records
# Now I want to achieve temporal and spatial homogeneity. To do so, I only need sites, not species. 

sites = orig_data %>% ungroup() %>% dplyr::select(-Species, swept_area, -Abundance, -Biomass) %>% distinct()

# Sites: 31 470

## GEOGRAPHICAL PARTITIONING ####
ggplot(sites) + geom_point(aes(x = Longitude, y = Latitude)) + facet_wrap(.~ Year)

## Here I consiedered 3 areas to guarantee a relatively homogeneous 
# spaito-temporal coverage in each region 

barents_sea_data = sites %>% filter(Latitude > 69) %>% mutate(Region = "Barents Sea")
ggplot(barents_sea_data) + geom_point(aes(x = Longitude, y = Latitude), size = 1.2) + facet_wrap(.~ Year)

norwegian_sea_data = sites %>% filter(Latitude <= 69 & Latitude > 62 & Longitude < 20 & Year > 1996) %>% mutate(Region = "Norwegian Sea")
ggplot(norwegian_sea_data) + geom_point(aes(x = Longitude, y = Latitude), size = 1.2) + facet_wrap(.~ Year)

north_data = sites %>% filter(Latitude <= 62) %>% mutate(Region = "North Sea")
ggplot(north_data) + 
  geom_point(aes(x = Longitude, y = Latitude, col = Source)) + facet_wrap(.~ Year)

# Considering both sources of data is a bit tricky because I will need to analyse it separately
# I will mostly use DATRAS data for the North Sea, and include IMR data in the SM

all_data = rbind(barents_sea_data, norwegian_sea_data, north_data)

dim(unique(all_data)) # 31 437
summary(all_data)
# Monthly distribution 
table(all_data$Year,all_data$Month,  all_data$Region)

# Here I could probably use 2008-2018, because of the different seasons it was sampled
table(all_data$Month[all_data$Source == "IMR" & all_data$Region == "North Sea"],
      all_data$Year[all_data$Source == "IMR" & all_data$Region == "North Sea"])

# Here everything looks quite ok. Even though July was only sampled after 2000, August was from the 
# beggining, and July is not a dominant month of sampling

table(all_data$Month[all_data$Source == "ICES-DATRAS" & all_data$Region == "North Sea"],
      all_data$Year[all_data$Source == "ICES-DATRAS" & all_data$Region == "North Sea"])

# North Sea is clearly okay, with January and February months
# Norwegian Sea is also okay if I take October and November

# Barents Sea is the problem --> There are some hauls on February and some on August-September-October
# If this has a spatial pattern, then it may be a problem, but otherwise i could take all of them, I think

ggplot(filter(all_data, Region == "Barents Sea" & Month %in% c(2,7,8,9))) +
  geom_point(aes(x = Longitude, y = Latitude, col = as.factor(Month))) + facet_wrap(.~ Year)


## I see that spatiotemporal variation in Month is clear. I do everything again with seasons


all_data = rbind(barents_sea_data, norwegian_sea_data, north_data) %>% 
  mutate(Season = case_when(Month %in% c(1,2,3)  ~ "Winter",
                            Month %in% c(4,5,6)  ~ "Spring",
                            Month %in% c(7,8,9) ~ "Summer",
                            Month %in% c(10,11,12) ~ "Autumn"))


ggplot(filter(all_data, Season == "Summer")) + geom_point(aes(x = Longitude, y = Latitude)) + facet_wrap(.~ Year)

dim(unique(all_data)) # 31 437

# Monthly distribution 
table(all_data$Year,all_data$Season,  all_data$Region)



# Again, Svalbard easy with Summer
# Norwegian Sea easy with Autumn 
# North Sea easy with Winter
# Barents Sea not easy at all. Either Winter or Summer . 
#I check which one is better in terms of standardisation 

ggplot(filter(all_data, Region == "Barents Sea" & Season == "Winter")) +
  geom_point(aes(x = Longitude, y = Latitude)) + facet_wrap(.~ Year)

ggplot(filter(all_data, Region == "Barents Sea" & Season == "Winter")) +
  geom_boxplot(aes(x = as.factor(Year), y = Latitude))

ggplot(filter(all_data, Region == "Barents Sea" & Season == "Summer")) +
  geom_point(aes(x = Longitude, y = Latitude)) + facet_wrap(.~ Year)

ggplot(filter(all_data, Region == "Barents Sea" & Season == "Summer")) +
  geom_boxplot(aes(x = as.factor(Year), y = Latitude))

# I could also include Winter in the analysis, but it is more temporally biased

# This is important to take into account whenever I want to include Winter in the analysis. 


all_data = all_data %>% filter((Region == "North Sea" & Season == "Summer" & Source == "ICES-DATRAS") |
                                 (Region == "Norwegian Sea" & Season == "Autumn") |
                                 (Region == "Barents Sea" & Season == "Summer"))



ggplot(all_data) +
  geom_point(aes(x = Longitude, y = Latitude, col = Region)) + facet_wrap(.~ Year)
unique(all_data$Region)
# 14 478

## Environmental variables addition
library(raster)

# This layer comes from Marine copernicus
SST <- raster::stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Envars/SST/SST_mean_area.tif")
#This layer is self made
dis <- stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Envars/Dist_coast/distance_to_coast.tif")
# These layers are from Bio Oracle
layers <- stack("C:/Users/06061016/Documents/Treball/PhD/IMR1980-2017/Envars/Bio-Oracle/layers.grd")

names(layers) = c("Bathymetry","Salinity","Currents","DO"," Sea.Ice","Iron","Nitrate","Prod.Bottom","Prod.Surface", "Temp.Bottom", "Surface.Temp","Temp.Range")

plot(dis)
bat = crop(layers[[1]], SST[[1]])
plot(bat)
coordinates = coordinates(SST[[1]])

SBT_d = na.omit(data.frame(raster::extract(layers, coordinates),
                           coordinates, Year = 2020, distance_to_coast = raster::extract(dis, coordinates)))

SBT_d = SBT_d %>% dplyr::rename(Longitude = x, Latitude = y)
SBT_d$Bathymetry[SBT_d$Bathymetry>0] = 0


head(SBT_d)
range(SBT_d$Bathymetry)

df_final = data.frame(all_data, 
                      raster::extract(layers, all_data[,c(3,4)]),
                      raster::extract(dis, all_data[,c(3,4)]))

colnames(df_final)[23] = "distance_to_coast"
df_final$distance_to_coast[is.na(df_final$distance_to_coast)] = 0
# Create 3 different columns with 3 different distances to coast (5,10 and 20 km)

df_f = df_final %>% 
  dplyr::select(Longitude, Latitude, Year, Bot_depth, swept_area, Region, Month, Season, distance_to_coast) %>%
  mutate(Site = paste(Longitude, Latitude, Year, Month, Bot_depth, swept_area),
         coast_off_5 = case_when(distance_to_coast <= 5000 ~ "coastal",
                                 distance_to_coast > 5000 ~ "offshore",
                                 is.na(distance_to_coast) ~ "coastal"),
         coast_off_10 = case_when(distance_to_coast <= 10000 ~ "coastal",
                                  distance_to_coast > 10000 ~ "offshore",
                                  is.na(distance_to_coast) ~ "coastal"),
         coast_off_20 =  case_when(distance_to_coast <= 20000 ~ "coastal",
                                   distance_to_coast > 20000 ~ "offshore",
                                   is.na(distance_to_coast) ~ "coastal")) %>% distinct()


unique(df_f$Region)
ggplot(df_f) + geom_point(aes(x = Longitude, y = Latitude, col = Region)) + facet_grid(.~ coast_off_20)

## How many years have a good coverage of the data? 

ggplot(filter(df_f, coast_off_20 == "coastal")) + 
  geom_point(aes(x = Longitude, y = Latitude, col = Region)) + facet_wrap(.~ Year)


nrow(df_f)# 14 478
colSums(is.na(df_f))
# Now I merge the environmental data in df_f with the 
# biological data in orig_data. 

data_fin = left_join(df_f, orig_data) # 181 685

data_fin = data_fin %>% filter(!is.na(Region))
unique(data_fin$Species) # 221 species

# Number of species from the IMR database after cleaning
# 146
data_fin %>% filter(Region != "North Sea") %>% ungroup() %>% dplyr::select(Species) %>% distinct()

# Number of species from the ICES database after cleaning
# 156
data_fin %>% filter(Region == "North Sea") %>% ungroup() %>% dplyr::select(Species) %>% distinct()

## Trait Addition ####

traits_db = read.csv("data/Species_traits.csv", header = T)
data_traits = left_join(data_fin, traits_db, by = "Species") 

unique(data_traits$Species) # 221 Species

# 178 669 records

# For now I remove the duplicates and only take a single value for each longitude and latitude and year.
# There are 15 sites that have different data for different depths and different swept_area, I remove them

all_data %>% mutate(Longitude = num(Longitude, digits = 3),
                    Latitude = num(Latitude, digits = 3)) %>%
  group_by(Longitude, Latitude, Year) %>% dplyr::count() %>% filter(n > 1)

sites_lly = all_data %>% mutate(Site_sim = paste(Longitude, Latitude, Year)) %>% group_by(Site_sim) %>% dplyr::count() %>% filter(n>1)

data_traits = data_traits %>% mutate(Site_sim = paste(Longitude, Latitude, Year)) %>%
  filter(!Site_sim %in% sites_lly$Site_sim)


colnames(data_traits)
nrow(data_traits) # 181 519
length(unique(data_traits$Site_sim)) # 14 459
# 181 519 records in 14 459 hauls
write.csv(data_traits, "data/data_ready.csv")

