

## 1. Load thermal envelope and geographic shifts data #### 

library(nlme)
library(effects)
library(tidyverse)
library(effects)
detach("package:raster", unload=TRUE) # In case it's loaded from previous script

realized_slope = read.csv("../data/species_lat_shifts_slopes.csv")
load("../data/slopes_thermal_envelope.RData")

slopes_realized = realized_slope %>% #filter(lat_p < 0.05) %>%
  dplyr::select(Species, lat_s, Region) %>%
  rename(real_lat_s = lat_s)

slopes = left_join(slopes_realized, slopes) 

# Add a variable for species records
rec_data = read.csv("../data/data_offshore_regional_realized_shift.csv")
rec_data = rec_data %>% dplyr::select(Longitude, Latitude, Year, Region, Species) %>% 
  distinct() %>% group_by(Region, Species) %>%
  summarise(records = n())

slopes = left_join(slopes, rec_data)

## 2. Create the correlation plot in each region ####

# Transform degrees of latitude to km, the real one is already changed
slopes$thermal_lat_s = slopes$thermal_lat_s * 111

# Shorten species names for better graphic representatin
slopes$Species_l = paste(substr(slopes$Species, 0, 1), ".", substr(slopes$Species, 2, 1),sub(".*? ", "", slopes$Species), sep = "")
# I create numbers to make raeding easier: 
code = data.frame(Species = unique(slopes$Species), ID = as.character(c(1:length(unique(slopes$Species)))))
slopes = left_join(slopes, code)
# Create each regional plot

a = ggplot(filter(slopes, Region == "Barents Sea")) + 
  ggtitle("Barents Sea + Svalbard") + 
  xlab("Thermal niche latitudinal shift (km/yr)") + ylab("Realized latitudinal shift") + 
  geom_point(aes(x = thermal_lat_s, y = real_lat_s), alpha = 0.6) + theme_bw() + 
  geom_abline(slope = 1, col = "red", linetype = "dashed") + 
  #scale_color_gradient2(low = "khaki", mid = "springgreen4", high = "navyblue", midpoint= range(bar_f$Var)[2]/2,  space ="Lab") + 
  geom_text(aes(x = thermal_lat_s, y = real_lat_s, label = Species_l),alpha = 0.6, size = 3) 

b = ggplot(filter(slopes, Region == "Norwegian Sea")) +
  ggtitle("Norwegian Sea") + xlab("Thermal niche latitudinal shift (km/yr)") +
  ylab("Realized latitudinal shift") + 
  geom_point(aes(x = thermal_lat_s, y = real_lat_s), alpha = 0.6) + theme_bw() +
  geom_abline(slope = 1, col = "red", linetype = "dashed") + 
  #scale_color_gradient2(low = "khaki", mid = "springgreen4", high = "navyblue", midpoint= range(bar_f$Var)[2]/2,  space ="Lab") + 
  geom_text(aes(x = thermal_lat_s, y = real_lat_s, label = Species_l),alpha = 0.6, size = 3) 

c = ggplot(filter(slopes, Region == "North Sea")) + 
  ggtitle("North Sea") + xlab("Thermal niche latitudinal shift (km/yr)") + 
  ylab("Realized latitudinal shift") + 
  geom_point(aes(x = thermal_lat_s, y = real_lat_s), alpha = 0.6) + theme_bw() + 
  geom_abline(slope = 1, col = "red", linetype = "dashed") + 
  #scale_color_gradient2(low = "khaki", mid = "springgreen4", high = "navyblue", midpoint= range(bar_f$Var)[2]/2,  space ="Lab") + 
  geom_text(aes(x = thermal_lat_s, y = real_lat_s, label = Species_l),alpha = 0.6, size = 3) 

# Arrange the plots ina  pdf file
pdf("../figures/Figure_4_corr_lat_therm.pdf", width = 4, height = 8)
gridExtra::grid.arrange(a,b,c, nrow = 3)
dev.off()

## 3. Statistical analysis of the correlations between shifts ####
slopes_bar = filter(slopes, slopes$Region == "Barents Sea")
slopes_now = filter(slopes, slopes$Region == "Norwegian Sea")
slopes_nor = filter(slopes, slopes$Region == "North Sea")

# Regional Correlations
cor.test(slopes_bar$thermal_lat_s, slopes_bar$real_lat_s) # Not significant   
cor.test(slopes_now$thermal_lat_s, slopes_now$real_lat_s) # Not significant
cor.test(slopes_nor$thermal_lat_s, slopes_nor$real_lat_s) # Not significant 

# Regional agreement in direction of change
slopes %>% group_by(Region) %>% 
  mutate(agreement = case_when(real_lat_s > 0 & thermal_lat_s > 0 ~ "positive agreement",
                               real_lat_s < 0 & thermal_lat_s < 0 ~ "negative agreement",
                               real_lat_s > 0 & thermal_lat_s < 0 ~ "disagreement pos neg",
                               real_lat_s < 0 & thermal_lat_s > 0 ~ "disgreement neg pos")) %>%
  group_by(Region, agreement) %>% count() %>%
  mutate(percent = case_when(Region == "Barents Sea" ~ n/35 * 100,
                             Region == "Norwegian Sea" ~ n/25 * 100,
                             Region == "North Sea" ~ n/57 *100))


# Regional analysis adding one species at a time, in order of abundance

slopes_bar = slopes_bar %>% arrange(desc(records))
slopes_now = slopes_now %>% arrange(desc(records))
slopes_nor = slopes_nor %>% arrange(desc(records))


res_bar = list()
for (i in 3:nrow(slopes_bar)){
  
  test = cor.test(slopes_bar$thermal_lat_s[1:i], slopes_bar$real_lat_s[1:i])
  res_bar[[i]] = data.frame(pval = test$p.value, cor = test$estimate, n = i,
                            region = "Barents")
  print(i)
  
  
}
res_bar  = do.call(rbind, res_bar)

res_now = list()
for (i in 3:nrow(slopes_now)){
  
  test = cor.test(slopes_now$thermal_lat_s[1:i], slopes_now$real_lat_s[1:i])
  res_now[[i]] = data.frame(pval = test$p.value, cor = test$estimate, n = i,
                            region = "Norwegian Sea")
  print(i)
  
  
}
res_now  = do.call(rbind, res_now)

res_nor = list()
for (i in 3:nrow(slopes_nor)){
  
  test = cor.test(slopes_nor$thermal_lat_s[1:i], slopes_nor$real_lat_s[1:i])
  res_nor[[i]] = data.frame(pval = test$p.value, cor = test$estimate, n = i,
                            region = "North Sea")
  print(i)
  
  
}
res_nor  = do.call(rbind, res_nor)

final_cor_tests = rbind(res_bar, res_now, res_nor)

pears = ggplot(data = final_cor_tests,aes(x = n, y = cor, col = region)) +
  geom_line() + xlab("") + scale_color_discrete(name="") +
  ylab("Pearson correlation")

pval = ggplot(data = final_cor_tests,aes(x = n, y = pval, col = region)) +
  geom_line() + geom_hline(yintercept = 0.05, linetype = "dashed") +
  xlab("Number of species") + scale_color_discrete(name="") + ylab("p-value")


pdf("../figures/Figure_S6_correlations_cumulative.pdf", width = 7, height = 5)
gridExtra::grid.arrange(pears,pval)
dev.off()


####
## 4. Trait analysis ####

# Load data, add traits and records per species

traits = read.csv("../data/Species_traits.csv")
rec_data = read.csv("../data/data_offshore_regional_realized_shift.csv")
rec_data = rec_data %>% select(Longitude, Latitude, Year, Region, Species) %>% 
  distinct() %>% group_by(Region, Species) %>%
  summarise(records = n())

slopes_traits = left_join(na.omit(slopes), traits)
slopes_traits = left_join(slopes_traits, rec_data)
rm(rec_data)
# Save 10 species most and less well sampled
slopes_traits %>% select(Species, Region, records) %>% arrange(records)


## 4.1 Trait analysis: Effect of records ####

# Effect of records in absolute shifting rates in each region

summary(lm(abs(real_lat_s) ~ records, data = slopes_traits[slopes_traits$Region == "Barents Sea",])) # Significant
summary(lm(abs(real_lat_s) ~ records, data = slopes_traits[slopes_traits$Region == "Norwegian Sea",])) # Significant
summary(lm(abs(real_lat_s) ~ records, data = slopes_traits[slopes_traits$Region == "North Sea",])) # Significant
# Significant in all study areas

# Effect of records in non-absolute shifting rate

summary(lm(real_lat_s ~ records, data = slopes_traits[slopes_traits$Region == "Barents Sea",])) 
summary(lm(real_lat_s ~ records, data = slopes_traits[slopes_traits$Region == "Norwegian Sea",])) 
summary(lm(real_lat_s ~ records, data = slopes_traits[slopes_traits$Region == "North Sea",])) 

# Only significant in the Barents Sea --> What about heteroscedasticity
# Heteroscedasticity is clearly violated
plot(lm(real_lat_s ~ records, data = slopes_traits[slopes_traits$Region == "Barents Sea",])) 

# But if taken into account: 

mod1 = gls(real_lat_s ~ records, data = slopes_traits[slopes_traits$Region == "Barents Sea",]) 

mod2 = gls(real_lat_s ~ records, weights = varPower(form = ~ records),
           data = slopes_traits[slopes_traits$Region == "Barents Sea",])

# standardized residuals versus fitted values
plot(mod1, resid(., type = "p") ~ fitted(.), abline = 0)
plot(mod2, resid(., type = "p") ~ fitted(.), abline = 0)
summary(mod2)
cor(slopes_traits[slopes_traits$Region == "Barents Sea",]$real_lat_s,predict(mod2))^2
#  Not corrected but still significant

## 4.2 Trait analysis: Effect of other traits ####
# Barents
# In the Barents Sea I have to add the records as explanatory variable, and 
# This introduces a problem with the residuals, because there is heteroscedasticity
# To correct for this, I apply gls() with a Power structure of the residuals


bar_trait = slopes_traits[slopes_traits$Region == "Barents Sea",]

summary(gls(real_lat_s ~ TempPrefMean + records, weights = varPower(form = ~ records + TempPrefMean), data = bar_trait)) #sig
summary(gls(real_lat_s ~ TempRange + records, weights = varPower(form = ~ records + TempRange), data = bar_trait)) 
summary(gls(real_lat_s ~ DepthMax + records, weights = varPower(form = ~ records), data = bar_trait))
summary(gls(real_lat_s ~ Troph + records, weights = varPower(form = ~ records + Troph), data = bar_trait)) # Almost
summary(gls(real_lat_s ~ DemersPelag + records, weights = varPower(form = ~ records), data = bar_trait)) # sig
summary(gls(real_lat_s ~ AgeMatMin + records, weights = varPower(form = ~ records), data = na.omit(bar_trait[,c("real_lat_s", "records", "AgeMatMin")])))
summary(gls(real_lat_s ~ Fecundity + records, weights = varPower(form = ~ records), data = na.omit(bar_trait[,c("real_lat_s", "records", "Fecundity")])))
summary(gls(real_lat_s ~ Zoo_final + records, weights = varPower(form = ~ records), data = bar_trait)) # significant

## In the Barents Sea, significant traits are TempPrefMean, DemersPelag and Zoo_final
## If considered together, only Zoo_final, the biogeographic trait, resists

summary(gls(real_lat_s ~ TempPrefMean + Zoo_final+ records, weights = varPower(form = ~ records), data = bar_trait))
summary(gls(real_lat_s ~ TempRange + Zoo_final+ records, weights = varPower(form = ~ records), data = bar_trait))
summary(gls(real_lat_s ~ Troph + Zoo_final+ records, weights = varPower(form = ~ records), data = bar_trait))
summary(gls(real_lat_s ~ DemersPelag + Zoo_final+ records, weights = varPower(form = ~ records), data = bar_trait))

## North Sea
nor_trait = slopes_traits[slopes_traits$Region == "North Sea",]

summary(lm(real_lat_s ~ TempPrefMean, data = nor_trait))
summary(lm(real_lat_s ~ DepthMax, data = nor_trait)) # Significant
summary(lm(real_lat_s ~ DepthMax, data = filter(nor_trait, nor_trait$DepthMax <2000))) # Significant

summary(gls(real_lat_s ~ DepthMax, data = nor_trait, weights = varPower(form = ~ DepthMax)))
plot(gls(real_lat_s ~ DepthMax, data = nor_trait, weights = varPower(form = ~ DepthMax)))

sjPlot::plot_model(gls(real_lat_s ~ DepthMax, data = nor_trait, weights = varPower(form = ~ DepthMax)), "pred", show.data = T)

summary(lm(real_lat_s ~ TempRange, data = nor_trait))
summary(lm(real_lat_s ~ MaxLengthTL, data = nor_trait))
summary(lm(real_lat_s ~ Troph, data = nor_trait))  
summary(lm(real_lat_s ~ Fecundity, data = nor_trait))
summary(lm(real_lat_s ~ AgeMatMin, data = nor_trait)) # Significant
summary(lm(real_lat_s ~ AgeMatMin, data = nor_trait[nor_trait$AgeMatMin < 8,])) # significant

sjPlot::plot_model(lm(real_lat_s ~ AgeMatMin, data = nor_trait), "pred", show.data = T)
summary(gls(real_lat_s ~ AgeMatMin, data = nor_trait[!is.na(nor_trait$AgeMatMin),], weights = varPower(form = ~ AgeMatMin)))
plot(gls(real_lat_s ~ AgeMatMin, data = nor_trait[!is.na(nor_trait$AgeMatMin),], weights = varPower(form = ~ AgeMatMin)))


summary(lm(real_lat_s ~ DemersPelag, data = nor_trait))
summary(lm(real_lat_s ~ Zoo_final, data = nor_trait)) 

#Norwegian Sea
now_trait = slopes_traits[slopes_traits$Region == "Norwegian Sea",]

summary(lm(real_lat_s ~ TempPrefMean, data = now_trait)) # Close to significance
summary(lm(real_lat_s ~ DepthMax, data = now_trait)) 
summary(lm(real_lat_s ~ TempRange, data = now_trait))
summary(lm(real_lat_s ~ MaxLengthTL, data = now_trait))
summary(lm(real_lat_s ~ Troph, data = now_trait))  
summary(lm(real_lat_s ~ Fecundity, data = now_trait))
summary(lm(real_lat_s ~ AgeMatMin, data = now_trait)) 
summary(lm(real_lat_s ~ DemersPelag, data = now_trait))
summary(lm(real_lat_s ~ Zoo_final, data = now_trait))

sjPlot::plot_model(lm(real_lat_s ~ Zoo_final, data = now_trait), 
                   type = "pred", terms = "Zoo_final")


## 4.3 Plots ####


# Barents 

model0 = gls(real_lat_s ~ TempPrefMean + records, weights = varPower(form = ~ records + TempPrefMean), data = bar_trait)
model1 = gls(real_lat_s ~ Zoo_final + records, weights = varPower(form = ~ records), data = na.omit(bar_trait[,c("real_lat_s", "records", "Zoo_final")]))
model2 = gls(real_lat_s ~ DemersPelag + records, weights = varPower(form = ~ records), data = na.omit(bar_trait[,c("real_lat_s", "records", "DemersPelag")]))
model3 = gls(real_lat_s ~ Troph + records, weights = varPower(form = ~ records), data = na.omit(bar_trait[,c("real_lat_s", "records", "Troph")]))

summary(model0)
summary(model1)
summary(model2)
summary(model3)

plot(predictorEffects(model0)) 
plot(predictorEffects(model1)) 
plot(predictorEffects(model2)) 
plot(predictorEffects(model3)) 

cor(slopes_traits[slopes_traits$Region == "Barents Sea",]$real_lat_s,predict(model0))^2
cor(slopes_traits[slopes_traits$Region == "Barents Sea",]$real_lat_s,predict(model1))^2
cor(slopes_traits[slopes_traits$Region == "Barents Sea",]$real_lat_s,predict(model2))^2
cor(slopes_traits[slopes_traits$Region == "Barents Sea",]$real_lat_s,predict(model3))^2


# Norwegian Sea 

model4 = gls(real_lat_s ~ TempPrefMean, data = now_trait)
model5 = gls(real_lat_s ~ DemersPelag,data = now_trait)
model6 = gls(real_lat_s ~ Zoo_final,data = now_trait)

summary(model4)
summary(model5)
summary(model6)

# North Sea 

model8 = lm(real_lat_s ~ AgeMatMin, data = nor_trait)
model9 = gls(real_lat_s ~ DepthMax, data = nor_trait)

summary(model8)
summary(model9)

a = effects::predictorEffects(model0)$TempPrefMean 
b = effects::predictorEffects(model1)$Zoo_final
c = effects::predictorEffects(model2)$DemersPelag
d = effects::predictorEffects(model3)$Troph

e = effects::predictorEffects(model4)$TempPrefMean
f = effects::predictorEffects(model5)$DemersPelag
g = effects::predictorEffects(model6)$Zoo_final

i = effects::predictorEffects(model8)$AgeMatMin
j = effects::predictorEffects(model9)$DepthMax


pdf("../figures/Figure_5_traits.pdf", width = 10, height = 9)
gridExtra::grid.arrange(plot(a,main = "Barents Sea" ,xlab = "Preferred temperature (°C)", ylab = "Latitudinal shift (km/yr)", col.line = "darkred", band.colors = "darkred", band.transparency = 0.2), 
                        plot(b,main = "Barents Sea" ,xlab = "Biogeography", ylab = "Latitudinal shift (km/yr)", col.line = "darkred", CI.color = "darkred", band.transparency = 0.2), 
                        plot(c,main = "Barents Sea" ,xlab = "Habitat", ylab = "Latitudinal shift (km/yr)", col.line = "darkred", confint = list(color = "darkred"), band.transparency = 0.2),
                        plot(d,main = "Barents Sea" ,xlab = "Trophic level", ylab = "Latitudinal shift (km/yr)", col.line = "darkred", band.colors = "darkred",band.transparency = 0.2),
                        
                        
                        plot(e, main = "Norwegian Sea",xlab = "Preferred temperature (°C)",ylab = "Latitudinal shift (km/yr)", col.line = "gold3", band.colors = "gold3", band.transparency = 0.2), 
                        plot(f, main = "Norwegian Sea",xlab = "Habitat",ylab = "Latitudinal shift (km/yr)", col.line = "gold3", band.colors = "gold3", band.transparency = 0.2), 
                        plot(g, main = "Norwegian Sea",xlab = "Biogeography",ylab = "Latitudinal shift (km/yr)", col.line = "gold3", band.colors = "gold3", band.transparency = 0.2), 

                        plot(i,main = "North Sea",xlab ="Minimum age at maturity (yr)" , ylab = "Latitudinal shift (km/yr)", col.line = "aquamarine3", band.colors = "aquamarine3", band.transparency = 0.2),
                        plot(j,main = "North Sea",xlab ="Maximum depth (m)" , ylab = "Latitudinal shift (km/yr)", col.line = "aquamarine3", band.colors = "aquamarine3", band.transparency = 0.2),
                        
                        nrow = 3)

dev.off()
