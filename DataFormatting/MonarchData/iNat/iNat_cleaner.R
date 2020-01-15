#-------------#
#-iNaturalist-#
#-------------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#

setwd("~/Monarchs/DataFormatting/MonarchData/iNat")

#-Load Data-#

files <- list.files(pattern = "iNat20", full.names = TRUE)

iNat <- read.csv(files[1])

for(i in 2:length(files)){
iNat <- rbind(iNat, read.csv(files[i]))
}

#-Format Data-#

iNat$longitude <- as.numeric(iNat$longitude)
iNat$latitude <- as.numeric(iNat$latitude)

iNat <- iNat %>% 
  filter(!(is.na(longitude)|is.na(latitude))) %>%
  filter(quality_grade == "research") %>%
  select(latitude, longitude, observed_on, scientific_name)

#-Seasonal Dynamics-#

#Format as date
iNat$observed_on <- as.Date(iNat$observed_on)

#Set temporal periods
iNat <- iNat %>% 
  mutate(yr = format(observed_on, "%Y"),
         mo_day = format(observed_on, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                  ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                  ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
iNat <- iNat %>% filter(period > 0)

#-Export butterflies-#
colnames(iNat) <- c("y", "x", "date", "scientific_name", "yr", "mo_day", "period")
save(iNat, file = "iNatData.Rdata")

#-Export Monarchs-#
Monarchs <- iNat %>% filter(scientific_name == "Danaus plexippus")

#-Standardize-#
Monarchs <- cbind(Monarchs[,c(2:1,5:7)], rep("iNat", length(Monarchs[,1])), 
                  rep(NA, length(Monarchs[,1])), rep(NA, length(Monarchs[,1])))
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
