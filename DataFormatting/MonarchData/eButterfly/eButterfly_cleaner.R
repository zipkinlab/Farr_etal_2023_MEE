#------------#
#-eButterfly-#
#------------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#

setwd("~/Monarchs/DataFormatting/MonarchData/eButterfly")

#-Load Data-#

eButt <- read.csv("all_records_2019_02_26_10_00_11.csv")

#-Format Data-#

eButt$longitude <- as.numeric(as.character(eButt$Longitude))
eButt$latitude <- as.numeric(as.character(eButt$Latitude))

eButt <- eButt %>% 
  filter(!(is.na(longitude)|is.na(latitude))) %>%
  filter(!(Adult.Unknown == "" & Adult.Male == "" & (Adult.Female == ""|is.na(Adult.Female)))) %>%
  select(latitude, longitude, Date.Observed, Genus, Species)

#Format as date
eButt$Date.Observed <- as.Date(eButt$Date.Observed, format = "%m/%d/%Y")

#-Seasonal Dynamics-#

eButt <- eButt %>% 
  mutate(yr = format(Date.Observed, "%Y"),
         mo_day = format(Date.Observed, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                         ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                                ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
eButt <- eButt %>% filter(period > 0 & yr >= 2016 & yr <=2018)

#-Export butterflies-#
colnames(eButt) <- c("y", "x", "date", "genus", "species", "yr", "mo_day", "period")
save(eButt, file = "eButterflyData.Rdata")

#-Export Monarchs-#
Monarchs <- eButt %>% filter(genus == "Danaus" & species == "plexippus")

#-Standardize-#
Monarchs <- cbind(Monarchs[,c(2:1,6:8)], rep("eButt", length(Monarchs[,1])), 
                  rep(NA, length(Monarchs[,1])), rep(NA, length(Monarchs[,1])))
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
