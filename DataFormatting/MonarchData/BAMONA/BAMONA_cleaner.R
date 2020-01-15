#--------#
#-BAMONA-#
#--------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#

setwd("~/Monarchs/DataFormatting/MonarchData/BAMONA")

#-Load Data-#

files <- list.files(pattern = "bamona_data", full.names = TRUE)

BAMONA <- read.csv(files[1])

for(i in 2:length(files)){
  BAMONA <- rbind(BAMONA, read.csv(files[i]))
}

#-Format Data-#

BAMONA$longitude <- as.numeric(BAMONA$Longitude)
BAMONA$latitude <- as.numeric(BAMONA$Lat.Long)

BAMONA <- BAMONA %>% 
  filter(!(is.na(longitude)|is.na(latitude))) %>%
  filter(Organism.Type == "butterfly") %>%
  select(latitude, longitude, Observation.Date, Scientific.Name)

#-Seasonal Dynamics-#

#Format as date
BAMONA$Observation.Date <- as.Date(BAMONA$Observation.Date, format = "%m/%d/%Y")

#Set temporal periods
BAMONA <- BAMONA %>% 
  mutate(yr = format(Observation.Date, "%Y"),
         mo_day = format(Observation.Date, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                         ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                                ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
BAMONA <- BAMONA %>% filter(period > 0)

#-Export butterflies-#
colnames(BAMONA) <- c("y", "x", "date", "scientific_name", "yr", "mo_day", "period")
save(BAMONA, file = "BAMONAData.Rdata")

#-Export Monarchs-#
Monarchs <- BAMONA %>% filter(scientific_name == "Danaus plexippus")

#-Standardize-#
Monarchs <- cbind(Monarchs[,c(2:1,5:7)], rep("BAMONA", length(Monarchs[,1])), 
                  rep(NA, length(Monarchs[,1])), rep(NA, length(Monarchs[,1])))
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
