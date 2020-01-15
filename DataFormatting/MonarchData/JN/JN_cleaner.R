#--------------#
#-JourneyNorth-#
#--------------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#
setwd("~/Monarchs/DataFormatting/MonarchData/JN")

#-Load Data-#

JN <- read.csv("JN.csv")

#-Format Data-#
JN$Lat <- as.numeric(as.character(JN$Lat))
JN$Long <- as.numeric(as.character(JN$Long))

JN <- JN %>% 
  filter(!(is.na(Lat)|is.na(Long)))

#-Seasonal Dynamics-#

#Format as date
JN$Date <- as.Date(JN$Date)

#Set temporal periods
JN <- JN %>% 
  mutate(yr = format(Date, "%Y"),
         mo_day = format(Date, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                         ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                                ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
Monarchs <- JN %>% filter(period > 0)

#-Standardize-#
Monarchs <- cbind(Monarchs[,c(3:2,4:6)], rep("JN", length(Monarchs[,1])), 
                  rep(NA, length(Monarchs[,1])),  rep(NA, length(Monarchs[,1])))
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
