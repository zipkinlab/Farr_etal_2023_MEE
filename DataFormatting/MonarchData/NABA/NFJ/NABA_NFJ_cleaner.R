#----------------------#
#-NABA: Fourth of July-#
#----------------------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#

setwd("~/Monarchs/DataFormatting/MonarchData/NABA/NFJ")

#-Load Data-#

NFJ1 <- read.csv("NFJ.csv")

NFJ2 <- tbl_df(read.table('Pollard_NABA.txt',sep="\t",header=T,quote="\"",na.strings=c('NA',''),
                             colClasses=c('character','numeric',rep('character',3),rep('numeric',2),rep('character',3),rep('numeric',4),rep('character',2)))) %>%
  filter(Program == "NABA")

#-Format Data-#

NFJ1$date <- as.Date(paste(NFJ1$ObsYear, "-", NFJ1$ObsMonth, "-", NFJ1$ObsDay, sep = ""))
NFJ2$date <- as.Date(NFJ2$SurveyDate)

NFJ1$x <- as.numeric(NFJ1$Lng)
NFJ1$y <- as.numeric(NFJ1$Lat)
NFJ2$x <- as.numeric(NFJ2$Longitude)
NFJ2$y <- as.numeric(NFJ2$Latitude)

NFJ1 <- NFJ1 %>% mutate(effort = Party_Hours, count = Monarch) %>%
  select(x, y, date, count, effort)

NFJ2 <- NFJ2 %>% mutate(effort = Duration.hr, count = Monarchs) %>%
  select(x, y, date, count, effort)

NFJ <- full_join(NFJ1, NFJ2, by = c("x", "y", "date", "count", "effort"))

#-Seasonal Dynamics-#

#Set temporal periods
NFJ <- NFJ %>% 
  mutate(yr = format(date, "%Y"),
         mo_day = format(date, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                         ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                                ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
NFJ <- NFJ %>% filter(period > 0 & yr >= 2016)

#-Spatial information-#

# #Observation circle radius in km (7.5 mi radius)
# radius <- 7.5 * 1.609344
# 
# #Area of observation circle
# Area <- pi*radius^2
# 
# NFJ$Area <- Area

#-Standardize-#

Monarchs <- cbind(NFJ[,c(1:2,6:8)], rep("NFJ", length(NFJ[,1])), NFJ[,c(4:5)])
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
