#---------#
#-Pollard-#
#---------#

#-Libraries-#

library(dplyr)
library(tidyr)

#-Set WD-#

setwd("~/Monarchs/DataFormatting/MonarchData/Pollard")

#-Load Data-#

Pollard1 <- tbl_df(read.table('Pollard_NABA.txt',sep="\t",header=T,quote="\"",na.strings=c('NA',''),
                      colClasses=c('character','numeric',rep('character',3),rep('numeric',2),rep('character',3),rep('numeric',4),rep('character',2)))) %>%
           filter(Program != "NABA")

Pollard2 <- read.csv("PollardData.csv")

#-Format Data-#

Pollard1$x <- as.numeric(Pollard1$Longitude)
Pollard1$y <- as.numeric(Pollard1$Latitude)
Pollard2$x <- as.numeric(Pollard2$Long)
Pollard2$y <- as.numeric(Pollard2$Lat)

Pollard1$date <- as.Date(Pollard1$SurveyDate)
Pollard2$date <- as.Date(Pollard2$Date, format = "%m/%d/%Y")

Pollard1 <- Pollard1 %>% mutate(effort = Duration.hr, count = Monarchs) %>%
  select(x, y, date, count, effort)

Pollard2 <- Pollard2 %>% mutate(effort = Duration, count = Total.monarch) %>%
  select(x, y, date, count, effort)

Pollard <- full_join(Pollard1, Pollard2, by = c("x", "y", "date", "count", "effort"))

#-Seasonal Dynamics-#

#Set temporal periods
Pollard <- Pollard %>% 
  mutate(yr = format(date, "%Y"),
         mo_day = format(date, "%m-%d")) %>%
  mutate(period = ifelse(mo_day >= "03-08" & mo_day < "04-04", 1, 
                         ifelse(mo_day >= "04-19" & mo_day < "05-02", 2,
                                ifelse(mo_day >= "05-03" & mo_day < "06-06", 3, 0)))) 

#Remove observation outside of periods
Pollard <- Pollard %>% filter(period > 0 & yr >= 2016 & yr < 2019)

#-Standardize-#

Monarchs <- cbind(Pollard[,c(1:2,6:8)], rep("Pollard", length(Pollard[,1])), Pollard[,c(4:5)])
colnames(Monarchs) <- c("x", "y", "yr", "mo_day", "period", "program", "count", "effort")

save(Monarchs, file = "Monarch.Rdata")
