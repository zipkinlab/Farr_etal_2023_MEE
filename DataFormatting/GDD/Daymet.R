#-------------------#
#-Daymet downloader-#
#-------------------#

#-----------#
#-Libraries-#
#-----------#

library(httr)

#-----------------#
#-Lat/long values-#
#-----------------#

#prj.daymet <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.31406705 +units=m +no_defs" 
prj.daymet <- "+init=epsg:4326"

#Data coords
d.coords <- st_coordinates(st_transform(Data, crs = st_crs(prj.daymet)))
d.lat <- as.numeric(d.coords[,2])
d.long <- as.numeric(d.coords[,1])
  
#Spring mesh coords
sp.nodexy <-  st_coordinates(st_transform(st_as_sf(data.frame(NA, X = mesh[[1]]$loc[,1], Y = mesh[[1]]$loc[,2]), 
         coords = c("X", "Y"), crs = prj), crs = st_crs(prj.daymet)))
sp.lat <- sp.nodexy[,2]
sp.long <- sp.nodexy[,1]
  
#Summer mesh coords
sm.nodexy <-  st_coordinates(st_transform(st_as_sf(data.frame(NA, X = mesh[[3]]$loc[,1], Y = mesh[[3]]$loc[,2]), 
                                                   coords = c("X", "Y"), crs = prj), crs = st_crs(prj.daymet)))
sm.lat <- sm.nodexy[,2]
sm.long <- sm.nodexy[,1]

#-------#
#-Dates-#
#-------#

#Data dates
d.start_date <- d.end_date <- paste(Data$yr, "-", Data$mo_day, sep = "")
  
#Spring mesh dates
sp.start_date <- c("2016-03-08", "2016-04-19", "2017-03-08", "2017-04-19", "2018-03-08", "2018-04-19")
sp.end_date <- c("2016-04-04", "2016-05-02", "2017-04-04", "2017-05-02", "2018-04-04", "2018-05-02")
  
#Summer mesh dates

#----------#
#-Function-#
#----------#

source("~/Monarchs/DataFormatting/GDD/read_Daymet.R")

daymet <- function(long,lat,start_date,end_date){
  query <- list("lat" = lat,
                "lon" = long,
                "vars" = "tmax,tmin,prcp",
                "start" = start_date,
                "end" = end_date)
  GET(url = "https://daymet.ornl.gov/single-pixel/api/data?lat={}&lon={}", 
      query = query, write_disk(path = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv", overwrite = TRUE))
  return(unlist(read_daymet(file = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv")[1:3,9]))
}

#---------------#
#-Download data-#
#---------------#

#Data
daymet.data <- mapply(FUN = daymet, long = d.long, lat = d.lat, start_date = d.start_date, end_date = d.end_date)
save(daymet.data, file = "~/Monarchs/DataFormatting/GDD/Data.Rdata")

#Spring mesh
daymet.sp <- mapply(FUN = daymet, long = sp.long[1:2], lat = sp.lat[1:2], start_date = sp.start_date[1], end_date = sp.end_date[1])

GET(url = "https://daymet.ornl.gov/single-pixel/api/data", 
    query = query, write_disk(path = "", overwrite = TRUE))

#Summer mesh
GET(url = "https://daymet.ornl.gov/single-pixel/api/data", 
    query = query, write_disk(path = "", overwrite = TRUE))


