#----------------------------#
#----Covariate extraction----#
#----------------------------#

#TODO: 
#fix raster stack saving issue
#change data loading code
#fix order issue of node data
#fix gdd

#-----------#
#-Libraries-#
#-----------#

library(tidyverse)
library(sf)
library(INLA)
library(TMB)
library(httr)

#---------------------#
#-Load spatial domain-#
#---------------------#

#Spring
domain_sp <- st_read("./SpatialDomains/Domain_Spring.shp")
#Summer
domain_su <- st_read("./SpatialDomains/Domain_Summer.shp")

#-------------------#
#-Load Monarch Data-#
#-------------------#

files <- list.files(path = "./MonarchData/", pattern = "Monarch.Rdata", 
                    recursive = TRUE, full.names = TRUE)

Data <- rep(NA, 8)

for(i in 1:length(files)){
  load(files[i])
  Data <- rbind(Data, Monarchs)
}

Data <- Data[-1,]

#----------------#
#-Load Bias Data-#
#----------------#

load("./MonarchData/iNat/iNatData.Rdata")

#-SPDF-#

Data <- st_as_sf(x = Data, coords = c("x", "y"), crs = st_crs(domain_su))

iNat <- st_as_sf(x = iNat, coords = c("x", "y"), crs = st_crs(domain_su))

#-Projection-#

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

domain_sp <- st_transform(domain_sp, crs = st_crs(prj))
domain_su <- st_transform(domain_su, crs = st_crs(prj))

Data <- st_transform(Data, crs = st_crs(prj))

iNat <- st_transform(iNat, crs = st_crs(prj))

#-Spatio-temporal domain-#

Data <- rbind(Data %>% filter(period < 3) %>% .[domain_sp,], 
              Data %>% filter(period == 3) %>% .[domain_su,])

iNat <- rbind(iNat[domain_sp,], iNat[domain_su,])

#------#
#-Mesh-#
#------#

#Spring
domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))
#Summer
domain_su_seg <- inla.sp2segment(as(domain_su, "Spatial"))

mesh <- list()

mesh[[1]] <- inla.mesh.2d(boundary=domain_sp_seg,
                          max.edge = 75,
                          cutoff=30)

mesh[[2]] <- inla.mesh.2d(boundary=domain_sp_seg,
                          max.edge = 75,
                          cutoff=30)

mesh[[3]] <- inla.mesh.2d(boundary=domain_su_seg,
                          max.edge = 75,
                          cutoff = 30)


nodes1 <- mesh[[1]]$n
nodes2 <- mesh[[2]]$n
nodes3 <- mesh[[3]]$n

loc1 <- st_coordinates(Data %>% filter(period == 1))
loc2 <- st_coordinates(Data %>% filter(period == 2))
loc3 <- st_coordinates(Data %>% filter(period == 3))

nodexy1 <- mesh[[1]]$loc[,1:2]
nodexy2 <- mesh[[2]]$loc[,1:2]
nodexy3 <- mesh[[3]]$loc[,1:2]

nodexy <- rbind(nodexy1, nodexy2, nodexy3)

NodeDF <- data.frame(Year  = as.factor(rep(2016:2018, sum(nodes1, nodes2, nodes3))), 
                     period = as.factor(rep(1:3, c(nodes1 * 3, nodes2 * 3, nodes3 * 3))),
                     longitude = (nodexy %x% rep(1, 3))[,1], 
                     latitude = (nodexy %x% rep(1, 3))[,2]) %>% 
  st_as_sf(., coords = c("longitude", "latitude"), crs = st_crs(prj))

A1 <- inla.spde.make.A(mesh[[1]], loc = loc1)
A2 <- inla.spde.make.A(mesh[[2]], loc = loc2)
A3 <- inla.spde.make.A(mesh[[3]], loc = loc3)

#-Covariates-#

#NDVI

#Code only needs run once
# NDVInames <- c(list.files(path = "./DataFormatting/NDVI/2016/Spring", pattern = "tif$", full.names = TRUE),
#                list.files(path = "./DataFormatting/NDVI/2017/Spring", pattern = "tif$", full.names = TRUE),
#                list.files(path = "./DataFormatting/NDVI/2018/Spring", pattern = "tif$", full.names = TRUE))
# NDVIspring <- raster::stack(NDVInames)
# NDVIspring <- raster::projectRaster(NDVIspring, crs = prj)
# raster::writeRaster(NDVIspring, file = "NDVIspring.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
# raster::removeTmpFiles(h=0)
# 
# NDVInames <- c(list.files(path = "./DataFormatting/NDVI/2016/Summer", pattern = "tif$", full.names = TRUE),
#                list.files(path = "./DataFormatting/NDVI/2017/Summer", pattern = "tif$", full.names = TRUE),
#                list.files(path = "./DataFormatting/NDVI/2018/Summer", pattern = "tif$", full.names = TRUE))
# NDVIsummer <- raster::stack(NDVInames)
# NDVIsummer <- raster::projectRaster(NDVIsummer, crs = prj)
# raster::writeRaster(NDVIsummer, file = "NDVIsummer.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
# raster::removeTmpFiles(h=0)

NDVIspring <- raster::stack("NDVIspring.tif")
NDVIsummer <- raster::stack("NDVIsummer.tif")

NDVIdate <- c("03-21", "04-06", "04-22", "03-22", "04-07", "04-23", "05-08", "05-24", "05-09", "05-25")

Data$NDVI <- NA

for(i in 1:dim(Data)[1]){
Data$NDVI[i] <- ifelse(Data$mo_day[i] < NDVIdate[1] & Data$period[i] < 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIspring[[1]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[2] & Data$period[i] < 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIspring[[2]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[3] & Data$period[i] < 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIspring[[3]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] >= NDVIdate[3] & Data$period[i] < 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIspring[[4]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[4] & Data$period[i] < 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIspring[[5]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[5] & Data$period[i] < 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIspring[[6]], Data[i,], buffer = 10)), na.rm = TRUE),       
ifelse(Data$mo_day[i] < NDVIdate[6] & Data$period[i] < 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIspring[[7]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] >= NDVIdate[6] & Data$period[i] < 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIspring[[8]], Data[i,], buffer = 10)), na.rm = TRUE),       
ifelse(Data$mo_day[i] < NDVIdate[4] & Data$period[i] < 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIspring[[9]], Data[i,], buffer = 10)), na.rm = TRUE),              
ifelse(Data$mo_day[i] < NDVIdate[5] & Data$period[i] < 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIspring[[10]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[6] & Data$period[i] < 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIspring[[11]], Data[i,], buffer = 10)), na.rm = TRUE),       
ifelse(Data$mo_day[i] >= NDVIdate[6] & Data$period[i] < 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIspring[[12]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[7] & Data$period[i] == 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIsummer[[1]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[8] & Data$period[i] == 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIsummer[[2]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] >= NDVIdate[8] & Data$period[i] == 3 & Data$yr[i] == 2016,
0.0001*mean(unlist(raster::extract(NDVIsummer[[3]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[9] & Data$period[i] == 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIsummer[[4]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[10] & Data$period[i] == 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIsummer[[5]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] >= NDVIdate[10] & Data$period[i] == 3 & Data$yr[i] == 2017,
0.0001*mean(unlist(raster::extract(NDVIsummer[[6]], Data[i,], buffer = 10)), na.rm = TRUE),  
ifelse(Data$mo_day[i] < NDVIdate[9] & Data$period[i] == 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIsummer[[7]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] < NDVIdate[10] & Data$period[i] == 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIsummer[[8]], Data[i,], buffer = 10)), na.rm = TRUE),
ifelse(Data$mo_day[i] >= NDVIdate[10] & Data$period[i] == 3 & Data$yr[i] == 2018,
0.0001*mean(unlist(raster::extract(NDVIsummer[[9]], Data[i,], buffer = 10)), na.rm = TRUE),                                        
NA)))))))))))))))))))))             
}

# Data <- Data %>% mutate(NDVI = ifelse(mo_day < NDVIdate[1] & period < 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIspring[[1]], .),
#                                       ifelse(mo_day < NDVIdate[2] & period < 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIspring[[2]], .),
#                                       ifelse(mo_day < NDVIdate[3] & period < 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIspring[[3]], .),
#                                       ifelse(mo_day >= NDVIdate[3] & period < 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIspring[[4]], .),
#                                       ifelse(mo_day < NDVIdate[4] & period < 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIspring[[5]], .),
#                                       ifelse(mo_day < NDVIdate[5] & period < 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIspring[[6]], .),       
#                                       ifelse(mo_day < NDVIdate[6] & period < 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIspring[[7]], .),
#                                       ifelse(mo_day >= NDVIdate[6] & period < 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIspring[[8]], .),       
#                                       ifelse(mo_day < NDVIdate[4] & period < 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIspring[[9]], .),              
#                                       ifelse(mo_day < NDVIdate[5] & period < 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIspring[[10]], .),
#                                       ifelse(mo_day < NDVIdate[6] & period < 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIspring[[11]], .),       
#                                       ifelse(mo_day >= NDVIdate[6] & period < 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIspring[[12]], .),
#                                       ifelse(mo_day < NDVIdate[7] & period == 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIsummer[[1]], .),
#                                       ifelse(mo_day < NDVIdate[8] & period == 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIsummer[[2]], .),
#                                       ifelse(mo_day >= NDVIdate[8] & period == 3 & yr == 2016,
#                                       0.0001*mean(raster::extract(NDVIsummer[[3]], .),
#                                       ifelse(mo_day < NDVIdate[9] & period == 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIsummer[[4]], .),
#                                       ifelse(mo_day < NDVIdate[10] & period == 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIsummer[[5]], .),
#                                       ifelse(mo_day >= NDVIdate[10] & period == 3 & yr == 2017,
#                                       0.0001*mean(raster::extract(NDVIsummer[[6]], .),  
#                                       ifelse(mo_day < NDVIdate[9] & period == 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIsummer[[7]], .),
#                                       ifelse(mo_day < NDVIdate[10] & period == 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIsummer[[8]], .),
#                                       ifelse(mo_day >= NDVIdate[10] & period == 3 & yr == 2018,
#                                       0.0001*mean(raster::extract(NDVIsummer[[9]], .),                                        
#                                       NA))))))))))))))))))))))             
                                      
 
#GDD

#prj.daymet <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +a=6378137 +b=6356752.31406705 +units=m +no_defs" 
prj.daymet <- "+init=epsg:4326"

#Data coords
d.coords <- st_coordinates(st_transform(Data, crs = st_crs(prj.daymet)))
d.lat <- as.numeric(d.coords[,2])
d.long <- as.numeric(d.coords[,1])

#Data dates
d.start_date <- paste(as.Date(paste0(Data$yr, "-", Data$mo_day))-14)

d.end_date <- paste0(Data$yr, "-", Data$mo_day)

source("./GDD/gdd_fun.R")

daymet <- function(long,lat,start_date,end_date){
  query <- list("lat" = lat,
                "lon" = long,
                # "vars" = "tmax,tmin,prcp",
                "vars" = "tmax,tmin",
                "start" = start_date,
                "end" = end_date)
  GET(url = "https://daymet.ornl.gov/single-pixel/api/data?lat={}&lon={}", 
      query = query, write_disk(path = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv", overwrite = TRUE))
  return(read_gdd(file = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv", date = end_date, lat = lat, long = long))
}

daymet.data <- do.call(rbind, mapply(FUN = daymet, long = d.long, lat = d.lat, 
                                     start_date = d.start_date, end_date = d.end_date, SIMPLIFY = FALSE))

Data$gdd1 <- daymet.data$gdd1
Data$gdd2 <- daymet.data$gdd2

# Data <- Data %>% mutate(lat = as.numeric(st_coordinates(st_transform(Data, crs = st_crs(prj.daymet)))[,2]), lon = as.numeric(st_coordinates(st_transform(Data, crs = st_crs(prj.daymet)))[,1]))
# 
# Data <- Data %>% mutate(date = paste0(yr, "-", mo_day)) %>% full_join(., daymet.data, by = c("lat", "lon", "date"))

#Population Density
PopDnames <- c("~/Monarchs/DataFormatting/PopulationDensity/2015/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif",
               "~/Monarchs/DataFormatting/PopulationDensity/2020/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
PopDstack <- raster::stack(PopDnames)

clipper <- st_sfc(st_polygon(list(rbind(c(-110,15),c(-110,50),c(-65,50),c(-65,15),c(-110,15)))))
clipper <- st_sf(geometry = clipper, crs = st_crs(PopDstack))

PopDstack <- raster::crop(PopDstack, clipper)
PopDstack <- raster::projectRaster(PopDstack, crs = prj)

raster::writeRaster(PopDstack, file = "PopDstack.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
raster::removeTmpFiles(h=0)
PopDstack <- raster::stack("PopDstack.tif")

Data$PopD <- NA

for(i in 1:dim(Data)[1]){
  Data$PopD[i] <- mean(unlist(raster::extract(PopDstack, Data[i,], buffer = 10)), na.rm = TRUE)
}

#Sampling Bias
#load(file = "./DataFormatting/Biasobj3.Rdata")
#load("~/Monarchs/DataFormatting/MonarchData/iNat/iNatBias.Rdata")
#load("~/Monarchs/DataFormatting/MonarchData/iNat/iNatBiasOutput.Rdata")

Data <- Data %>% mutate(Bias = ifelse(yr == "2016", lengths(st_intersects(x = st_buffer(., dist = 10), y = iNat %>% filter(yr == "2016"))),
                               ifelse(yr == "2017", lengths(st_intersects(x = st_buffer(., dist = 10), y = iNat %>% filter(yr == "2017"))),
                                                    lengths(st_intersects(x = st_buffer(., dist = 10), y = iNat %>% filter(yr == "2018")))))) 



#-Extract covariates @ nodes-#

nodes <- dim(NodeDF)[1]

for(i in 1:dim(NodeDF)[1]){
  if(NodeDF$Year[i] == "2016"){
    NodeDF$Bias[i] <- lengths(st_intersects(x = st_buffer(NodeDF[i,], dist = 10), y = iNat %>% filter(yr == "2016")))
  }
  if(NodeDF$Year[i] == "2017"){
    NodeDF$Bias[i] <- lengths(st_intersects(x = st_buffer(NodeDF[i,], dist = 10), y = iNat %>% filter(yr == "2017")))
  }
  if(NodeDF$Year[i] == "2018"){
    NodeDF$Bias[i] <- lengths(st_intersects(x = st_buffer(NodeDF[i,], dist = 10), y = iNat %>% filter(yr == "2018")))
  }
  if(i %in% seq(1,nodes,1000)){
    print(paste0("iteration", i))
  }
}



#dyn.load(dynlib("./DataAnalysis/Bias"))
#B.out <- sdreport(obj3)

# year <- c("2016", "2017", "2018")
# nBias <- NULL
# # Predictor <- function(x,t){
# #   A <- as(inla.spde.make.A(mesh[[p]], st_coordinates(x)), "dgTMatrix")
# #   Pred <- exp(beta0 + A %*% omg)
# #   return(Pred)
# # }
# 
# # Data$Bias <- NA
# # j <- 1
# 
# for(p in 1:3){
#   for(t in 1:3){
#     # beta0 <- output[[j]]$par.fixed
#     # omg <- output[[j+1]]$par.random
#     # tau <- exp(output[[j+1]]$par.fixed[2])
#     # omg <- omg/tau
#     # data <- Data %>% filter(yr == year[t] & period == p & (program == "iNat"|program == "JN"))
#     # A <- as(inla.spde.make.A(mesh[[p]], st_coordinates(data)), "dgTMatrix")
#     # Pred <- exp(beta0 + A %*% omg)
#     # Data <- Data %>% mutate(Bias = ifelse(yr == year[t] & period == p & (program == "iNat"|program == "JN"), Predictor(.,t), Data$Bias))
#     # #Data[Data$yr == year[t] & Data$period == p & Data$program == "iNat", 9] <- Pred
#     # nBias <- c(nBias, exp(beta0 + omg))
#     # j <- j+2
#     nBias <- c(nBias, lengths(st_intersects(x = st_buffer(st_sfc(st_multipoint(meshxy[[p]])), dist = sqrt(10/pi)), y = iNat %>% filter(yr == year[t]))))
#   }
# }

#sp.nBias <- cbind(nodes[1:(sp.nodes*3)], nodes[(sp.nodes*3+1):(sp.nodes*3+sp.nodes*3)])
#su.nBias <- nodes[(sp.nodes*3+sp.nodes*3+1):(sp.nodes*3+sp.nodes*3+su.nodes*3)]



#Spring NDVI
sp1.nNDVI16 <- sp1.nNDVI17 <- sp1.nNDVI18 <- sp2.nNDVI16 <- sp2.nNDVI17 <- sp2.nNDVI18 <- sp1.popd  <- sp2.popd <- NULL
for(k in 1:nodes1){
  tmp <- NULL
  for(j in c(1,2,5,6,9,10)){
    tmp[j] <- 0.0001*mean(unlist(raster::extract(NDVIspring[[j]], st_coordinates(st_point(nodexy1[k,])), buffer = 10)), na.rm = TRUE)
    # if(tmp[j] < -2000){
    #   tmp[j] <- -1
    # }else{tmp[j] <- tmp[j]*0.0001}
  }
  sp1.nNDVI16[k] <- mean(tmp[1:2])
  sp1.nNDVI17[k] <- mean(tmp[5:6])
  sp1.nNDVI18[k] <- mean(tmp[9:10])
  sp1.popd[k] <- mean(unlist(raster::extract(PopDstack, st_coordinates(st_point(nodexy1[k,])), buffer = 10)), na.rm = TRUE)
}

for(k in 1:nodes2){
  tmp <- NULL
  for(j in c(3,4,7,8,11,12)){
    tmp[j] <- 0.0001*mean(unlist(raster::extract(NDVIspring[[j]], st_coordinates(st_point(nodexy2[k,])), buffer = 10)), na.rm = TRUE)
    # if(tmp[j] < -2000){
    #   tmp[j] <- -1
    # }else{tmp[j] <- tmp[j]*0.0001}
  }
  sp2.nNDVI16[k] <- mean(tmp[3:4])
  sp2.nNDVI17[k] <- mean(tmp[7:8])
  sp2.nNDVI18[k] <- mean(tmp[11:12])
  sp2.popd[k] <- mean(unlist(raster::extract(PopDstack, st_coordinates(st_point(nodexy2[k,])), buffer = 10)), na.rm = TRUE)
}

#Summer NDVI
su.nNDVI16 <- su.nNDVI17 <- su.nNDVI18 <- su.popd <- NULL
for(k in 1:nodes3){
  tmp <- NULL
  for(j in 1:9){
    tmp[j] <- 0.0001*mean(unlist(raster::extract(NDVIsummer[[j]], st_coordinates(st_point(nodexy3[k,])), buffer = 10)), na.rm = TRUE)
    # if(tmp[j] < -2000){
    #   tmp[j] <- -1
    # }else{tmp[j] <- tmp[j]*0.0001}
  }
  su.nNDVI16[k] <- mean(tmp[1:3])
  su.nNDVI17[k] <- mean(tmp[4:6])
  su.nNDVI18[k] <- mean(tmp[7:9])
  su.popd[k] <- mean(unlist(raster::extract(PopDstack, st_coordinates(st_point(nodexy3[k,])), buffer = 10)), na.rm = TRUE)
}

nNDVI <- c(sp1.nNDVI16, sp1.nNDVI17, sp1.nNDVI18, 
           sp2.nNDVI16, sp2.nNDVI17, sp2.nNDVI18,
           su.nNDVI16, su.nNDVI17, su.nNDVI18)

nPopD <- c(rep(sp1.popd, 3), rep(sp2.popd, 3), rep(su.popd, 3))

NodeDF$NDVI <- nNDVI

NodeDF$PopD <- nPopD

#Bias
# nBias <- NULL
# t_k <- rep(1:t_n, each = nodes)
# k_k <- rep(1:nodes, t_n)
# for(k in 1:(3*nodes)){
#   nBias[k] <- as.numeric(exp(beta[t_k[k]] + omg[k_k[k]] + eps[k_k[k], t_k[k]]))
# }


# NodesDF <- st_as_sf(data.frame(Bias = nBias, NDVI = nNDVI, PopD = nPopD,
#                                yr = as.factor(rep(rep(2016:2018, 3), c(rep(nodes1, 3), rep(nodes2, 3), rep(nodes3, 3)))),
#                                period = as.factor(rep(1:3, c(nodes1 * 3, nodes2 * 3, nodes3 * 3))),
#                                X = c(nodexy1[,1], nodexy2[,1], nodexy3[,1]),
#                                Y = c(nodexy1[,2], nodexy2[,2], nodexy3[,2])), 
#                     coords = c("X", "Y"), crs = prj)

#No sampling outside domain
outside <- sapply(st_intersects(NodeDF, st_union(x = domain_su, y = domain_sp)), function(x){length(x)==0})
NodeDF$Bias2 <- NodeDF$Bias
NodeDF$Bias2[outside] <- 0

#-Arrange Data-#

Data <- Data %>% mutate(type = ifelse(program == "NFJ"|program == "Pollard", "count", "obs")) %>%
  arrange(desc(type), period, yr)

#-Save data-#

save(Data, file = "FormattedData.Rdata")
save(NodeDF, file = "FormattedNodeData.Rdata")
