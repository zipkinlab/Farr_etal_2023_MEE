#----------------------------#
#----Covariate extraction----#
#----------------------------#

#TODO: 
#fix raster stack saving issue
#change data loading code
#fix order issue of node data


#-Libraries-#
library(tidyverse)
library(sf)
library(INLA)
library(TMB)

#-Load spatial domain-#

#Spring
domain_sp <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Spring.shp")
#Summer
domain_sm <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Summer.shp")

#-Load Monarch Data-#

files <- list.files(path = "~/Monarchs/DataFormatting/MonarchData/", pattern = "Monarch.Rdata", 
                    recursive = TRUE, full.names = TRUE)

Data <- rep(NA, 8)

for(i in 1:length(files)){
  load(files[i])
  Data <- rbind(Data, Monarchs)
}

Data <- Data[-1,]

#iNaturalist
load("~/Monarchs/DataFormatting/MonarchData/iNat/Monarch.Rdata")
Data <- Monarchs

#Journey North
#load("~/Monarchs/DataFormatting/MonarchData/JN/Monarch.Rdata")
#Data <- rbind(Data, Monarchs)

#eButterfly
#NULL

#NABA: Butterflies I've Seen
#NULL

#Butterflies and Moths of North America
#NULL

#NBA: NFJ
load("~/Monarchs/DataFormatting/MonarchData/NABA/NFJ/Monarch.Rdata")
Data <- merge(Data, Monarchs, all = TRUE)

#NABA: Sightings

#Pollard
load("~/Monarchs/DataFormatting/MonarchData/Pollard/Monarch.Rdata")
Data <- merge(Data, Monarchs, all = TRUE)

#Other structured data??


#-SPDF-#

Data <- st_as_sf(x = Data, coords = c("x", "y"), crs = st_crs(domain_sm))

#-Projection-#

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

domain_sp <- st_transform(domain_sp, crs = st_crs(prj))
domain_sm <- st_transform(domain_sm, crs = st_crs(prj))

Data <- st_transform(Data, crs = st_crs(prj))

#-Spatio-temporal domain-#

Data <- rbind(Data %>% filter(period < 3) %>% .[domain_sp,], 
              Data %>% filter(period == 3) %>% .[domain_sm,])

#-Visualize-#

ggplot() + 
  geom_sf(data = domain_sp, alpha = 0.5) +
  geom_sf(data = domain_sm, alpha = 0.5) +
  geom_sf(data = Data, aes(col = program)) +
  facet_wrap(~yr + period) +
  theme_bw()

#-Mesh-#

#Spring
domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))
#Summer
domain_sm_seg <- inla.sp2segment(as(domain_sm, "Spatial"))

mesh <- list()

mesh[[2]] <- mesh[[1]] <- inla.mesh.2d(boundary=domain_sp_seg,
                                       max.edge=c(36, 108),
                                       cutoff=10)

mesh[[3]] <- inla.mesh.2d(boundary=domain_sm_seg,
                          max.edge=c(36, 108),
                          cutoff=25)

sp.nodes <- mesh[[1]]$n
su.nodes <- mesh[[3]]$n

loc <- st_coordinates(Data)
#A <- inla.spde.make.A(mesh, loc = loc) #FIX ERROR

#-Covariates-#

#NDVI
NDVInames <- c(list.files(path = "./DataFormatting/NDVI/2016/Spring", pattern = "tif$", full.names = TRUE),
               list.files(path = "./DataFormatting/NDVI/2017/Spring", pattern = "tif$", full.names = TRUE),
               list.files(path = "./DataFormatting/NDVI/2018/Spring", pattern = "tif$", full.names = TRUE))
NDVIspring <- raster::stack(NDVInames)
NDVIspring <- raster::projectRaster(NDVIspring, crs = prj)
raster::writeRaster(NDVIspring, file = "NDVIspring.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
raster::removeTmpFiles(h=0)

NDVInames <- c(list.files(path = "./DataFormatting/NDVI/2016/Summer", pattern = "tif$", full.names = TRUE),
               list.files(path = "./DataFormatting/NDVI/2017/Summer", pattern = "tif$", full.names = TRUE),
               list.files(path = "./DataFormatting/NDVI/2018/Summer", pattern = "tif$", full.names = TRUE))
NDVIsummer <- raster::stack(NDVInames)
NDVIsummer <- raster::projectRaster(NDVIsummer, crs = prj)
raster::writeRaster(NDVIsummer, file = "NDVIsummer.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
raster::removeTmpFiles(h=0)

NDVIdate <- c("03-21", "04-06", "04-22", "03-22", "04-07", "04-23", "05-08", "05-24", "05-09", "05-25")

Data <- Data %>% mutate(NDVI = ifelse(mo_day < NDVIdate[1] & period < 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIspring[[1]], .),
                                      ifelse(mo_day < NDVIdate[2] & period < 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIspring[[2]], .),
                                      ifelse(mo_day < NDVIdate[3] & period < 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIspring[[3]], .),
                                      ifelse(mo_day >= NDVIdate[3] & period < 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIspring[[4]], .),
                                      ifelse(mo_day < NDVIdate[4] & period < 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIspring[[5]], .),
                                      ifelse(mo_day < NDVIdate[5] & period < 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIspring[[6]], .),       
                                      ifelse(mo_day < NDVIdate[6] & period < 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIspring[[7]], .),
                                      ifelse(mo_day >= NDVIdate[6] & period < 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIspring[[8]], .),       
                                      ifelse(mo_day < NDVIdate[4] & period < 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIspring[[9]], .),              
                                      ifelse(mo_day < NDVIdate[5] & period < 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIspring[[10]], .),
                                      ifelse(mo_day < NDVIdate[6] & period < 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIspring[[11]], .),       
                                      ifelse(mo_day >= NDVIdate[6] & period < 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIspring[[12]], .),
                                      ifelse(mo_day < NDVIdate[7] & period == 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIsummer[[1]], .),
                                      ifelse(mo_day < NDVIdate[8] & period == 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIsummer[[2]], .),
                                      ifelse(mo_day >= NDVIdate[8] & period == 3 & yr == 2016,
                                      0.0001*raster::extract(NDVIsummer[[3]], .),
                                      ifelse(mo_day < NDVIdate[9] & period == 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIsummer[[4]], .),
                                      ifelse(mo_day < NDVIdate[10] & period == 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIsummer[[5]], .),
                                      ifelse(mo_day >= NDVIdate[10] & period == 3 & yr == 2017,
                                      0.0001*raster::extract(NDVIsummer[[6]], .),  
                                      ifelse(mo_day < NDVIdate[9] & period == 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIsummer[[7]], .),
                                      ifelse(mo_day < NDVIdate[10] & period == 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIsummer[[8]], .),
                                      ifelse(mo_day >= NDVIdate[10] & period == 3 & yr == 2018,
                                      0.0001*raster::extract(NDVIsummer[[9]], .),                                        
                                      NA))))))))))))))))))))))             
                                      
                                      
#Population Density
PopDnames <- c("~/Monarchs/DataFormatting/PopulationDensity/2015/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif",
               "~/Monarchs/DataFormatting/PopulationDensity/2020/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
PopDstack <- raster::stack(PopDnames)

clipper <- st_sfc(st_polygon(list(rbind(c(-110,15),c(-110,50),c(-65,50),c(-65,15),c(-110,15)))))
clipper <- st_sf(geometry = clipper, crs = st_crs(PopDstack))

PopDstack <- raster::crop(PopDstack, clipper)
PopDstack <- raster::projectRaster(PopDstack, crs = prj)

raster::writeRaster(PopDstack, file = "PopDstack.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
#save(PopDstack, file = "PopDstack.Rdata")
#load(file = "./DataFormatting/NDVIstack.Rdata")

#Sampling Bias
#load(file = "./DataFormatting/Biasobj3.Rdata")
load("~/Monarchs/DataFormatting/MonarchData/iNat/iNatBias.Rdata")
#dyn.load(dynlib("./DataAnalysis/Bias"))
#B.out <- sdreport(obj3)

year <- c("2016", "2017", "2018")
nBias <- NULL
Predictor <- function(x,t){
  A <- as(inla.spde.make.A(mesh[[p]], st_coordinates(x)), "dgTMatrix")
  Pred <- exp(beta[t] + A %*% omg + A %*% eps[,t])
  return(Pred)
}

Data$Bias <- NA

for(p in 1:3){
  beta <- Bias[[p]]$par.fixed[1:3]
  omg <- Bias[[p]]$par.random[1:mesh[[p]]$n]
  eps <- Bias[[p]]$par.random[(mesh[[p]]$n+1):(mesh[[p]]$n + mesh[[p]]$n*3)]
  tauO <- exp(Bias[[p]]$par.fixed[5])
  tauE <- exp(Bias[[p]]$par.fixed[6])
  dim(eps) <- c(mesh[[p]]$n, 3)
  omg <- omg/tauO
  eps <- eps/tauE
  for(t in 1:3){
    #data <- Data %>% filter(yr == year[t] & period == p & program == "iNat")
    #A <- as(inla.spde.make.A(mesh[[p]], st_coordinates(data)), "dgTMatrix")
    #Pred <- exp(beta[t] + A %*% omg + A %*% eps[,t])
    Data <- Data %>% mutate(Bias = ifelse(yr == year[t] & period == p & program == "iNat", Predictor(.,t), Data$Bias))
    #Data[Data$yr == year[t] & Data$period == p & Data$program == "iNat", 9] <- Pred
    nBias <- c(nBias, exp(beta[t] + omg + eps[,t]))
  }
}

#sp.nBias <- cbind(nodes[1:(sp.nodes*3)], nodes[(sp.nodes*3+1):(sp.nodes*3+sp.nodes*3)])
#su.nBias <- nodes[(sp.nodes*3+sp.nodes*3+1):(sp.nodes*3+sp.nodes*3+su.nodes*3)]

#-Extract covariates @ nodes-#

sp.nodexy <- mesh[[1]]$loc[,1:2]
su.nodexy <- mesh[[3]]$loc[,1:2]

#Spring NDVI
sp1.nNDVI16 <- sp1.nNDVI17 <- sp1.nNDVI18 <- sp2.nNDVI16 <- sp2.nNDVI17 <- sp2.nNDVI18 <- NULL
for(k in 1:sp.nodes){
  tmp <- NULL
  for(j in 1:12){
    tmp[j] <- 0.0001*raster::extract(NDVIspring[[j]], st_coordinates(st_point(sp.nodexy[k,])))
    # if(tmp[j] < -2000){
    #   tmp[j] <- -1
    # }else{tmp[j] <- tmp[j]*0.0001}
  }
  sp1.nNDVI16[k] <- mean(tmp[1:2])
  sp1.nNDVI17[k] <- mean(tmp[5:6])
  sp1.nNDVI18[k] <- mean(tmp[9:10])
  sp2.nNDVI16[k] <- mean(tmp[3:4])
  sp2.nNDVI17[k] <- mean(tmp[7:8])
  sp2.nNDVI18[k] <- mean(tmp[11:12])
}

sp.nNDVI <- cbind(c(sp1.nNDVI16, sp1.nNDVI17, sp1.nNDVI18), 
                  c(sp2.nNDVI16, sp2.nNDVI17, sp2.nNDVI18))

#Summer NDVI
su.nNDVI16 <- su.nNDVI17 <- su.nNDVI18 <- NULL
for(k in 1:su.nodes){
  tmp <- NULL
  for(j in 1:9){
    tmp[j] <- 0.0001*raster::extract(NDVIsummer[[j]], st_coordinates(st_point(su.nodexy[k,])))
    # if(tmp[j] < -2000){
    #   tmp[j] <- -1
    # }else{tmp[j] <- tmp[j]*0.0001}
  }
  su.nNDVI16[k] <- mean(tmp[1:3])
  su.nNDVI17[k] <- mean(tmp[4:6])
  su.nNDVI18[k] <- mean(tmp[7:9])
}

su.nNDVI <- c(su.nNDVI16, su.nNDVI17, su.nNDVI18)

nNDVI <- c(sp1.nNDVI16, sp1.nNDVI17, sp1.nNDVI18, 
           sp2.nNDVI16, sp2.nNDVI17, sp2.nNDVI18,
           su.nNDVI16, su.nNDVI17, su.nNDVI18)

#Bias
# nBias <- NULL
# t_k <- rep(1:t_n, each = nodes)
# k_k <- rep(1:nodes, t_n)
# for(k in 1:(3*nodes)){
#   nBias[k] <- as.numeric(exp(beta[t_k[k]] + omg[k_k[k]] + eps[k_k[k], t_k[k]]))
# }


NodesDF <- st_as_sf(data.frame(Bias = nBias, NDVI = nNDVI, 
                               yr = as.factor(rep(rep(2016:2018, 3), c(rep(sp.nodes, 6), rep(su.nodes, 3)))),
                               period = as.factor(rep(1:3, c(sp.nodes * 3, sp.nodes* 3, su.nodes * 3))),
                               X = c(rep(sp.nodexy[,1], 2), su.nodexy[,1]),
                               Y = c(rep(sp.nodexy[,2], 2), su.nodexy[,2])), 
                    coords = c("X", "Y"), crs = prj)

#No sampling outside domain
outside <- sapply(st_intersects(NodesDF, st_union(x = domain_sm, y = domain_sp)), function(x){length(x)==0})
NodesDF$Bias[outside] <- 0

#-Arrange Data-#

Data <- Data %>% mutate(type = ifelse(program == "NFJ"|program == "Pollard", "count", "obs")) %>%
  arrange(desc(type), period, yr)

#-Save data-#

save(Data, file = "./DataFormatting/FormattedData.Rdata")
save(NodesDF, file = "./DataFormatting/FormattedNodeData.Rdata")
