#----Create intensity field for
#----iNat butterflies for bias
#----correction.
#----

#-Libraries-#

library(tidyverse)
library(sf)
library(INLA)
library(TMB)

#Function to create dual mesh (original code from xxx)
source("~/Monarchs/DataAnalysis/rDualMesh.R")

#-Load spatial domain-#

#Spring
domain_sp <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Spring.shp")
#Summer
domain_sm <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Summer.shp")

#-Load Bias Data-#

load("~/Monarchs/DataFormatting/MonarchData/iNat/iNatData.Rdata")
Data <- iNat

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
  geom_sf(data = Data) +
  facet_wrap(~yr + period) +
  theme_bw()

#-Build mesh-#

#Spring
domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))
#Summer
domain_sm_seg <- inla.sp2segment(as(domain_sm, "Spatial"))


mesh_sp <- inla.mesh.2d(boundary=domain_sp_seg,
                        max.edge=c(36, 108),
                        cutoff=10)

mesh_sm <- inla.mesh.2d(boundary=domain_sm_seg,
                        max.edge=c(36, 108),
                        cutoff=25)

#Create SPDE objects (Matern covariance)
spde_sp <- inla.spde2.matern(mesh_sp)

spde_sm <- inla.spde2.matern(mesh_sm)

#-Compile data-#

#Indices
nodes_sp <- mesh_sp$n
nodes_sm <- mesh_sm$n
nobs_sp1 <- sum(Data$period == 1)
nobs_sp2 <- sum(Data$period == 2)
nobs_sm <- sum(Data$period == 3)
t_n <- Data %>% summarize(out = n_distinct(yr)) %>% select(out) %>% .$out
t_i_sp1 <- Data %>% filter(period == 1) %>% group_by(yr) %>% summarize(out = n()) %>% select(out) %>% .$out
t_i_sp2 <- Data %>% filter(period == 2) %>% group_by(yr) %>% summarize(out = n()) %>% select(out) %>% .$out
t_i_sm <- Data %>% filter(period == 3) %>% group_by(yr) %>% summarize(out = n()) %>% select(out) %>% .$out
t_i_sp1 <- rep(0:(t_n-1), t_i_sp1)
t_i_sp2 <- rep(0:(t_n-1), t_i_sp2)
t_i_sm <- rep(0:(t_n-1), t_i_sm)

#Create dual mesh for weights
dmesh_sp <- book.mesh.dual(mesh_sp)
dmesh_sm <- book.mesh.dual(mesh_sm)
dmesh_sp <- st_as_sf(dmesh_sp)
dmesh_sm <- st_as_sf(dmesh_sm)
st_crs(dmesh_sp) <- st_crs(domain_sp)
st_crs(dmesh_sm) <- st_crs(domain_sm)

#Weights of each node (area of dual mesh)
weight_sp <- sapply(1:dim(dmesh_sp)[1], function(i) {
  if (st_intersects(dmesh_sp[i,], domain_sp, sparse = FALSE))
    return(st_area(st_intersection(dmesh_sp[i,], domain_sp)))
  else return(0)
})

weight_sm <- sapply(1:dim(dmesh_sm)[1], function(i) {
  if (st_intersects(dmesh_sm[i,], domain_sm, sparse = FALSE))
    return(st_area(st_intersection(dmesh_sm[i,], domain_sm)))
  else return(0)
})

#Projection matrix of observations
A_sp1 <- as(inla.spde.make.A(mesh_sp, st_coordinates(Data %>% filter(period == 1))), "dgTMatrix")
A_sp2 <- as(inla.spde.make.A(mesh_sp, st_coordinates(Data %>% filter(period == 2))), "dgTMatrix")
A_sm <- as(inla.spde.make.A(mesh_sm, st_coordinates(Data %>% filter(period == 3))), "dgTMatrix")


#TMB data
data_sp1 <- list("nodes" = mesh_sp$n, "nobs" = nobs_sp1, "t_n" = t_n, "t_i" = t_i_sp1, "weight" = weight_sp,
                 "A" = A_sp1, "spde" = spde_sp$param.inla[c("M0","M1","M2")])

data_sp2 <- list("nodes" = mesh_sp$n, "nobs" = nobs_sp2, "t_n" = t_n, "t_i" = t_i_sp2, "weight" = weight_sp,
                 "A" = A_sp2, "spde" = spde_sp$param.inla[c("M0","M1","M2")])

data_sm <- list("nodes" = mesh_sm$n, "nobs" = nobs_sm, "t_n" = t_n, "t_i" = t_i_sm, "weight" = weight_sm,
                "A" = A_sm, "spde" = spde_sm$param.inla[c("M0","M1","M2")])



#-TMB-#

#Compile TMB code (only once)
#compile("./DataAnalysis/Bias.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/Bias"))

#-Early Spring-#

#Parameters to estimate
params1 <- list("beta" = rep(0, t_n), "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
                "omega" = rep(0, nodes_sp), "epsilon" = matrix(0, nrow = nodes_sp, ncol = t_n))

#PHASE 1: Fit fixed parameters
map.eps <- factor(rep(NA, nodes_sp*t_n))
dim(map.eps) <- c(nodes_sp, t_n)
map1 <- list("log_kappa" = as.factor(NA), "log_tau_O" = as.factor(NA), "log_tau_E" = as.factor(NA),
             "omega" = factor(rep(NA, nodes_sp)), "epsilon" = map.eps)

#Make AD objective function
obj1 <- MakeADFun(data = data_sp1, parameters = params1, map = map1, DLL="Bias")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#PHASE 2: Fit random effects
params2 <- list("beta" = opt1$par[1:3], "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
                "omega" = rep(0, nodes_sp), "epsilon" = matrix(0, nrow = nodes_sp, ncol = t_n))

map2 <- list("beta" = factor(rep(NA, t_n)))

#Define random effects
random = c("omega", "epsilon")

#Make AD objective function
obj2 <- MakeADFun(data = data_sp1, parameters = params2, map = map2, random = random, DLL="Bias")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
out2 <- sdreport(obj2)

#PHASE 3: Final run
par.omg <- out2$par.random[1:nodes_sp]
par.eps <- out2$par.random[(nodes_sp+1):(nodes_sp + nodes_sp*t_n)]
dim(par.eps) <- c(nodes_sp, t_n)
params3 <- list("beta" = opt1$par[1:3], "log_kappa" = opt2$par[1], 
                "log_tau_O" = opt2$par[2], "log_tau_E" = opt2$par[3],
                "omega" = par.omg, "epsilon" = par.eps)

#Make AD objective function
obj3 <- MakeADFun(data = data_sp1, parameters = params3, random = random, DLL="Bias")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
out3 <- sdreport(obj3)

#Parameter values
beta <- out3$par.fixed[1:3]
omg <- out3$par.random[1:nodes_sp]
eps <- out3$par.random[(nodes_sp+1):(nodes_sp + nodes_sp*t_n)]
tauO <- exp(out3$par.fixed[5])
tauE <- exp(out3$par.fixed[6])
dim(eps) <- c(nodes_sp, t_n)

omg <- omg/tauO
eps <- eps/tauE

#Predicted values
intensity <- raster::raster(domain_sp, resolution = c(1, 1))
Grid <- coordinates(intensity)
GridA <- inla.spde.make.A(mesh_sp, loc = Grid)

df <- data.frame(cbind(t(rep(NA, t_n + 1))))
colnames(df) <- c("value", "x", "y", "yr")

year <- c("2016", "2017", "2018")

for(t in 1:t_n){
  Pred <- exp(beta[t] + GridA %*% omg + GridA %*% eps[,t])
  raster::values(intensity) <- as.vector(Pred)
  intensity <- as(intensity, "SpatialPixelsDataFrame")
  intensity <- st_as_sf(intensity)
  intensity <- intensity[domain_sp,]
  intensity <- as.data.frame(as(intensity, "Spatial"))
  intensity$yr <- year[t]
  colnames(intensity) <- c("value", "x", "y", "yr")
  df <- rbind(df, intensity)
  rm(intensity)
  intensity <- raster::raster(domain_sp, resolution = c(1, 1))
}

df <- df[-1,]

Figure1 <- ggplot(data = df) + 
  geom_sf(data = Data %>% filter(period == 1), size = 0.25, col = "white") + 
  geom_tile(data = df, aes(x = x, y = y, fill = value), alpha = 0.75) +
  #geom_sf(data = domain, alpha = 0, linetype = "dashed", lwd = 1) +
  scale_fill_gradient2(low = "#EDEDED", mid = "#FFCCCC", high = "#660000",
                       labels = c(0,0.05,0.15,0.2), breaks = c(0,0.05,0.15,0.2),
                       limits=c(0, 0.2), oob = scales::squish) +
  facet_wrap(~yr) +
  theme_bw()

Figure1

EarlySpring <- out3

#-Late Spring-#

#Parameters to estimate
params1 <- list("beta" = rep(0, t_n), "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
                "omega" = rep(0, nodes_sp), "epsilon" = matrix(0, nrow = nodes_sp, ncol = t_n))

#PHASE 1: Fit fixed parameters
map.eps <- factor(rep(NA, nodes_sp*t_n))
dim(map.eps) <- c(nodes_sp, t_n)
map1 <- list("log_kappa" = as.factor(NA), "log_tau_O" = as.factor(NA), "log_tau_E" = as.factor(NA),
             "omega" = factor(rep(NA, nodes_sp)), "epsilon" = map.eps)

#Make AD objective function
obj1 <- MakeADFun(data = data_sp2, parameters = params1, map = map1, DLL="Bias")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#PHASE 2: Fit random effects
params2 <- list("beta" = opt1$par[1:3], "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
                "omega" = rep(0, nodes_sp), "epsilon" = matrix(0, nrow = nodes_sp, ncol = t_n))

map2 <- list("beta" = factor(rep(NA, t_n)))

#Define random effects
random = c("omega", "epsilon")

#Make AD objective function
obj2 <- MakeADFun(data = data_sp2, parameters = params2, map = map2, random = random, DLL="Bias")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
out2 <- sdreport(obj2)

#PHASE 3: Final run
par.omg <- out2$par.random[1:nodes_sp]
par.eps <- out2$par.random[(nodes_sp+1):(nodes_sp + nodes_sp*t_n)]
dim(par.eps) <- c(nodes_sp, t_n)
params3 <- list("beta" = opt1$par[1:3], "log_kappa" = opt2$par[1], 
                "log_tau_O" = opt2$par[2], "log_tau_E" = opt2$par[3],
                "omega" = par.omg, "epsilon" = par.eps)

#Make AD objective function
obj3 <- MakeADFun(data = data_sp2, parameters = params3, random = random, DLL="Bias")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
out3 <- sdreport(obj3)

#Parameter values
beta <- out3$par.fixed[1:3]
omg <- out3$par.random[1:nodes_sp]
eps <- out3$par.random[(nodes_sp+1):(nodes_sp + nodes_sp*t_n)]
tauO <- exp(out3$par.fixed[5])
tauE <- exp(out3$par.fixed[6])
dim(eps) <- c(nodes_sp, t_n)

omg <- omg/tauO
eps <- eps/tauE

#Predicted values
intensity <- raster::raster(domain_sp, resolution = c(1, 1))
Grid <- coordinates(intensity)
GridA <- inla.spde.make.A(mesh_sp, loc = Grid)

df <- data.frame(cbind(t(rep(NA, t_n + 1))))
colnames(df) <- c("value", "x", "y", "yr")

year <- c("2016", "2017", "2018")

for(t in 1:t_n){
  Pred <- exp(beta[t] + GridA %*% omg + GridA %*% eps[,t])
  raster::values(intensity) <- as.vector(Pred)
  intensity <- as(intensity, "SpatialPixelsDataFrame")
  intensity <- st_as_sf(intensity)
  intensity <- intensity[domain_sp,]
  intensity <- as.data.frame(as(intensity, "Spatial"))
  intensity$yr <- year[t]
  colnames(intensity) <- c("value", "x", "y", "yr")
  df <- rbind(df, intensity)
  rm(intensity)
  intensity <- raster::raster(domain_sp, resolution = c(1, 1))
}

df <- df[-1,]

Figure2 <- ggplot(data = df) + 
  geom_sf(data = Data %>% filter(period == 2), size = 0.25, col = "white") + 
  geom_tile(data = df, aes(x = x, y = y, fill = value), alpha = 0.75) +
  #geom_sf(data = domain, alpha = 0, linetype = "dashed", lwd = 1) +
  scale_fill_gradient2(low = "#EDEDED", mid = "#FFCCCC", high = "#660000",
                       labels = c(0,0.05,0.15,0.2), breaks = c(0,0.05,0.15,0.2),
                       limits=c(0, 0.2), oob = scales::squish) +
  facet_wrap(~yr) +
  theme_bw()

Figure2

LateSpring <- out3

#-Early Summer-#

#Parameters to estimate
params1 <- list("beta" = rep(0, t_n), "log_kappa" = -3, "log_tau_O" = 2, "log_tau_E" = 2,
                "omega" = rep(0, nodes_sm), "epsilon" = matrix(0, nrow = nodes_sm, ncol = t_n))

#PHASE 1: Fit fixed parameters
map.eps <- factor(rep(NA, nodes_sm*t_n))
dim(map.eps) <- c(nodes_sm, t_n)
map1 <- list("log_kappa" = as.factor(NA), "log_tau_O" = as.factor(NA), "log_tau_E" = as.factor(NA),
             "omega" = factor(rep(NA, nodes_sm)), "epsilon" = map.eps)

#Make AD objective function
obj1 <- MakeADFun(data = data_sm, parameters = params1, map = map1, DLL="Bias")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#PHASE 2: Fit random effects
params2 <- list("beta" = opt1$par[1:3], "log_kappa" = -3, "log_tau_O" = 2, "log_tau_E" = 2,
                "omega" = rep(0, nodes_sm), "epsilon" = matrix(0, nrow = nodes_sm, ncol = t_n))

map2 <- list("beta" = factor(rep(NA, t_n)))

#Define random effects
random = c("omega", "epsilon")

#Make AD objective function
obj2 <- MakeADFun(data = data_sm, parameters = params2, map = map2, random = random, DLL="Bias")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
out2 <- sdreport(obj2)

#PHASE 3: Final run
par.omg <- out2$par.random[1:nodes_sm]
par.eps <- out2$par.random[(nodes_sm+1):(nodes_sm + nodes_sm*t_n)]
dim(par.eps) <- c(nodes_sm, t_n)
params3 <- list("beta" = opt1$par[1:3], "log_kappa" = out2$par.fixed[1], 
                "log_tau_O" = out2$par.fixed[2], "log_tau_E" = out2$par.fixed[3],
                "omega" = par.omg, "epsilon" = par.eps)

#Make AD objective function
obj3 <- MakeADFun(data = data_sm, parameters = params3, random = random, DLL="Bias")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
out3 <- sdreport(obj3)

#Parameter values
beta <- out3$par.fixed[1:3]
omg <- out3$par.random[1:nodes_sm]
eps <- out3$par.random[(nodes_sm+1):(nodes_sm + nodes_sm*t_n)]
tauO <- exp(out3$par.fixed[5])
tauE <- exp(out3$par.fixed[6])
dim(eps) <- c(nodes_sm, t_n)

omg <- omg/tauO
eps <- eps/tauE

#Predicted values
intensity <- raster::raster(domain_sm, resolution = c(1, 1))
Grid <- coordinates(intensity)
GridA <- inla.spde.make.A(mesh_sm, loc = Grid)

df <- data.frame(cbind(t(rep(NA, t_n + 1))))
colnames(df) <- c("value", "x", "y", "yr")

year <- c("2016", "2017", "2018")

for(t in 1:t_n){
  Pred <- exp(beta[t] + GridA %*% omg + GridA %*% eps[,t])
  raster::values(intensity) <- as.vector(Pred)
  intensity <- as(intensity, "SpatialPixelsDataFrame")
  intensity <- st_as_sf(intensity)
  intensity <- intensity[domain_sm,]
  intensity <- as.data.frame(as(intensity, "Spatial"))
  intensity$yr <- year[t]
  colnames(intensity) <- c("value", "x", "y", "yr")
  df <- rbind(df, intensity)
  rm(intensity)
  intensity <- raster::raster(domain_sm, resolution = c(1, 1))
}

df <- df[-1,]

Figure3 <- ggplot(data = df) + 
  geom_sf(data = Data %>% filter(period == 3), size = 0.25, col = "white") + 
  geom_tile(data = df, aes(x = x, y = y, fill = value), alpha = 0.75) +
  #geom_sf(data = domain, alpha = 0, linetype = "dashed", lwd = 1) +
  scale_fill_gradient2(low = "#EDEDED", mid = "#FFCCCC", high = "#660000",
                       labels = c(0,0.05,0.1), breaks = c(0,0.05,0.1),
                       limits=c(0, 0.1), oob = scales::squish) +
  facet_wrap(~yr) +
  theme_bw()

Figure3

EarlySummer <- out3

#-Save Output-#

Bias <- list(EarlySpring, LateSpring, EarlySummer)

save(Bias, file = "iNatBias.Rdata")
