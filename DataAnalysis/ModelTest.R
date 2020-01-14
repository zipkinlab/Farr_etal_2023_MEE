#------------------------#
#-Preliminary model test-#
#------------------------#

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

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

domain_sp <- st_transform(domain_sp, crs = st_crs(prj))
domain_sm <- st_transform(domain_sm, crs = st_crs(prj))

#-Load Monarch data-#

load(file = "~/Monarchs/DataFormatting/FormattedData.Rdata")

#-Load Node data-#
load(file = '~/Monarchs/DataFormatting/FormattedNodeData.Rdata')

#-Visualize-#

ggplot() + 
  geom_sf(data = domain_sp, alpha = 0.5) +
  geom_sf(data = domain_sm, alpha = 0.5) +
  geom_sf(data = Data, aes(col = program)) +
  facet_wrap(~yr + period) +
  theme_bw()

#-Build mesh-#
#Spring
domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))
#Summer
domain_sm_seg <- inla.sp2segment(as(domain_sm, "Spatial"))

sp.mesh <- inla.mesh.2d(boundary=domain_sp_seg,
                                       max.edge=c(36, 108),
                                       cutoff=10)

su.mesh <- inla.mesh.2d(boundary=domain_sm_seg,
                          max.edge=c(36, 108),
                          cutoff=25)

#Create SPDE objects (Matern covariance)
sp.spde <- inla.spde2.matern(sp.mesh)

su.spde <- inla.spde2.matern(su.mesh)

#-Compile data-#

#Indices
t_n <- Data %>% summarize(out = n_distinct(yr)) %>% select(out) %>% .$out
p_n <- Data %>% summarize(out = n_distinct(period)) %>% select(out) %>% .$out

sp.nodes <- sp.mesh$n
su.nodes <- su.mesh$n
nodes <- as.integer(t_n*(sp.nodes*2 + su.nodes))
tp_k <- as.integer(fct_cross(as.factor(NodesDF$yr), as.factor(NodesDF$period))) - 1

nobs <- dim(Data %>% filter(type == "obs"))[1]
nobs_sp1 <- dim(Data %>% filter(type == "obs" & period == 1))[1]
nobs_sp2 <- dim(Data %>% filter(type == "obs" & period == 2))[1]
nobs_su <- dim(Data %>% filter(type == "obs" & period == 3))[1]
ncount <- dim(Data %>% filter(type == "count"))[1]
ncount_sp1 <- dim(Data %>% filter(type == "count" & period == 1))[1]
ncount_sp2 <- dim(Data %>% filter(type == "count" & period == 2))[1]
ncount_su <- dim(Data %>% filter(type == "count" & period == 3))[1]
tp_i <- as.integer(fct_cross(as.factor(Data$yr), as.factor(Data$period))) - 1

#Create dual mesh for weights
sp.dmesh <- book.mesh.dual(sp.mesh)
sp.dmesh <- st_as_sf(sp.dmesh)
st_crs(sp.dmesh) <- st_crs(domain_sp)

su.dmesh <- book.mesh.dual(su.mesh)
su.dmesh <- st_as_sf(su.dmesh)
st_crs(su.dmesh) <- st_crs(domain_sm)

#Weights of each node (area of dual mesh)
sp.weight <- sapply(1:dim(sp.dmesh)[1], function(i) {
  if (st_intersects(sp.dmesh[i,], domain_sp, sparse = FALSE))
    return(st_area(st_intersection(sp.dmesh[i,], domain_sp)))
  else return(0)
})

su.weight <- sapply(1:dim(su.dmesh)[1], function(i) {
  if (st_intersects(su.dmesh[i,], domain_sm, sparse = FALSE))
    return(st_area(st_intersection(su.dmesh[i,], domain_sm)))
  else return(0)
})

weight <- rep(c(sp.weight, sp.weight, su.weight), t_n)# * 1e6

#Projection matrix of observations to nodes
A1 <- as(inla.spde.make.A(sp.mesh, st_coordinates(Data %>% filter(period == 1))), "dgTMatrix")
A2 <- as(inla.spde.make.A(sp.mesh, st_coordinates(Data %>% filter(period == 2))), "dgTMatrix")
A3 <- as(inla.spde.make.A(su.mesh, st_coordinates(Data %>% filter(period == 3))), "dgTMatrix")

#Counts
counts <- as.numeric(Data %>% filter(type == "count") %>% select(count) %>% .$count)

#Area/effort
area <- as.numeric(Data %>% filter(type == "count") %>% select(effort) %>% .$effort)
area <- area*6*1000*1e-6

#Observation covariates
NDVI <- Data$NDVI
Bias <- as.numeric(Data %>% filter(program == "iNat") %>% select(Bias) %>% .$Bias)

#Node covariates
nNDVI <- NodesDF$NDVI
nBias <- NodesDF$Bias

#Scale covaraites
tmp1 <- (NDVI - mean(c(NDVI, nNDVI)))/sd(c(NDVI, nNDVI))
tmp2 <- (Bias - mean(c(Bias, nBias)))/sd(c(Bias, nBias))

tmp3 <- (nNDVI - mean(c(NDVI, nNDVI)))/sd(c(NDVI, nNDVI))
tmp4 <- (nBias - mean(c(Bias, nBias)))/sd(c(Bias, nBias))

NDVI <- tmp1
Bias <- tmp2
nNDVI <- tmp3
nBias <- tmp4

#Format node covariates
#dim(nNDVI) <- c(nodes, t_n)
#dim(nBias) <- c(nodes, t_n)

#-TMB-#

#PHASE 1: Fit intensity function intercepts

#TMB data
data <- list("nodes" = nodes, "nobs" = nobs, "ncount" = ncount, "counts" = counts,
             "t_n" = t_n, "p_n" = p_n, "tp_k" = tp_k, "tp_i" = tp_i, "weight" = weight, "area" = area,
             #"A" = A, "spde" = spde$param.inla[c("M0","M1","M2")],
             "nBias" = nBias, "nCov" = nNDVI, "Bias" = Bias, "Cov" = NDVI)

#Compile TMB code (only once)
#compile("./DataAnalysis/IM_full.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/IM_full"))

#Parameters to estimate
# params1 <- list("beta0" = rep(0, t_n), "beta1" = 0, "alpha0" = 0, "alpha1" = 0,
#                 "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
#                 "omega" = rep(0, nodes), "epsilon" = matrix(0, nrow = nodes, ncol = t_n))

params1 <- list("beta0" = 0, "beta1" = 0, 
                "delta1" = 0, "delta2" = 0,
                "gamma1" = 0, "gamma2" = 0,
                "alpha0" = 0, "alpha1" = 0)

# map.eps <- factor(rep(NA, nodes*t_n))
# dim(map.eps) <- c(nodes, t_n)
# map1 <- list(
#              "alpha0" = as.factor(NA), 
#              "alpha1" = as.factor(NA), 
#              #"beta0" = factor(rep(NA, t_n)), 
#              "beta1" = as.factor(NA),
#              "log_kappa" = as.factor(NA), 
#              "log_tau_O" = as.factor(NA), 
#              "log_tau_E" = as.factor(NA),
#              "omega" = factor(rep(NA, nodes)), 
#              "epsilon" = map.eps
#              )

#Make AD objective function
obj1 <- MakeADFun(data = data, parameters = params1, DLL="IM_full")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

out1 <- sdreport(obj1)
summary(out1)

#PHASE 2: Fit random effects

#TMB data
data2 <- list("nodes_sp" = sp.nodes, "nodes_su" = su.nodes, 
             "nobs_sp1" = nobs_sp1, "nobs_sp2" =nobs_sp2, "nobs_su" = nobs_su,
             "ncount" = ncount, "ncount_sp1" = ncount_sp1, "ncount_sp2" = ncount_sp2, 
             "ncount_su" = ncount_su, "counts" = counts,
             "t_n" = t_n, "p_n" = p_n, "tp_k" = tp_k, "tp_i" = tp_i, "weight" = weight, "area" = area,
             "A1" = A1, "A2" = A2, "A3" = A3,
             "spde_sp" = sp.spde$param.inla[c("M0","M1","M2")], "spde_su" = su.spde$param.inla[c("M0","M1","M2")],
             "nBias" = nBias, "nCov" = nNDVI, "Bias" = Bias, "Cov" = NDVI)

#Compile TMB code (only once)
#compile("./DataAnalysis/IM_random_separate.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/IM_random_separate"))

#Parameters to estimate
params2 <- list("beta0" = opt1$par[1], "beta1" = opt1$par[2],
                "delta1" = opt1$par[3], "delta2" = opt1$par[4],
                "gamma1" = opt1$par[5], "gamma2" = opt1$par[6],
                "alpha0" = opt1$par[7], "alpha1" = opt1$par[8],
                "log_kappa_sp1" = 0, "log_tau_sp1" = 1, "omega_sp1" = rep(0, sp.nodes),
                "log_kappa_sp2" = 0, "log_tau_sp2" = 1, "omega_sp2" = rep(0, sp.nodes),
                "log_kappa_su" = 0, "log_tau_su" = 1, "omega_su" = rep(0, su.nodes))

map2 <- list(
             "alpha0" = as.factor(NA), 
             "alpha1" = as.factor(NA),
             "beta0" = as.factor(NA), 
             "beta1" = as.factor(NA),
             "delta1" =  as.factor(NA),
             "delta2" = as.factor(NA),
             "gamma1" = as.factor(NA),
             "gamma2" = as.factor(NA)
             #"log_kappa_sp1" = as.factor(NA),
             #"log_tau_sp1" = as.factor(NA),
             #"omega_sp1" = factor(rep(NA, sp.nodes)),
             #"log_kappa_sp2" = as.factor(NA),
             #"log_tau_sp2" = as.factor(NA),
             #"omega_sp2" = factor(rep(NA, sp.nodes)),
             #"log_kappa_su" = as.factor(NA),
             #"log_tau_su" = as.factor(NA),
             #"omega_su" = factor(rep(NA, su.nodes))
             )

#Random effects
random <- c("omega_sp1", "omega_sp2", "omega_su")

#Make AD objective function
obj2 <- MakeADFun(data = data2, parameters = params2, random = random, map = map2, DLL="IM_random_separate")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
out2 <- sdreport(obj2)

#PHASE 3: Refit fixed parameters (not related to random effects)

#Parameters to estimate
par.omg_sp1 <- out2$par.random[1:sp.nodes]
par.omg_sp2 <- out2$par.random[(sp.nodes+1):(sp.nodes*2)]
par.omg_su <- out2$par.random[(sp.nodes*2+1):(sp.nodes*2+su.nodes)]
params3 <- list("beta0" = opt1$par[1], "beta1" = opt1$par[2],
                "delta1" = opt1$par[3], "delta2" = opt1$par[4],
                "gamma1" = opt1$par[5], "gamma2" = opt1$par[6],
                "alpha0" = opt1$par[7], "alpha1" = opt1$par[8],
                "log_kappa_sp1" = out2$par.fixed[1], "log_tau_sp1" = out2$par.fixed[4], "omega_sp1" = par.omg_sp1,
                "log_kappa_sp2" = out2$par.fixed[2], "log_tau_sp2" = out2$par.fixed[5], "omega_sp2" = par.omg_sp2,
                "log_kappa_su" = out2$par.fixed[3], "log_tau_su" = out2$par.fixed[6], "omega_su" = par.omg_su)

map3 <- list(
  # "alpha0" = as.factor(NA), 
  # "alpha1" = as.factor(NA),
  # "beta0" = as.factor(NA), 
  # "beta1" = as.factor(NA),
  # "delta1" =  as.factor(NA),
  # "delta2" = as.factor(NA),
  # "gamma1" = as.factor(NA),
  # "gamma2" = as.factor(NA)
  "log_kappa_sp1" = as.factor(NA),
  "log_tau_sp1" = as.factor(NA),
  "omega_sp1" = factor(rep(NA, sp.nodes)),
  "log_kappa_sp2" = as.factor(NA),
  "log_tau_sp2" = as.factor(NA),
  "omega_sp2" = factor(rep(NA, sp.nodes)),
  "log_kappa_su" = as.factor(NA),
  "log_tau_su" = as.factor(NA),
  "omega_su" = factor(rep(NA, su.nodes))
)

#Make AD objective function
obj3 <- MakeADFun(data = data2, parameters = params3, random = random, map = map3, DLL="IM_random_separate")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
out3 <- sdreport(obj3)

#-Parameter values-#

beta0 <- out3$par.fixed[1]
beta1 <- out3$par.fixed[2]
delta1 <- out3$par.fixed[3]
delta2 <- out3$par.fixed[4]
gamma1 <- out3$par.fixed[5]
gamma2 <- out3$par.fixed[6]
tau_sp1 <- exp(out2$par.fixed[4])
tau_sp2 <- exp(out2$par.fixed[5])
tau_su <- exp(out2$par.fixed[6])
omg_sp1 <- par.omg_sp1/tau_sp1
omg_sp2 <- par.omg_sp2/tau_sp2
omg_su <- par.omg_su/tau_su

intercept <- NULL
intercept[1] <- beta0 + delta1 + delta2
intercept[2] <- beta0 + delta1
intercept[3] <- beta0 + delta1 + delta2 + gamma1
intercept[4] <- beta0 + delta1 + gamma1
intercept[5] <- beta0 + delta1 + delta2 + gamma1 + gamma2
intercept[6] <- beta0 + delta1 + gamma1 + gamma2
intercept[7] <- beta0
intercept[8] <- beta0 + gamma1
intercept[9] <- beta0 + gamma1 + gamma2

#-Visualization-#

sp.intensity <- raster::raster(domain_sp, resolution = c(1, 1))
sp.Grid <- coordinates(sp.intensity)
sp.GridA <- inla.spde.make.A(sp.mesh, loc = sp.Grid)

su.intensity <- raster::raster(domain_sm, resolution = c(1, 1))
su.Grid <- coordinates(su.intensity)
su.GridA <- inla.spde.make.A(su.mesh, loc = su.Grid)

df <- data.frame(cbind(t(rep(NA, 5))))
colnames(df) <- c("value", "x", "y", "period", "year")

NDVIspring <- raster::raster("~/Monarchs/DataFormatting/NDVIspring.tif")
NDVIsummer <- raster::raster("~/Monarchs/DataFormatting/NDVIsummer.tif")

NDVIspring <- raster::projectRaster(NDVIspring, sp.intensity, method = 'ngb')
NDVIsummer <- raster::projectRaster(NDVIsummer, su.intensity, method = 'ngb')

sp.NDVI <- matrix(NA, nrow = dim(sp.Grid)[1], ncol = t_n*2)
su.NDVI <- matrix(NA, nrow = dim(su.Grid)[1], ncol = t_n)

sp.NDVI[,1] <- 0.0001*raster::values(raster::overlay(NDVIspring[[1]], NDVIspring[[2]], fun = mean))
sp.NDVI[,2] <- 0.0001*raster::values(raster::overlay(NDVIspring[[3]], NDVIspring[[4]], fun = mean))
sp.NDVI[,3] <- 0.0001*raster::values(raster::overlay(NDVIspring[[5]], NDVIspring[[6]], fun = mean))
sp.NDVI[,4] <- 0.0001*raster::values(raster::overlay(NDVIspring[[7]], NDVIspring[[8]], fun = mean))
sp.NDVI[,5] <- 0.0001*raster::values(raster::overlay(NDVIspring[[9]], NDVIspring[[10]], fun = mean))
sp.NDVI[,6] <- 0.0001*raster::values(raster::overlay(NDVIspring[[11]], NDVIspring[[12]], fun = mean))

su.NDVI[,1] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[1]], NDVIsummer[[2]], NDVIsummer[[3]], fun = mean))
su.NDVI[,2] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[4]], NDVIsummer[[5]], NDVIsummer[[6]], fun = mean))
su.NDVI[,3] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[7]], NDVIsummer[[8]], NDVIsummer[[9]], fun = mean))

sp.NDVI <- (sp.NDVI - mean(sp.NDVI))/sd(sp.NDVI)
su.NDVI <- (su.NDVI - mean(su.NDVI))/sd(su.NDVI)


sp.Pred <- cbind(exp(intercept[1] + beta1 * sp.NDVI[,1] + sp.GridA %*% omg_sp1), 
              exp(intercept[2] + beta1 * sp.NDVI[,2] + sp.GridA %*% omg_sp2),
              exp(intercept[3] + beta1 * sp.NDVI[,3] + sp.GridA %*% omg_sp1), 
              exp(intercept[4] + beta1 * sp.NDVI[,4] + sp.GridA %*% omg_sp2),
              exp(intercept[5] + beta1 * sp.NDVI[,5] + sp.GridA %*% omg_sp1), 
              exp(intercept[6] + beta1 * sp.NDVI[,6] + sp.GridA %*% omg_sp2))

su.Pred <- cbind(exp(intercept[7] + beta1 * su.NDVI[,1] + su.GridA %*% omg_su),
              exp(intercept[8] + beta1 * su.NDVI[,2] + su.GridA %*% omg_su),
              exp(intercept[9] + beta1 * su.NDVI[,3] + su.GridA %*% omg_su))


tp <- cbind(c(1,2,1,2,1,2,3,3,3), c(1,1,2,2,3,3,1,2,3))

for(t in 1:6){
  raster::values(sp.intensity) <- as.vector(sp.Pred[,t])
  sp.intensity <- as(sp.intensity, "SpatialPixelsDataFrame")
  sp.intensity <- st_as_sf(sp.intensity)
  sp.intensity <- sp.intensity[domain_sp,]
  sp.intensity <- as.data.frame(as(sp.intensity, "Spatial"))
  sp.intensity$period <- tp[t,1]
  sp.intensity$year <- tp[t,2]
  colnames(sp.intensity) <- c("value", "x", "y", "period", "year")
  df <- rbind(df, sp.intensity)
  rm(sp.intensity)
  sp.intensity <- raster::raster(domain_sp, resolution = c(1, 1))
}

for(t in 1:3){
  raster::values(su.intensity) <- as.vector(su.Pred[,t])
  su.intensity <- as(su.intensity, "SpatialPixelsDataFrame")
  su.intensity <- st_as_sf(su.intensity)
  su.intensity <- su.intensity[domain_sm,]
  su.intensity <- as.data.frame(as(su.intensity, "Spatial"))
  su.intensity$period <- tp[t+6,1]
  su.intensity$year <- tp[t+6,2]
  colnames(su.intensity) <- c("value", "x", "y", "period", "year")
  df <- rbind(df, su.intensity)
  rm(su.intensity)
  su.intensity <- raster::raster(domain_sm, resolution = c(1, 1))
}

df <- df[-1,]

df1 <- st_as_sf(df, coords = c("x", "y"), crs = prj)

pal3 <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(100)

Figure1 <- ggplot(data = df) + 
  #geom_sf(data = Data, size = 0.25) + 
  geom_raster(data = df, aes(x = x, y = y, fill = value), alpha = 0.75) +
  #geom_sf(data = domain, alpha = 0, linetype = "dashed", lwd = 1) +
  # scale_fill_gradient2(low = "#EDEDED", mid = "#FFCCCC", high = "#660000",
  #                      labels = c(0,50,312,625), breaks = c(0,50,312,625),
  #                      limits=c(0, 625), oob = scales::squish) +
  scale_fill_gradientn(colors = pal3) +
  facet_wrap(~period + year) +
  #coord_sf(crs = st_crs("+proj=utm +zone=14 +ellps=GRS80 +units=km +no_defs")) + 
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

Figure1

#-Spatial random effect-#

re.raster <- raster::raster(domain_sp, resolution = c(1, 1))

for(t in 1:3){
  raster::values(su.intensity) <- as.vector(su.Pred[,t])
  su.intensity <- as(su.intensity, "SpatialPixelsDataFrame")
  su.intensity <- st_as_sf(su.intensity)
  su.intensity <- su.intensity[domain_sm,]
  su.intensity <- as.data.frame(as(su.intensity, "Spatial"))
  su.intensity$period <- tp[t+6,1]
  su.intensity$year <- tp[t+6,2]
  colnames(su.intensity) <- c("value", "x", "y", "period", "year")
  df <- rbind(df, su.intensity)
  rm(su.intensity)
  su.intensity <- raster::raster(domain_sm, resolution = c(1, 1))
}

df <- df[-1,]

df1 <- st_as_sf(df, coords = c("x", "y"), crs = prj)

pal3 <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(100)

Figure1 <- ggplot(data = df) + 
  #geom_sf(data = Data, size = 0.25) + 
  geom_raster(data = df, aes(x = x, y = y, fill = value), alpha = 0.75) +
  #geom_sf(data = domain, alpha = 0, linetype = "dashed", lwd = 1) +
  # scale_fill_gradient2(low = "#EDEDED", mid = "#FFCCCC", high = "#660000",
  #                      labels = c(0,50,312,625), breaks = c(0,50,312,625),
  #                      limits=c(0, 625), oob = scales::squish) +
  scale_fill_gradientn(colors = pal3) +
  facet_wrap(~period + year) +
  #coord_sf(crs = st_crs("+proj=utm +zone=14 +ellps=GRS80 +units=km +no_defs")) + 
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

Figure1

