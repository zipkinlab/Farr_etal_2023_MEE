#-------------------------------------------#
#----Spatiotemporal integrated model--------#
#----A case study on monarch butterflies----#
#-------------------------------------------#

#-----------#
#-Libraries-#
#-----------#

library(tidyverse)
library(sf)
library(INLA)
library(TMB)

#Function to create dual mesh (original code from Krainski et al. 2018)
source("~/Monarchs/DataAnalysis/rDualMesh.R")

#-----------#
#-Load data-#
#-----------#

#Spring domain
domain_sp <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Spring.shp")

#Summer domain
domain_su <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Summer.shp")

#Monarch data
load(file = "~/Monarchs/DataFormatting/FormattedData_2022.Rdata")

#Mesh node data
load(file = '~/Monarchs/DataFormatting/FormattedNodeData.Rdata')

#-------------#
#-Format data-#
#-------------#

#Geographic projection
prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

#Set projection
domain_sp <- st_transform(domain_sp, crs = st_crs(prj))
domain_su <- st_transform(domain_su, crs = st_crs(prj))

#Spring domain boundary
domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))

#Summer domain boundary
domain_su_seg <- inla.sp2segment(as(domain_su, "Spatial"))


#Construct triangulated mesh
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

#Create SPDE objects (Matern covariance structure)
spde1 <- inla.spde2.matern(mesh[[1]])
spde2 <- inla.spde2.matern(mesh[[2]])
spde3 <- inla.spde2.matern(mesh[[3]])

#Create dual mesh for weights (see Krainski et al. 2018)
dmesh1 <- book.mesh.dual(mesh[[1]])
dmesh1 <- st_as_sf(dmesh1)
st_crs(dmesh1) <- st_crs(domain_sp)

dmesh2 <- book.mesh.dual(mesh[[2]])
dmesh2 <- st_as_sf(dmesh2)
st_crs(dmesh2) <- st_crs(domain_sp)

dmesh3 <- book.mesh.dual(mesh[[3]])
dmesh3 <- st_as_sf(dmesh3)
st_crs(dmesh3) <- st_crs(domain_su)

#Weights of each node (area of dual mesh)
weight1 <- sapply(1:dim(dmesh1)[1], function(i) {
  if (st_intersects(dmesh1[i,], domain_sp, sparse = FALSE))
    return(st_area(st_intersection(dmesh1[i,], domain_sp)))
  else return(0)
})

weight2 <- sapply(1:dim(dmesh2)[1], function(i) {
  if (st_intersects(dmesh2[i,], domain_sp, sparse = FALSE))
    return(st_area(st_intersection(dmesh2[i,], domain_sp)))
  else return(0)
})

weight3 <- sapply(1:dim(dmesh3)[1], function(i) {
  if (st_intersects(dmesh3[i,], domain_su, sparse = FALSE))
    return(st_area(st_intersection(dmesh3[i,], domain_su)))
  else return(0)
})

#Projection matrix of observations to nodes
A1 <- as(inla.spde.make.A(mesh[[1]], st_coordinates(Data %>% filter(period == 1))), "dgTMatrix")
A2 <- as(inla.spde.make.A(mesh[[2]], st_coordinates(Data %>% filter(period == 2))), "dgTMatrix")
A3 <- as(inla.spde.make.A(mesh[[3]], st_coordinates(Data %>% filter(period == 3))), "dgTMatrix")

#--------------#
#-Compile data-#
#--------------#

#Number of years
t_n <- Data %>% summarize(out = n_distinct(yr)) %>% select(out) %>% .$out

#Number of stages
p_n <- Data %>% summarize(out = n_distinct(period)) %>% select(out) %>% .$out

#Number of nodes in early spring
nodes1 <- mesh[[1]]$n

#Number of nodes in late spring
nodes2 <- mesh[[2]]$n

#Number of nodes in early summer
nodes3 <- mesh[[3]]$n

#Number of nodes across stages and years
nodes <- as.integer(t_n*(nodes1 + nodes2 + nodes3))

#MOVE TO NODE FORMATTING
NodeDF$period <- as.factor(rep(1:3, c(nodes1 * 3, nodes2 * 3, nodes3 * 3)))

#Year-stage identifier for nodes
tp_k <- as.integer(fct_cross(NodeDF$period, NodeDF$Year)) - 1

#Stage identifier for nodes
p_k <- as.integer(NodeDF$period) - 1

#Number of presence-only observations
nobs <- dim(Data %>% filter(type == "obs"))[1]

#Number of presence-only observations in early spring
nobs1 <- dim(Data %>% filter(type == "obs" & period == 1))[1]

#Number of presence-only observations in late spring
nobs2 <- dim(Data %>% filter(type == "obs" & period == 2))[1]

#Number of presence-only observations in early summer
nobs3 <- dim(Data %>% filter(type == "obs" & period == 3))[1]

#Number of single-visit sites
ncount <- dim(Data %>% filter(type == "count"))[1]

#Number of single-visit sites in early spring
ncount1 <- dim(Data %>% filter(type == "count" & period == 1))[1]

#Number of single-visit sites in late spring
ncount2 <- dim(Data %>% filter(type == "count" & period == 2))[1]

#Number of single-visit sites in early summer
ncount3 <- dim(Data %>% filter(type == "count" & period == 3))[1]

#Single-visit counts @ each site
counts <- as.numeric(Data %>% filter(type == "count") %>% select(count) %>% .$count)

#Year-stage identifier for data
tp_i <- as.integer(fct_cross(as.factor(Data$period), as.factor(Data$yr))) - 1 #CHECK THIS

#Stage identifier for data
p_i <- as.integer(as.factor(Data$period)) - 1

#Weight of each node @ 100 m^2
weight <- rep(c(weight1, weight2, weight3), t_n) * 1e4

#Area/effort offset
area <- as.numeric(Data %>% filter(type == "count") %>% select(effort) %>% .$effort)

#Covert to 100 m^2
area <- area*6*1000*1e-2

#Observation covariates
Data <- Data %>% mutate(daysaccum = as.numeric(as.Date(paste0(yr, "-", mo_day)) - 
                                                 as.Date(paste0(yr, ifelse(period == 1, "-03-01",
                                                                           ifelse(period == 2, "-03-22", "-04-26"))))) + 1,
                        gdd.avg = gdd2/daysaccum)

#Normalized Difference Vegetation Index @ data
NDVI <- Data$NDVI

#Averaged daily growing degree days @ data
#GDD <- Data$gdd.avg
GDD <- Data$gdd2

#Number of other butterfly observations @ data
Bias <- as.numeric(Data %>% filter(type == "obs") %>% select(Bias) %>% .$Bias)

#Population density @ data
PopD <- as.numeric(Data %>% filter(type == "obs") %>% select(PopD) %>% .$PopD)


#Node covariates
NodeDF <- NodeDF %>% mutate(daysaccum = ifelse(period == 1, 35, 42),
                            gdd.avg = gdd2/daysaccum)

#Normalized Difference Vegetation Index @ nodes
nNDVI <- NodeDF$NDVI

#Averaged daily growing degree days @ nodes summed to 14 day window
nGDD <- NodeDF$gdd.avg * 14
#nGDD <- NodeDF$gdd2

#Number of other butterfly observations @ nodes
nBias <- NodeDF$Bias

#Population density @ nodes
nPopD <- NodeDF$PopD

#Scale covaraites for both data & nodes
tmp1 <- (NDVI - mean(c(NDVI, nNDVI)))/sd(c(NDVI, nNDVI))
tmp2 <- (GDD - mean(c(GDD, nGDD)))/sd(c(GDD, nGDD))
tmp3 <- (Bias - mean(c(Bias, nBias)))/sd(c(Bias, nBias))
tmp4 <- (PopD - mean(c(PopD, nPopD)))/sd(c(PopD, nPopD))

tmp5 <- (nNDVI - mean(c(NDVI, nNDVI)))/sd(c(NDVI, nNDVI))
tmp6 <- (nGDD - mean(c(GDD, nGDD)))/sd(c(GDD, nGDD))
tmp7 <- (nBias - mean(c(Bias, nBias)))/sd(c(Bias, nBias))
tmp8 <- (nPopD - mean(c(PopD, nPopD)))/sd(c(PopD, nPopD))

NDVI <- tmp1
GDD <- tmp2
Bias <- tmp3
PopD <- tmp4

nNDVI <- tmp5
nGDD <- tmp6
nBias <- tmp7
nPopD <- tmp8

#--------------#
#-Optimize nll-#
#--------------#

#Compile TMB code (only once)
# compile("./DataAnalysis/IntegratedModel_Final.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/IntegratedModel_Final"))

#PHASE 1: Fit fixed effects

#Compile data for TMB
data <- list("nodes1" = nodes1, "nodes2" = nodes2, "nodes3" = nodes3,
             "nobs1" = nobs1, "nobs2" = nobs2, "nobs3" = nobs3, 
             "ncount" = ncount, "ncount1" = ncount1, "ncount2" = ncount2, "ncount3" = ncount3,
             "counts" = counts,
             "t_n" = t_n, "p_n" = p_n, "tp_k" = tp_k, "tp_i" = tp_i, "p_i" = p_i, "p_k" = p_k,
             "weight" = weight, "area" = area,
             "A1" = A1, "spde1" = spde1$param.inla[c("M0","M1","M2")],             
             "A2" = A2, "spde2" = spde2$param.inla[c("M0","M1","M2")],
             "A3" = A3, "spde3" = spde3$param.inla[c("M0","M1","M2")],
             "nBias" = nBias, "nPopD" = nPopD,
             "nNDVI" = nNDVI, "nGDD" = nGDD, 
             "Bias" = Bias, "PopD" = PopD,
             "NDVI" = NDVI, "GDD" = GDD)

#Parameters to estimate
params <- list("beta0" = 0, 
               "beta1" = rep(0,3), "beta2" = rep(0,3),
               "beta3" = rep(0,3), "beta4" = rep(0,3),
               "delta1" = 0, "delta2" = 0,
               "gamma1" = 0, "gamma2" = 0,
               "alpha0" = 0, "alpha1" = rep(0,9), "alpha2" = rep(0,9),
               "log_kappa" = 0, "log_tau" = 1,
               "omega1" = rep(0, nodes1),
               "omega2" = rep(0, nodes2),
               "omega3" = rep(0, nodes3))

#Hold random effects constant
map <- list(
  "log_kappa" = as.factor(NA),
  "log_tau" = as.factor(NA),
  "omega1" = factor(rep(NA, nodes1)),
  "omega2" = factor(rep(NA, nodes2)),
  "omega3" = factor(rep(NA, nodes3))
)

#Random effects
random <- c("omega1", "omega2", "omega3")

#Make AD objective function
obj1 <- MakeADFun(data = data, parameters = params, random = random, map = map, DLL="IntegratedModel_Final")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#AIC
2 * length(obj1$par) - 2 * (-1 * sum(obj1$report()$jnll_comp[3:12]))

#Output
out1 <- sdreport(obj1)
summary(out1)

#Extract parameter values
beta0 <- out1$par.fixed[1]
beta1 <- out1$par.fixed[2:4]
beta2 <- out1$par.fixed[5:7]
beta3 <- out1$par.fixed[8:10]
beta4 <- out1$par.fixed[11:13]
delta1 <- out1$par.fixed[14]
delta2 <- out1$par.fixed[15]
gamma1 <- out1$par.fixed[16] 
gamma2 <- out1$par.fixed[17]

#Expected population densities per stage and year
beta <- NULL

#2016 stage 1 (early spring)
beta[1] <- beta0 + gamma1 + gamma2 + delta1 + delta2

#2016 stage 2 (late spring)
beta[2] <- beta0 + gamma1 + delta1 + delta2 

#2016 stage 3 (early summer)
beta[3] <- beta0 + delta1 + delta2

#2017 stage 1 (early spring)
beta[4] <- beta0 + gamma1 + gamma2 + delta1 

#2017 stage 2 (late spring)
beta[5] <- beta0 + gamma1 + delta1

#2017 stage 3 (early summer)
beta[6] <- beta0 + delta1

#2018 stage 1 (early spring)
beta[7] <- beta0 + gamma1 + gamma2

#2018 stage 2 (late spring)
beta[8] <- beta0 + gamma1

#2018 stage 3 (early summer)
beta[9] <- beta0

#Compile output
Output <- data.frame(year = factor(rep(2016:2018, each = 3)), 
                     period = factor(rep(1:3, 3)),
                     mean.den = beta)

Output$lower.den <- NA
Output$upper.den <- NA

#Profile confidence intervals
Output[1,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1","gamma2","delta1","delta2"))))
Output[2,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1","delta1","delta2"))))
Output[3,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","delta1","delta2"))))
Output[4,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1","gamma2","delta1"))))
Output[5,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1","delta1"))))
Output[6,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","delta1"))))
Output[7,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1","gamma2"))))
Output[8,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0","gamma1"))))
Output[9,4:5] <- confint(tmbprofile(obj1, lincomb = as.numeric(names(obj1$par) %in% c("beta0"))))

#PHASE 2: Fit random effects

#Parameters to estimate
params2 <- list("beta0" = beta0, "beta1" = beta1,
                "beta2" = beta2, "beta3" = beta3,
                "beta4" = beta4,
                "delta1" = delta1, "delta2" = delta2,
                "gamma1" = gamma1, "gamma2" = gamma2,
                "alpha0" = out1$par.fixed[18], 
                "alpha1" = out1$par.fixed[19:27],
                "alpha2" = out1$par.fixed[28:36],
                "log_kappa" = 0, "log_tau" = 1, 
                "omega1" = rep(0, nodes1), "omega2" = rep(0, nodes2), "omega3" = rep(0, nodes3))

#Hold fixed effects constant
map2 <- list(
  "alpha0" = as.factor(NA), 
  "alpha1" = factor(rep(NA, 9)),
  "alpha2" = factor(rep(NA, 9)),
  "beta0" = as.factor(NA), 
  "beta1" = factor(rep(NA, 3)),
  "beta2" = factor(rep(NA, 3)),
  "beta3" = factor(rep(NA, 3)),
  "beta4" = factor(rep(NA, 3)),
  "delta1" =  as.factor(NA),
  "delta2" = as.factor(NA),
  "gamma1" = as.factor(NA),
  "gamma2" = as.factor(NA)
)

#Make AD objective function
obj2 <- MakeADFun(data = data, parameters = params2, map = map2, random = random, DLL="IntegratedModel_Final")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

#Output
out2 <- sdreport(obj2)

#Estimated parameter values
Effect <- data.frame(param = rep(c("beta1", "beta2", "beta3", "beta4"), each = 3),
                     period = factor(rep(1:3, 4)),
                     mean.den = c(beta1, beta2, beta3, beta4))

Effect$lower <- NA
Effect$upper <- NA

Effect[1,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,1,rep(0,34))))
Effect[2,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,1,rep(0,33))))
Effect[3,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,1,rep(0,32))))
Effect[4,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,0,1,rep(0,31))))
Effect[5,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,0,0,1,rep(0,30))))
Effect[6,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,6),1,rep(0,29))))
Effect[7,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,7),1,rep(0,28))))
Effect[8,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,8),1,rep(0,27))))
Effect[9,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,9),1,rep(0,26))))
Effect[10,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,10),1,rep(0,25))))
Effect[11,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,11),1,rep(0,24))))
Effect[12,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,12),1,rep(0,23))))


omg_sp1 <- out2$par.random[1:nodes1]/exp(out2$par.fixed[2])
omg_sp2 <- out2$par.random[(nodes1+1):(nodes1+nodes2)]/exp(out2$par.fixed[2])
omg_su <- out2$par.random[(nodes1+nodes2+1):(nodes1+nodes2+nodes3)]/exp(out2$par.fixed[2])

#-----------------------#
#-Supplemental Material-#
#-----------------------#

#Appendix D Table D.1
write.csv(Output, file = "TableE.1.csv")

#Appendix D Table D.2
write.csv(Effect, file = "TableE.2.csv")

#Fixed effects output



