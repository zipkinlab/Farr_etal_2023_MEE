#---------------------------------------------------#
#----Estimating a log Gaussian Cox process using----#
#----stochastic partial differential equations------#
#----via Template Model Builder---------------------#
#----Create by Matthew T. Farr----------------------#
#---------------------------------------------------#

#-----------#
#-Libraries-#
#-----------#

library(spatstat)
library(RandomFields) #Package is deprecated from CRAN
library(INLA)
library(rgeos)
library(TMB)
library(gridExtra)

#-----------#
#-Functions-#
#-----------#

#Function to create dual mesh (original code from Krainski et al. 2018)
source("rDualMesh.R")

#Logit function
logit <- function(pp) 
{ 
  log(pp) - log(1-pp) 
}

#Inverse logit
expit <- function(eta) 
{
  1/(1+exp(-eta))
}

#-------------------------#
#-Spatio-temporal domains-#
#-------------------------#

#Full spatio-temporal extent
win <- owin(c(0, 3), c(0, 3)) #Window
loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)) #Coordinates

#Number of pixels
npix <- 300
spatstat.options(npixel = npix)

#Create triangulated mesh over domain
mesh <- inla.mesh.2d(loc.domain = loc.d,
                     offset = c(0.3, 1),
                     max.edge = c(0.3, 0.7),
                     cutoff = 0.05)

#X range of mesh
x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = npix)

#Y range of mesh
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = npix)

#Mesh window
Mwin <- owin(c(min(x0), max(x0)), c(min(y0), max(y0)))

#---------------#
#-Simulate LGCP-#
#---------------#

#Parameters of intensity function
beta <- c(NA, NA)
beta[1] <- 6
beta[2] <- runif(1,0.5,1.5)

#Change in intensity
delta <- 1

#Parameters of thinning function
alpha <- c(logit(0.01), 0.85)

#Ecological covariate
RandomFields::RFoptions(spConform=FALSE)
xcov1 <- RFsimulate(model = RMgauss(scale = 1.25), x = x0, y = y0, grid = TRUE)
xcov2 <- xcov1 + RFsimulate(model = RMgauss(var = 0.5, scale = 9.5), x = x0, y = y0, grid = TRUE)
# load("xcov1.Rds")
# load("xcov2.Rds")
x1 <- (xcov1 - mean(c(xcov1, xcov2))/sd(c(xcov1, xcov2)))
x2 <- (xcov2 - mean(c(xcov1, xcov2))/sd(c(xcov1, xcov2)))

#Thinning covariate
pcov <- RFsimulate(model = RMgauss(scale = 1), x = x0, y = y0, grid = TRUE)
#Alternatively load thinning covariate
# load("pcov.Rds")
#Standardize covariate
pcov <- (pcov - min(pcov))/sd(pcov)
#Range of Matern covaraince
range <- 1.2
#Scale of Matern covariance
kappa <- sqrt(8)/range
#Variance of Matern covaraince
sigma2 <- 0.2
#Precision of Matern covariance
tau <- sqrt(1/(sigma2*4*pi*kappa*kappa))
# tau <- exp(0.5 * log(1/(4*pi)) - log(sqrt(sigma2)) - log(kappa))

#Smoothness of Matern covaraince
nu <- 1

#Simulate spatial covariance

zcov <- RFsimulate(model = RMgauss(var = sigma2, scale = 1/kappa), x = x0, y = y0, grid = TRUE)
# zcov <- RFsimulate(model = RMmatern(var = sigma2, scale = 1/kappa, nu = nu), x = x0, y = y0, grid = TRUE)
# load("zcov.Rds")

#Intensity function
X1 <- as.im(x1, W = Mwin)[win]
X2 <- as.im(x2, W = Mwin)[win]
z <- as.im(zcov, W = Mwin)[win]

lambda1 <- exp(beta[1] + beta[2] * X1 + z)
lambda2 <- exp(beta[1] + delta[1] + beta[2] * X2 + z)

#Format intensity function
lambda1 <- as.im(lambda1, W=win)
lambda2 <- as.im(lambda2, W=win)

#Simulate LGCP
PP1 <- rpoispp(lambda1)[win]
PP2 <- rpoispp(lambda2)[win]

#---------------------------------#
#-Simulate opportunistic sampling-#
#---------------------------------#

#Thinning function
thin <- alpha[1] + alpha[2] * pcov
thin <- expit(thin)

#Format thinning function
thin <- as.im(thin, W=Mwin)[win]

#Latent observations
Z1 <- rbinom(n = PP1$n, size = 1, prob = thin[PP1])
PO1 <- PP1[Z1==1]
# PO1 <- PO1[win1]

Z2 <- rbinom(n = PP2$n, size = 1, prob = thin[PP2])
PO2 <- PP2[Z2==1]
# PO2 <- PO2[win2]

#-----------------#
#-Simulate counts-#
#-----------------#

#Number of counts
# ncount <- 25
ncount <- 100

#Count locations
# u.loc <- expand.grid(seq(1 + 0.175, 2.75 - 0.175, 0.35), seq(0.25 + 0.175, 2 - 0.175, 0.35))

# u.loc <- expand.grid(seq(1 + 0.0875, 2.75 - 0.0875, 0.175), seq(0.25 + 0.0875, 2 - 0.0875, 0.175))

u.loc <- expand.grid(seq(0.15, 2.85, 0.3), seq(0.15, 2.85, 0.3))

#Latent abundance
N <- NULL

#Counts
count <- rep(0, ncount)

#Covariate value @ location
Cov <- NULL

#Unit area
unit.A <- NULL

for(j in 1:ncount){
  unit <- disc(radius = 0.15, centre = as.numeric(u.loc[j,]))
  # unit <- disc(radius = 0.0875, centre = as.numeric(u.loc[j,]))
  count[j] <- PP1[unit]$n
  # count[j] <- rbinom(1, count[j], 0.9)
  Cov[j] <- mean(as.im(x1, W=Mwin)[unit])
  unit.A[j] <- area(unit)
}

#Count dataframe
countdf <- data.frame(count, Cov, unit.A, u.loc)

#---------------#
#-SPDE approach-#
#---------------#

#Update mesh
# mesh <- inla.mesh.2d(boundary = loc.d, max.edge = 0.3, cutoff = 0.05)

#Create SPDE objects (Matern covariance)
spde <- inla.spde2.matern(mesh)

#Create dual mesh for weights
dmesh <- book.mesh.dual(mesh)

#Format location domain
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

#Weights of each node (area of dual mesh)
weight <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

rm(dmesh)

#--------------#
#-Compile Data-#
#--------------#

#Node index
nodes <- mesh$n

#Observation index
nobs <- sum(PO1$n, PO2$n)
ncount <- dim(countdf)[1]

#Period index
t_i <- rep(0:1, c(PO1$n, PO2$n))
t_n <- 2

#Counts
counts <- countdf$count

#Area @ each count location
Area <- countdf$unit.A

#Covaraites @ nodes
nloc <- as.ppp(mesh$loc[,1:2], W = Mwin)
nCov <- matrix(NA, nrow = nodes, ncol = t_n)
nCov[,1] <- as.im(x1, W = Mwin)[nloc]
nCov[,2] <- as.im(x2, W = Mwin)[nloc]
nBias <- as.im(pcov, W = Mwin)[nloc]

#Covariates @ obs
Cov <- c(as.im(x1, W = Mwin)[PO1],
         as.im(x2, W = Mwin)[PO2],
         countdf$Cov)

Bias <- c(as.im(pcov, W = Mwin)[PO1],
          as.im(pcov, W = Mwin)[PO2])

#Location of obsevations
locxy <- rbind(cbind(PO1$x, PO1$y)[,2:1],
               cbind(PO2$x, PO2$y)[,2:1],
               as.matrix(countdf[,4:5]))

#Projection matrix of observations
A <- as(inla.spde.make.A(mesh, locxy), "dgTMatrix")

#------------#
#-Prediction-#
#------------#

Grid <- as.matrix(expand.grid(seq(0.01,2.99,0.01), seq(0.01,2.99,0.01)))
npred <- dim(Grid)[1]
Apred <- as(inla.spde.make.A(mesh, Grid), "dgTMatrix")
predX <- matrix(NA, ncol = t_n, nrow = npred)
predX[,1] <- interp.im(X1, x = Grid[,1], y = Grid[,2])
predX[,2] <- interp.im(X2, x = Grid[,1], y = Grid[,2])

#------------------#
#-Integrated Model-#
#------------------#

#Compile TMB code (only once)
# compile("Simulation.cpp")

#Load TMB code
dyn.load(dynlib("Simulation"))

#Compile data
data <- list("nodes" = nodes, "nobs" = nobs, "ncount" = ncount,
             "t_i" = t_i, "t_n" = t_n,
             "weight" = weight, "area" = Area, "A" = A, "counts" = counts,
             "nBias" = nBias, "nCov" = nCov, "Bias" = Bias, "Cov" = Cov,
             "spde" = spde$param.inla[c("M0","M1","M2")],
             "npred" = npred, "Apred" = Apred, "predX" = predX)

#PHASE 1: Fit fixed effects
#Parameters to estimate
params <- list("beta0" = 0, "beta1" = 1, "delta1" = 1,
               "alpha0" = logit(0.01), "alpha1" = 0, 
               "log_kappa" = log(kappa), "log_tau_O" = 1,
               "omega" = rep(0, mesh$n))

#Define random effects
random = c("omega")

#Map
map <- list(
  "log_kappa" = as.factor(NA), 
  "log_tau_O" = as.factor(NA), 
  "omega" = factor(rep(NA, mesh$n)))


#Make AD objective function
obj1 <- MakeADFun(data = data, parameters = params, random = random, map = map, DLL="Simulation")

#Trace parameters
#obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#Calculate standard deviations
out1 <- sdreport(obj1)

#PHASE 2: Fit random effects
params <- list("beta0" = opt1$par[1], "beta1" = opt1$par[2], "delta1" = opt1$par[3],
               "alpha0" = opt1$par[4], "alpha1" = opt1$par[5],
               "log_kappa" = log(kappa), "log_tau_O" = log(tau),
               "omega" = rep(0, mesh$n))

#Map
map <- list(
  "beta0" = as.factor(NA),
  "beta1" = as.factor(NA),
  "delta1" = as.factor(NA),
  "alpha0" = as.factor(NA),
  "alpha1" = as.factor(NA))


#Make AD objective function
obj2 <- MakeADFun(data = data, parameters = params, random = random, map = map, DLL="Simulation")

#Trace parameters
#obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

#Calculate standard deviations
out2 <- sdreport(obj2)

#Store output

output <- as.data.frame(round(cbind(Truth = c(beta, delta, alpha, PP1$n, PP2$n, kappa, tau),  
                                    rbind(summary(out1)[!grepl("omega|sigma_O|\\bkappa\\b|\\btau_O\\b", rownames(summary(out1))),], 
                                          summary(out2)[grepl("\\bkappa\\b|\\btau_O\\b", rownames(summary(out2))),])),
                 digits = 2))
colnames(output)[2:3] <- c("IM Estimate", "IM Std. Error")

#---------------------#
#-Presence-only model-#
#---------------------#

#Compile TMB code (only once)
# compile("Simulation_PO.cpp")

#Load TMB code
dyn.load(dynlib("Simulation_PO"))

#Compile data
data <- list("nodes" = nodes, "nobs" = nobs, "t_i" = t_i, "t_n" = t_n, 
             "weight" = weight, "A" = A[1:nobs,],
             "nBias" = nBias, "nCov" = nCov, "Bias" = Bias, "Cov" = Cov[1:nobs],
             "spde" = spde$param.inla[c("M0","M1","M2")],
             "npred" = npred, "Apred" = Apred, "predX" = predX)

#PHASE 1: Fit fixed effects
#Parameters to estimate
params <- list("beta0" = 0, "beta1" = 1, "delta1" = 1,
               "alpha0" = logit(0.01), "alpha1" = 0, 
               "log_kappa" = log(sqrt(8)/range), "log_tau_O" = 1,
               "omega" = rep(0, mesh$n))

#Map
map <- list(
  "log_kappa" = as.factor(NA), 
  "log_tau_O" = as.factor(NA), 
  "omega" = factor(rep(NA, mesh$n)))

#Make AD objective function
obj3 <- MakeADFun(data = data, parameters = params, random = random, map = map, DLL="Simulation_PO")

#Trace parameters
#obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)

#Calculate standard deviations
out3 <- sdreport(obj3)

params <- list("beta0" = opt1$par[1], "beta1" = opt1$par[2], "delta1" = opt1$par[3],
               "alpha0" = opt1$par[4], "alpha1" = opt1$par[5],
               "log_kappa" = log(sqrt(8)/range), "log_tau_O" = 1,
               "omega" = rep(0, mesh$n))

#Map
map <- list(
  "beta0" = as.factor(NA),
  "beta1" = as.factor(NA),
  "delta1" = as.factor(NA),
  "alpha0" = as.factor(NA),
  "alpha1" = as.factor(NA))

#Define random effects
random = c("omega")

#Make AD objective function
obj4 <- MakeADFun(data = data, parameters = params, map = map, random = random, DLL="Simulation_PO")

#Trace parameters
#obj4$env$tracepar <- TRUE

#Minimize objective function
opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr)

#Calculate standard deviations
out4 <- sdreport(obj4)

#Store output

output <- merge(output,  
                as.data.frame(round(rbind(summary(out3)[!grepl("omega|sigma_O|\\bkappa\\b|\\btau_O\\b", rownames(summary(out3))),],
                                          summary(out4)[grepl("\\bkappa\\b|\\btau_O\\b", rownames(summary(out4))),]), digits = 2)),
                by='row.names', all=TRUE)

colnames(output)[5:6] <- c("PO Estimate", "PO Std. Error")
output$Row.names[7:8] <- c("N[1]", "N[2]")
rownames(output) <- output$Row.names
output <- output[,-1]

convergence <- data.frame("pdHess" = c(out1$pdHess, out2$pdHess, out3$pdHess, out4$pdHess),
                          "grd.size" = c(all(abs(out1$gradient.fixed) < 0.01),
                                         all(abs(out2$gradient.fixed) < 0.01),
                                         all(abs(out3$gradient.fixed) < 0.01),
                                         all(abs(out4$gradient.fixed) < 0.01)))
rownames(convergence) <- c("IM_fixed", "IM_random", "PO_fixed", "PO_random")

out <- list(output, convergence)

message("finished")

#-----------#
#-Save file-#
#-----------#

ID <- length(list.files("Z:/Monarch/DataAnalysis/Simulation/Output/")) + 1
save(out, file = paste("Z:/Monarch/DataAnalysis/Simulation/Output/output", ID, ".Rds", sep=""))

message("saved")
print(iter)
