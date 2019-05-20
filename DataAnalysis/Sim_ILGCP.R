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
library(RandomFields)
library(INLA)
library(rgeos)
library(TMB)

#-------------------#
#-Working directory-#
#-------------------#

setwd("~/Users/farrm/Documents/GitHub/Monarchs/DataAnalysis")

#-----------#
#-Functions-#
#-----------#

#Function to simulate LGCP (original code from XXX)
source("rthinLGCP.R")

#Function to create dual mesh (original code from xxx)
source("rDualMesh.R")

#Logit function
logit <- function(pp) 
{ 
  log(pp) - log(1-pp) 
}

#----------------#
#-Spatial Domain-#
#----------------#

#Spatial domain
win <- owin(c(0, 3), c(0, 3)) #Window
loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0)) #Coordinates

#Number of pixels
npix <- 300
spatstat.options(npixel = npix)

#Pixelize spatial domain
mask <- as.mask(win)

#Create triangulated mesh over domain
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1),
                     max.edge = c(0.3, 0.7), cutoff = 0.05)

#X range of mesh
x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = npix)

#Y range of mesh
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = npix)

#---------------#
#-Simulate LGCP-#
#---------------#

#Parameters of intensity function
beta <- c(4, 0.85)

#Parameters of thinning function
alpha_1 <- c(logit(0.5), 0.5)
alpha_2 <- c(logit(0.2), 0.5)

param <- c(beta, alpha_1, alpha_2)

#Ecological covariate
RandomFields::RFoptions(spConform=FALSE)
# xcov <- RFsimulate(model = RMgauss(scale = 0.8), x = x0, y = y0, grid = TRUE)
xcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))

#Thinning covariate
# pcov1 <- outer(x0, y0, function(x,y) sin(x) - y)
pcov1 <- RFsimulate(model = RMgauss(scale = 1.5), x = x0, y = y0, grid = TRUE)
# pcov2 <- outer(x0, y0, function(x,y) cos(y) - x)
pcov2 <- RFsimulate(model = RMgauss(scale = 2), x = x0, y = y0, grid = TRUE)

#Intensity function
mu <- beta[1] + beta[2] * xcov

#Format intensity function
mu <- as.im(mu, W=win)

#Thinning function
thin1 <- alpha_1[1] + alpha_1[2] * pcov1
thin2 <- alpha_2[1] + alpha_2[2] * pcov2

#Format thinning function
thin1 <- as.im(thin1, W=win)
thin2 <- as.im(thin2, W=win)

#Variance of Matern covaraince
sigma2x <- 0.2

#Range of Matern covaraince
range <- 1.2

#Smoothness of Matern covaraince (neighborhood size?)
nu <- 1

#Simulate LGCP
TPP1 <- rthinLGCP(model = 'matern', mu = mu, thin = thin1, var = sigma2x,
                 scale = range/sqrt(8), nu = nu, win = win, saveLambda = TRUE)

TPP2 <- rthinLGCP(model = 'matern', mu = mu, thin = thin2, var = sigma2x,
                  scale = range/sqrt(8), nu = nu, win = win, saveLambda = TRUE)

#---------------#
#-SPDE approach-#
#---------------#

#Create SPDE objects (Matern covariance)
spde <- inla.spde2.pcmatern(mesh, prior.range = c(0.05, 0.01), prior.sigma = c(1, 0.01))

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

#--------------#
#-Compile Data-#
#--------------#

#Node index
nodes <- NULL
nodes[1] <- mesh$n
nodes[2] <- npix*npix

#Observation index
nobs1 <- TPP1$n
nobs2 <- TPP2$n

#Location of obsevations
locxy1 <- cbind(TPP1$x, TPP1$y)[,2:1]
locxy2 <- cbind(TPP2$x, TPP2$y)[,2:1]

#Projection matrix of nodes
A1 <- as(as.matrix(Diagonal(nodes[1], rep(1, nodes[1]))), "dgTMatrix")

#Projection matrix of observations
A2 <- as(inla.spde.make.A(mesh, locxy1), "dgTMatrix")
A3 <- as(inla.spde.make.A(mesh, locxy2), "dgTMatrix")

#---------------#
#-Design matrix-#
#---------------#

#Convert covariate to image
xcov.im <- im(xcov, x0, y0)

#Design matrix for integration nodes
Xmat0 <- cbind(rep(1, nodes[1]), interp.im(xcov.im, x = mesh$loc[,1], y = mesh$loc[,2]))

#Design matrix for observation locations
Xmat1 <- cbind(rep(1, nobs1), interp.im(xcov.im, x = locxy1[,1], y = locxy1[,2]))
Xmat2 <- cbind(rep(1, nobs2), interp.im(xcov.im, x = locxy2[,1], y = locxy2[,2]))

#Convert covariate to image
pcov1.im <- im(pcov1, x0, y0)
pcov2.im <- im(pcov2, x0, y0)

#Design matrix for integration nodes
Pmat1_1 <- cbind(rep(1, nodes[1]), interp.im(pcov1.im, x = mesh$loc[,1], y = mesh$loc[,2]))

#Design matrix for observation locations
Pmat2_1 <- cbind(rep(1, nobs1), interp.im(pcov1.im, x = locxy1[,1], y = locxy1[,2]))

#Design matrix for integration nodes
Pmat1_2 <- cbind(rep(1, nodes[1]), interp.im(pcov2.im, x = mesh$loc[,1], y = mesh$loc[,2]))

#Design matrix for observation locations
Pmat2_2 <- cbind(rep(1, nobs2), interp.im(pcov2.im, x = locxy2[,1], y = locxy2[,2]))

#-----#
#-TMB-#
#-----#

#Compile TMB code (only once)
#compile("C:/Users/farrm/Documents/GitHub/Monarchs/DataAnalysis/ILGTCP.cpp")

#Load TMB code
dyn.load(dynlib("C:/Users/farrm/Documents/GitHub/Monarchs/DataAnalysis/ILGTCP"))

#Compile data
data <- list("nodes" = nodes[1], "nobs1" = nobs1, "nobs2" = nobs2, "weight" = weight,
             "A1" = A1, "A2" = A2, "A3" = A3, "Xmat0" = Xmat0, "Xmat1" = Xmat1, "Xmat2" = Xmat2,
             "Pmat1_1" = Pmat1_1, "Pmat2_1" = Pmat2_1, "Pmat1_2" = Pmat1_2, "Pmat2_2" = Pmat2_2,
             "spde" = spde$param.inla[c("M0","M1","M2")])

#Parameters to estimate
params1 <- list("beta" = c(0, 1), "alpha_1" = c(0, 0), "alpha_2" = c(0, 0), "log_kappa" = 0, "omega" = rep(0, mesh$n))

#PHASE 1: Fit fixed parameters
map1 <- list("log_kappa" = as.factor(NA), "omega" = factor(rep(NA, mesh$n)))

#Make AD objective function
obj1 <- MakeADFun(data = data, parameters = params1, map = map1, DLL="ILGTCP")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

#PHASE 2: Fit random effects
params2 <- list("beta" = opt1$par[1:2], "alpha_1" = opt1$par[3:4], "alpha_2" = opt1$par[5:6], "log_kappa" = 0, "omega" = rep(0, mesh$n))

map2 <- list("beta" = factor(rep(NA, 2)), "alpha_1" = factor(rep(NA, 2)), "alpha_2" = factor(rep(NA, 2)))

#Define random effects
random = "omega"

#Make AD objective function
obj2 <- MakeADFun(data = data, parameters = params2, map = map2, random = random, DLL="ILGTCP")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

#PHASE 3: Final run
params3 <- list("beta" = opt1$par[1:2], "alpha_1" = opt1$par[3:4], "alpha_2" = opt1$par[5:6], "log_kappa" = opt2$par[1], "omega" = rep(0, mesh$n))

#Make AD objective function
obj3 <- MakeADFun(data = data, parameters = params3, map = map1, random = random, DLL="ILGTCP")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)

#--------#
#-Output-#
#--------#

report <- obj3$report()

#------------#
#-Projection-#
#------------#

projection <- inla.mesh.projector(mesh = mesh, xlim = 0:3, ylim = 0:3, dims = c(301, 301))
proj.mean <- inla.mesh.project(projection, unlist(report[2]))
image(proj.mean)
