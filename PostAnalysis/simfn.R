simfn <- function(){

#-----------#
#-Libraries-#
#-----------#

library(spatstat)
library(RandomFields)
library(INLA)
library(maptools)
library(sf)

#-----------#
#-Functions-#
#-----------#

#Function to create dual mesh (original code from xxx)
source("Z:/Monarch/DataAnalysis/Simulation/rDualMesh.R")

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

#Period 1
# win1 <- owin(c(1, 2.75), c(0.25, 2))

#Period 2
# win2 <- owin(c(0, 3), c(0, 3), poly = list(x=c(0.25,0.75,0.75,2.75,2.75,0.25), y=c(0.25,0.25,2.25,2.25,2.75,2.75)))

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
beta <- c(6, 0.5)

#Change in intensity
delta <- 1

#Parameters of thinning function
alpha <- c(logit(0.01), 0.85)

#Ecological covariate
set.seed(100)
RandomFields::RFoptions(spConform=FALSE)
xcov1 <- RFsimulate(model = RMgauss(scale = 1.25), x = x0, y = y0, grid = TRUE)
xcov2 <- xcov1 + RFsimulate(model = RMgauss(var = 0.5, scale = 9.5), x = x0, y = y0, grid = TRUE)
x1 <- (xcov1 - mean(c(xcov1, xcov2))/sd(c(xcov1, xcov2)))
x2 <- (xcov2 - mean(c(xcov1, xcov2))/sd(c(xcov1, xcov2)))

#Thinning covariate
pcov <- RFsimulate(model = RMgauss(scale = 1), x = x0, y = y0, grid = TRUE)
pcov <- (pcov - min(pcov))/sd(pcov)
set.seed(NULL)

#Variance of Matern covaraince
sigma2x <- 0.2
#Range of Matern covaraince
range <- 1.2

#Smoothness of Matern covaraince (neighborhood size?)
nu <- 1

#Simulate spatial covariance
zcov <- RFsimulate(model = RMmatern(var = sigma2x, scale = range, nu = nu), x = x0, y = y0, grid = TRUE)

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

#Distances
dist <- NULL

for(j in 1:ncount){
  unit <- disc(radius = 0.15, centre = as.numeric(u.loc[j,]))
  count[j] <- PP1[unit]$n
  count[j] <- rbinom(1, count[j], 0.9)
  Cov[j] <- mean(as.im(x1, W=Mwin)[unit])
  unit.A[j] <- area(unit)
}

#Count dataframe
countdf <- data.frame(count, Cov, unit.A, u.loc)

#Compile data
npixels <- prod(lambda1$dim)/area(win)
intensity1 <- as.SpatialGridDataFrame.im(lambda1[win]/npixels)
intensity2 <- as.SpatialGridDataFrame.im(lambda2[win]/npixels)
domain <- st_as_sfc(win)
# domain1 <- st_as_sfc(win1)
# domain2 <- st_as_sfc(win2)
# pop1 <- st_as_sf(PP1[win1])
# pop2 <- st_as_sf(PP2[win2])
pop1 <- st_as_sf(PP1[win])
pop2 <- st_as_sf(PP2[win])
# po1 <- st_as_sf(PO1[win1])
# po2 <- st_as_sf(PO2[win2])
po1 <- st_as_sf(PO1[win])
po2 <- st_as_sf(PO2[win])
sites <- st_as_sfc(discs(as.ppp(u.loc, win), radii = 0.15))

sim.data <- list(intensity1 = intensity1, intensity2 = intensity2,
                 domain = domain, # domain1 = domain1, domain2 = domain2,
                 pop1 = pop1, pop2 = pop2, po1 = po1, po2 = po2, sites = sites,
                 countdf = countdf)

return(sim.data)

}
