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
library(gridExtra)

#-----------#
#-Functions-#
#-----------#

#Function to create dual mesh (original code from xxx)
source("~/Monarchs/DataAnalysis/rDualMesh.R")

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

#Mesh window
Mwin <- owin(c(min(x0), max(x0)), c(min(y0), max(y0)))

#---------------#
#-Simulate LGCP-#
#---------------#

#Parameters of intensity function
beta <- c(4.5, 0.5)

#Change in intensity
delta <- c(0.5, -0.75)

#Parameters of thinning function
alpha <- c(logit(0.005), 0.85)

#Ecological covariate
RandomFields::RFoptions(spConform=FALSE)
xcov1 <- RFsimulate(model = RMgauss(scale = 1.25), x = x0, y = y0, grid = TRUE)
xcov2 <- xcov1 + RFsimulate(model = RMgauss(var = 0.5, scale = 9.5), x = x0, y = y0, grid = TRUE)
xcov3 <- xcov2 + RFsimulate(model = RMgauss(var = 0.5, scale = 9.5), x = x0, y = y0, grid = TRUE)
x1 <- (xcov1 - mean(c(xcov1, xcov2, xcov3))/sd(c(xcov1, xcov2, xcov3)))
x2 <- (xcov2 - mean(c(xcov1, xcov2, xcov3))/sd(c(xcov1, xcov2, xcov3)))
x3 <- (xcov3 - mean(c(xcov1, xcov2, xcov3))/sd(c(xcov1, xcov2, xcov3)))

#Thinning covariate
pcov <- RFsimulate(model = RMgauss(scale = 1), x = x0, y = y0, grid = TRUE)
pcov <- (pcov - min(pcov))/sd(pcov)

#Variance of Matern covaraince
sigma2x <- 0.2

#Range of Matern covaraince
range <- 1.2

#Smoothness of Matern covaraince (neighborhood size?)
nu <- 1

#Simulate spatial covariance
zcov <- RFsimulate(model = RMmatern(var = sigma2x, scale = sqrt(8)/range, nu = nu), x = x0, y = y0, grid = TRUE)

#Simulate spatio-temporal covariance
ztcov1 <- RFsimulate(model = RMgauss(var = 0.1, scale = 1), x = x0, y = y0, grid = TRUE)
ztcov2 <- RFsimulate(model = RMgauss(var = 0.1, scale = 1), x = x0, y = y0, grid = TRUE)
ztcov3 <- RFsimulate(model = RMgauss(var = 0.1, scale = 1), x = x0, y = y0, grid = TRUE)

#Intensity function
X1 <- as.im(x1, W = Mwin)[win]
X2 <- as.im(x2, W = Mwin)[win]
X3 <- as.im(x3, W = Mwin)[win]
z <- as.im(zcov, W = Mwin)[win]
zt1 <- as.im(ztcov1, W = Mwin)[win]
zt2 <- as.im(ztcov2, W = Mwin)[win]
zt3 <- as.im(ztcov3, W = Mwin)[win]
lambda1 <- exp(beta[1] + beta[2] * X1 + z + zt1)
lambda2 <- exp(beta[1] + delta[1] + beta[2] * X2 + z + zt2)
lambda3 <- exp(beta[1] + delta[1] + delta[2] + beta[2] * X3 + z + zt3)

#Format intensity function
lambda1 <- as.im(lambda1, W=win)
lambda2 <- as.im(lambda2, W=win)
lambda3 <- as.im(lambda3, W=win)

#Simulate LGCP
PP1 <- rpoispp(lambda1)[win]
PP2 <- rpoispp(lambda2)[win]
PP3 <- rpoispp(lambda3)[win]

# #Lambda full
# xF <- as.im(x1, W = Mwin)
# zF <- as.im(zcov, W = Mwin)
# lambdaF <- exp(beta[1] + beta[2] * xF + zF)
# lambdaF <- as.im(lambdaF, W=Mwin)
# 
# #Plot
# par(mar=c(0, 0, 0, 0))
# plot(lambdaF, main = "")
# plot(mesh, add = T)
# points(PP1, pch = 19, col = "white")
# points(PP2, pch = 19, col = "black")
# points(PP3, pch = 19, col = "red")

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

Z2 <- rbinom(n = PP2$n, size = 1, prob = thin[PP2])
PO2 <- PP2[Z2==1]

Z3 <- rbinom(n = PP3$n, size = 1, prob = thin[PP3])
PO3 <- PP3[Z3==1]

#points(PO, pch = 19, col = "red")

par(mar = c(0,1,1,1))
layout(matrix(c(3,6,7,2,5,9,1,4,8), nrow = 3, ncol = 3, byrow = TRUE))

plot(lambda3, main = paste("Pop3 =", PP3$n, "po =", PO3$n))
points(PP3, pch = 19, col = "black")
points(PO3, pch = 19, col = "red")
plot(lambda2, main = paste("Pop2 =", PP2$n, "po =", PO2$n))
points(PP2, pch = 19, col = "black")
points(PO2, pch = 19, col = "red")
plot(lambda1, main = "")
points(PP1, pch = 19, col = "black")
points(PO1, pch = 19, col = "red")

#-----------------#
#-Simulate counts-#
#-----------------#

#Number of counts
ncount <- 25

#Count locations
u.loc <- expand.grid(seq(0.5, 2.5, 0.5), seq(0.5, 2.5, 0.5))

#Counts
count <- NULL

#Covariate value @ location
Cov <- NULL

#Unit area
unit.A <- NULL

for(j in 1:ncount){
  unit <- disc(radius = 0.15, centre = as.numeric(u.loc[j,]))
  count[j] <- PP1[unit]$n
  count[j] <- rbinom(1, count[j], 0.9)
  Cov[j] <- mean(as.im(x1, W=Mwin)[unit])
  unit.A[j] <- area(unit)
  plot(unit, border = "white", add = T)
  text(x = u.loc[j,1], y = u.loc[j,2], label = paste(count[j]), col = "white", font = 2, cex=1.5)
}

#Count dataframe
countdf <- data.frame(count, Cov, unit.A, u.loc)

title(main = paste("Pop1 =", PP1$n, "po =", PO1$n, "count =", sum(countdf$count)))

#countdf <- countdf[-c(2,4,6:10,12,14,16:20,21,24),]

#---------------#
#-SPDE approach-#
#---------------#

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
nobs <- sum(PO1$n, PO2$n, PO3$n)
ncount <- dim(countdf)[1]

#Period index
t_i <- rep(0:2, c(PO1$n, PO2$n, PO3$n))
t_n <- 3

#Counts
counts <- countdf$count

#Area @ each count location
Area <- countdf$unit.A

#Covaraites @ nodes
nloc <- as.ppp(mesh$loc[,1:2], W = Mwin)
nCov <- matrix(NA, nrow = nodes, ncol = t_n)
nCov[,1] <- as.im(x1, W = Mwin)[nloc]
nCov[,2] <- as.im(x2, W = Mwin)[nloc]
nCov[,3] <- as.im(x3, W = Mwin)[nloc]
nBias <- as.im(pcov, W = Mwin)[nloc]

#Covariates @ obs
Cov <- c(as.im(x1, W = Mwin)[PO1],
         as.im(x2, W = Mwin)[PO2],
         as.im(x3, W = Mwin)[PO3],
         countdf$Cov)

Bias <- c(as.im(pcov, W = Mwin)[PO1],
          as.im(pcov, W = Mwin)[PO2],
          as.im(pcov, W = Mwin)[PO3])

#Location of obsevations
locxy <- rbind(cbind(PO1$x, PO1$y)[,2:1],
               cbind(PO2$x, PO2$y)[,2:1],
               cbind(PO3$x, PO3$y)[,2:1],
               as.matrix(countdf[,4:5]))

#Projection matrix of observations
A <- as(inla.spde.make.A(mesh, locxy), "dgTMatrix")

#------------------#
#-Integrated Model-#
#------------------#

#Compile TMB code (only once)
#compile("IM_time.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/IM_time"))

#Compile data
data <- list("nodes" = nodes, "nobs" = nobs, "ncount" = ncount, "t_i" = t_i, "t_n" = t_n,
             "weight" = weight, "area" = Area, "A" = A, "counts" = counts,
             "nBias" = nBias, "nCov" = nCov, "Bias" = Bias, "Cov" = Cov,
             "spde" = spde$param.inla[c("M0","M1","M2")])

#Parameters to estimate
params <- list("beta0" = 0, "beta1" = 1, 
               "delta1" = 0, "delta2" = 0,
               "alpha0" = 0, "alpha1" = 0, 
               "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
               "omega" = rep(0, mesh$n),
               "epsilon" = matrix(0, nrow = nodes, ncol = t_n))

#PHASE 1: Fit fixed parameters

map.eps <- factor(rep(NA, nodes*t_n))

map <- list(
            #"beta0" = as.factor(NA),
            #"beta1" = as.factor(NA),
            #"delta1" = as.factor(NA),
            #"delta2" = as.factor(NA),
            #"alpha0" = as.factor(NA), 
            #"alpha1" = as.factor(NA),
            "log_kappa" = as.factor(NA), 
            "log_tau_O" = as.factor(NA), 
            "log_tau_E" = as.factor(NA),
            "omega" = factor(rep(NA, mesh$n)),
            "epsilon" = map.eps)


#Make AD objective function
obj1 <- MakeADFun(data = data, parameters = params, map = map, DLL="IM_time")

#Trace parameters
obj1$env$tracepar <- TRUE

#Minimize objective function
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
out1 <- sdreport(obj1)

#PHASE 2: Fit random effects
params <- list("beta0" = opt1$par[1], "beta1" = opt1$par[2],
               "delta1" = opt1$par[3], "delta2" = opt1$par[4],
               "alpha0" = opt1$par[5], "alpha1" = opt1$par[6], 
               "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
               "omega" = rep(0, mesh$n),
               "epsilon" = matrix(0, nrow = nodes, ncol = t_n))

map <- list(
  "beta0" = as.factor(NA),
  "beta1" = as.factor(NA),
  "delta1" = as.factor(NA),
  "delta2" = as.factor(NA),
  "alpha0" = as.factor(NA), 
  "alpha1" = as.factor(NA))

#Define random effects
random = c("omega", "epsilon")

#Make AD objective function
obj2 <- MakeADFun(data = data, parameters = params, map = map, random = random, DLL="IM_time")

#Trace parameters
obj2$env$tracepar <- TRUE

#Minimize objective function
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
out2 <- sdreport(obj2)

#------------#
#-Prediction-#
#------------#

#Beta0
beta0 <- opt1$par[1]
beta0.SE <- summary(out1, "fixed")[1, "Std. Error"]
beta0.low <- beta0 - beta0.SE
beta0.high <- beta0 + beta0.SE

#Beta1
beta1 <- opt1$par[2]
beta1.SE <- summary(out1, "fixed")[2, "Std. Error"]
beta1.low <- beta1 - beta1.SE
beta1.high <- beta1 + beta1.SE

#Delta1
delta1 <- opt1$par[3]
delta1.SE <- summary(out1, "fixed")[3, "Std. Error"]
delta1.low <- delta1 - delta1.SE
delta1.high <- delta1 + delta1.SE

#Delta2
delta2 <- opt1$par[4]
delta2.SE <- summary(out1, "fixed")[4, "Std. Error"]
delta2.low <- delta2 - delta2.SE
delta2.high <- delta2 + delta2.SE

#Alpha0
alpha0 <- opt1$par[5]
alpha0.SE <- summary(out1, "fixed")[5, "Std. Error"]

#Alpha1
alpha1 <- opt1$par[6]
alpha1.SE <- summary(out1, "fixed")[6, "Std. Error"]

#TauO
log.tauO <- out2$par.fixed[2]
tauO <- exp(log.tauO)
log.tauO.SE <- summary(out2, "fixed")[2, "Std. Error"]
tauO.low <- exp(log.tauO - log.tauO.SE)
tauO.high <- exp(log.tauO + log.tauO.SE)

#Omega
omg <- out2$par.random[1:nodes]
omg.SE <- summary(out2, "random")[1:nodes, "Std. Error"]
omg.low <- omg - omg.SE
omg.high <- omg + omg.SE
omg <- omg/tauO
omg.low <- omg.low/tauO.low
omg.high <- omg.high/tauO.high

#TauE
log.tauE <- out2$par.fixed[3]
tauE <- exp(log.tauE)
log.tauE.SE <- summary(out2, "fixed")[3, "Std. Error"]
tauE.low <- exp(log.tauE - log.tauE.SE)
tauE.high <- exp(log.tauE + log.tauE.SE)

#Epsilon
eps <- out2$par.random[(nodes+1):(nodes + nodes*t_n)]
dim(eps) <- c(nodes, t_n)
eps.SE <- summary(out2, "random")[(nodes+1):(nodes + nodes*t_n), "Std. Error"]
dim(eps.SE) <- c(nodes, t_n)
eps.low <- eps - eps.SE
eps.high <- eps + eps.SE
eps <- eps/tauE
eps.low <- eps.low/tauE.low
eps.high <- eps.high/tauE.high

#Predicted values
Grid <- as.matrix(expand.grid(seq(0.01,2.99,0.01), seq(0.01,2.99,0.01)))
GridA <- inla.spde.make.A(mesh, loc = Grid)

x1.cov <- interp.im(X1, x = Grid[,1], y = Grid[,2])
x2.cov <- interp.im(X2, x = Grid[,1], y = Grid[,2])
x3.cov <- interp.im(X3, x = Grid[,1], y = Grid[,2])


Pred <- cbind(exp(beta0 + beta1 * x1.cov + GridA %*% omg + GridA %*% eps[,1]), 
              exp(beta0 + beta1 * x2.cov + delta1 + GridA %*% omg + GridA %*% eps[,2]),
              exp(beta0 + beta1 * x3.cov + delta1 + delta2 + GridA %*% omg + GridA %*% eps[,3]),
              exp(beta0.low + beta1.low * x1.cov + GridA %*% omg.low + GridA %*% eps.low[,1]), 
              exp(beta0.low + beta1.low * x2.cov + delta1.low + GridA %*% omg.low + GridA %*% eps.low[,2]),
              exp(beta0.low + beta1.low * x3.cov + delta1.low + delta2.low + GridA %*% omg.low + GridA %*% eps.low[,3]),
              exp(beta0.high + beta1.high * x1.cov + GridA %*% omg.high + GridA %*% eps.high[,1]), 
              exp(beta0.high + beta1.high * x2.cov + delta1.high + GridA %*% omg.high + GridA %*% eps.high[,2]),
              exp(beta0.high + beta1.high * x3.cov + delta1.high + delta2.high + GridA %*% omg.high + GridA %*% eps.high[,3]))

Pred.df <- data.frame(T1 = as.vector(Pred[,1]), T2 = as.vector(Pred[,2]), T3 = as.vector(Pred[,3]),
                      T1.low = as.vector(Pred[,4]), T2.low = as.vector(Pred[,5]), T3.low = as.vector(Pred[,6]),
                      T1.high = as.vector(Pred[,7]), T2.high = as.vector(Pred[,8]), T3.high = as.vector(Pred[,9]),
                      X = Grid[,1]*100, Y = Grid[,2]*100)

Pred.T1 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1"), W = win)
Pred.T2 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2"), W = win)
Pred.T3 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3"), W = win)
Pred.T1.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1.low"), W = win)
Pred.T2.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2.low"), W = win)
Pred.T3.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3.low"), W = win)
Pred.T1.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1.high"), W = win)
Pred.T2.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2.high"), W = win)
Pred.T3.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3.high"), W = win)

plot(Pred.T3, main = "")
plot(Pred.T2, main = "")
plot(Pred.T1, main = "")

output1 <- round(cbind(c(beta, delta, alpha), 
                c(beta0, beta1, delta1, delta2, alpha0, alpha1), 
                c(beta0.SE, beta1.SE, delta1.SE, delta2.SE, alpha0.SE, alpha1.SE)),
                 digits = 2)

colnames(output1) <- c("Truth", "Est. Mean", "Est. SE")

output2 <- round(cbind(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9),
                       c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9)), digits = 2)
rownames(output2) <- c("Pop1", "Pop2", "Pop3")
colnames(output2) <- c("Truth", "Estimated")

frame()
vps <- gridBase::baseViewports()
grid::pushViewport(vps$inner, vps$figure, vps$plot)  
grid::grid.draw(tableGrob(output1))
grid::popViewport(3)

frame()
vps <- gridBase::baseViewports()
grid::pushViewport(vps$inner, vps$figure, vps$plot)  
grid::grid.draw(tableGrob(output2))
grid::popViewport(3)

par(mar = c(0,3,3,1))
plot(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9),
     type = "l", xaxt = "n", lwd = 1.75, col = "blue",
     xlab = "Time period", ylab = "Abundance",
     ylim = c(0, max(c(mean(Pred.T1.high),
                       mean(Pred.T2.high),
                       mean(Pred.T3.high)))*9*1.1))
axis(1, at = seq(1,3,1))
points(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9),
       type = "p", col = "white", pch = 19, cex = 2.5)
points(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9),
       type = "l", lwd = 1.75, col = "black")
points(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9),
       type = "p", col = "white", pch = 19, cex = 2.5)
points(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9),
       type = "p", col = "blue", pch = 19)
arrows(seq(1,3,1), c(mean(Pred.T1.low)*9, mean(Pred.T2.low)*9, mean(Pred.T3.low)*9), 
       seq(1,3,1), c(mean(Pred.T1.high)*9, mean(Pred.T2.high)*9, mean(Pred.T3.high)*9), 
       length=0.1, angle=90, code=3, col = "blue", lwd = 1.75)
points(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9),
       type = "p", col = "black", pch = 19)

#---------------------#
#-Presence-only model-#
#---------------------#

#Compile TMB code (only once)
#compile("PO_time.cpp")

#Load TMB code
dyn.load(dynlib("./DataAnalysis/PO_time"))

#Compile data
data <- list("nodes" = nodes, "nobs" = nobs, "t_i" = t_i, "t_n" = t_n,
             "weight" = weight, "A" = A[1:nobs,],
             "nBias" = nBias, "nCov" = nCov, "Bias" = Bias, "Cov" = Cov[1:nobs],
             "spde" = spde$param.inla[c("M0","M1","M2")])

#Parameters to estimate
params <- list("beta0" = 0, "beta1" = 1, 
               "delta1" = 0, "delta2" = 0,
               "alpha0" = 0, "alpha1" = 0, 
               "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
               "omega" = rep(0, mesh$n),
               "epsilon" = matrix(0, nrow = nodes, ncol = t_n))

#PHASE 1: Fit fixed parameters

map.eps <- factor(rep(NA, nodes*t_n))

map <- list(
  #"beta0" = as.factor(NA),
  #"beta1" = as.factor(NA),
  #"delta1" = as.factor(NA),
  #"delta2" = as.factor(NA),
  #"alpha0" = as.factor(NA), 
  #"alpha1" = as.factor(NA),
  "log_kappa" = as.factor(NA), 
  "log_tau_O" = as.factor(NA), 
  "log_tau_E" = as.factor(NA),
  "omega" = factor(rep(NA, mesh$n)),
  "epsilon" = map.eps)


#Make AD objective function
obj3 <- MakeADFun(data = data, parameters = params, map = map, DLL="PO_time")

#Trace parameters
obj3$env$tracepar <- TRUE

#Minimize objective function
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
out3 <- sdreport(obj3)

#PHASE 2: Fit random effects
params <- list("beta0" = opt3$par[1], "beta1" = opt3$par[2],
               "delta1" = opt3$par[3], "delta2" = opt3$par[4],
               "alpha0" = opt3$par[5], "alpha1" = opt3$par[6], 
               "log_kappa" = 0, "log_tau_O" = 1, "log_tau_E" = 1,
               "omega" = rep(0, mesh$n),
               "epsilon" = matrix(0, nrow = nodes, ncol = t_n))

map <- list(
  "beta0" = as.factor(NA),
  "beta1" = as.factor(NA),
  "delta1" = as.factor(NA),
  "delta2" = as.factor(NA),
  "alpha0" = as.factor(NA), 
  "alpha1" = as.factor(NA))

#Define random effects
random = c("omega", "epsilon")

#Make AD objective function
obj4 <- MakeADFun(data = data, parameters = params, map = map, random = random, DLL="PO_time")

#Trace parameters
obj4$env$tracepar <- TRUE

#Minimize objective function
opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr)
out4 <- sdreport(obj4)

#------------#
#-Prediction-#
#------------#

#Beta0
beta0 <- opt3$par[1]
beta0.SE <- summary(out3, "fixed")[1, "Std. Error"]
beta0.low <- beta0 - beta0.SE
beta0.high <- beta0 + beta0.SE

#Beta1
beta1 <- opt3$par[2]
beta1.SE <- summary(out3, "fixed")[2, "Std. Error"]
beta1.low <- beta1 - beta1.SE
beta1.high <- beta1 + beta1.SE

#Delta1
delta1 <- opt3$par[3]
delta1.SE <- summary(out3, "fixed")[3, "Std. Error"]
delta1.low <- delta1 - delta1.SE
delta1.high <- delta1 + delta1.SE

#Delta2
delta2 <- opt3$par[4]
delta2.SE <- summary(out3, "fixed")[4, "Std. Error"]
delta2.low <- delta2 - delta2.SE
delta2.high <- delta2 + delta2.SE

#Alpha0
alpha0 <- opt3$par[5]
alpha0.SE <- summary(out3, "fixed")[5, "Std. Error"]

#Alpha1
alpha1 <- opt3$par[6]
alpha1.SE <- summary(out3, "fixed")[6, "Std. Error"]

#TauO
log.tauO <- out4$par.fixed[2]
tauO <- exp(log.tauO)
log.tauO.SE <- summary(out4, "fixed")[2, "Std. Error"]
tauO.low <- exp(log.tauO - log.tauO.SE)
tauO.high <- exp(log.tauO + log.tauO.SE)

#Omega
omg <- out4$par.random[1:nodes]
omg.SE <- summary(out4, "random")[1:nodes, "Std. Error"]
omg.low <- omg - omg.SE
omg.high <- omg + omg.SE
omg <- omg/tauO
omg.low <- omg.low/tauO.low
omg.high <- omg.high/tauO.high

#TauE
log.tauE <- out4$par.fixed[3]
tauE <- exp(log.tauE)
log.tauE.SE <- summary(out4, "fixed")[3, "Std. Error"]
tauE.low <- exp(log.tauE - log.tauE.SE)
tauE.high <- exp(log.tauE + log.tauE.SE)

#Epsilon
eps <- out4$par.random[(nodes+1):(nodes + nodes*t_n)]
dim(eps) <- c(nodes, t_n)
eps.SE <- summary(out4, "random")[(nodes+1):(nodes + nodes*t_n), "Std. Error"]
dim(eps.SE) <- c(nodes, t_n)
eps.low <- eps - eps.SE
eps.high <- eps + eps.SE
eps <- eps/tauE
eps.low <- eps.low/tauE.low
eps.high <- eps.high/tauE.high

#Predicted values
Grid <- as.matrix(expand.grid(seq(0.01,2.99,0.01), seq(0.01,2.99,0.01)))
GridA <- inla.spde.make.A(mesh, loc = Grid)

x1.cov <- interp.im(X1, x = Grid[,1], y = Grid[,2])
x2.cov <- interp.im(X2, x = Grid[,1], y = Grid[,2])
x3.cov <- interp.im(X3, x = Grid[,1], y = Grid[,2])


Pred <- cbind(exp(beta0 + beta1 * x1.cov + GridA %*% omg + GridA %*% eps[,1]), 
              exp(beta0 + beta1 * x2.cov + delta1 + GridA %*% omg + GridA %*% eps[,2]),
              exp(beta0 + beta1 * x3.cov + delta1 + delta2 + GridA %*% omg + GridA %*% eps[,3]),
              exp(beta0.low + beta1.low * x1.cov + GridA %*% omg.low + GridA %*% eps.low[,1]), 
              exp(beta0.low + beta1.low * x2.cov + delta1.low + GridA %*% omg.low + GridA %*% eps.low[,2]),
              exp(beta0.low + beta1.low * x3.cov + delta1.low + delta2.low + GridA %*% omg.low + GridA %*% eps.low[,3]),
              exp(beta0.high + beta1.high * x1.cov + GridA %*% omg.high + GridA %*% eps.high[,1]), 
              exp(beta0.high + beta1.high * x2.cov + delta1.high + GridA %*% omg.high + GridA %*% eps.high[,2]),
              exp(beta0.high + beta1.high * x3.cov + delta1.high + delta2.high + GridA %*% omg.high + GridA %*% eps.high[,3]))

Pred.df <- data.frame(T1 = as.vector(Pred[,1]), T2 = as.vector(Pred[,2]), T3 = as.vector(Pred[,3]),
                      T1.low = as.vector(Pred[,4]), T2.low = as.vector(Pred[,5]), T3.low = as.vector(Pred[,6]),
                      T1.high = as.vector(Pred[,7]), T2.high = as.vector(Pred[,8]), T3.high = as.vector(Pred[,9]),
                      X = Grid[,1]*100, Y = Grid[,2]*100)

Pred.T1 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1"), W = win)
Pred.T2 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2"), W = win)
Pred.T3 <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3"), W = win)
Pred.T1.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1.low"), W = win)
Pred.T2.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2.low"), W = win)
Pred.T3.low <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3.low"), W = win)
Pred.T1.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T1.high"), W = win)
Pred.T2.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T2.high"), W = win)
Pred.T3.high <- as.im(reshape2::acast(Pred.df, X ~ Y, value.var = "T3.high"), W = win)

par(mar = c(0,1,1,1))
layout(matrix(c(3,6,7,2,5,9,1,4,8), nrow = 3, ncol = 3, byrow = TRUE))

plot(lambda3, main = paste("Pop3 =", PP3$n, "po =", PO3$n))
points(PP3, pch = 19, col = "black")
points(PO3, pch = 19, col = "red")
plot(lambda2, main = paste("Pop2 =", PP2$n, "po =", PO2$n))
points(PP2, pch = 19, col = "black")
points(PO2, pch = 19, col = "red")
plot(lambda1, main = paste("Pop1 =", PP1$n, "po =", PO1$n))
points(PP1, pch = 19, col = "black")
points(PO1, pch = 19, col = "red")

plot(Pred.T3, main = "")
plot(Pred.T2, main = "")
plot(Pred.T1, main = "")

output1 <- round(cbind(c(beta, delta, alpha), 
                       c(beta0, beta1, delta1, delta2, alpha0, alpha1), 
                       c(beta0.SE, beta1.SE, delta1.SE, delta2.SE, alpha0.SE, alpha1.SE)),
                 digits = 2)

colnames(output1) <- c("Truth", "Est. Mean", "Est. SE")

output2 <- round(cbind(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9),
                       c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9)), digits = 2)
rownames(output2) <- c("Pop1", "Pop2", "Pop3")
colnames(output2) <- c("Truth", "Estimated")

frame()
vps <- gridBase::baseViewports()
grid::pushViewport(vps$inner, vps$figure, vps$plot)  
grid::grid.draw(tableGrob(output1))
grid::popViewport(3)

frame()
vps <- gridBase::baseViewports()
grid::pushViewport(vps$inner, vps$figure, vps$plot)  
grid::grid.draw(tableGrob(output2))
grid::popViewport(3)

par(mar = c(0,3,3,1))
plot(log(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9)),
     type = "l", xaxt = "n", lwd = 1.75, col = "blue",
     xlab = "Time period", ylab = "Abundance",
     ylim = c(0, log(max(c(mean(Pred.T1),
                       mean(Pred.T2),
                       mean(Pred.T3)))*9)*1.1))
axis(1, at = seq(1,3,1))
points(log(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9)),
       type = "p", col = "white", pch = 19, cex = 2.5)
points(log(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9)),
       type = "l", lwd = 1.75, col = "black")
points(log(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9)),
       type = "p", col = "white", pch = 19, cex = 2.5)
points(log(c(mean(Pred.T1)*9, mean(Pred.T2)*9, mean(Pred.T3)*9)),
       type = "p", col = "blue", pch = 19)
points(log(c(mean(lambda1)*9, mean(lambda2)*9, mean(lambda3)*9)),
       type = "p", col = "black", pch = 19)
