rthinLGCP <- local({
  
  rthinLGCP <- function(model="exp", mu = 0, thin = 1, param = NULL, ...,
                    win=NULL, saveLambda=TRUE, nsim=1, drop=TRUE) {
    ## validate
    if (!(is.numeric(mu) || is.function(mu) || is.im(mu))) 
      stop(paste(sQuote("mu"), "must be a constant, a function or an image"))
    if (is.numeric(mu) && !(length(mu) == 1)) 
      stop(paste(sQuote("mu"), "must be a single number"))
    if (!(is.numeric(thin) || is.function(thin) || is.im(thin))) 
      stop(paste(sQuote("thin"), "must be a constant, a function or an image"))
    if (is.numeric(thin) && !(length(thin) == 1)) 
      stop(paste(sQuote("thin"), "must be a single number"))
    ## check for outdated usage
    if(!all(nzchar(names(param))))
      stop("Outdated syntax of argument 'param' to rLGCP", call.=FALSE)
    ## 
    do.rLGCP(model=model, mu=mu, thin=thin, param=param, ...,
             win=win, saveLambda=saveLambda, nsim=nsim, drop=drop)
  }
  
  do.rLGCP <- function(model="exp", mu = 0, thin = 1, param = NULL, ...,
                       win=NULL, saveLambda=TRUE,
                       eps = NULL, dimyx = NULL, xy = NULL,
                       modelonly=FALSE, nsim=1, drop=TRUE) {
    ## make RF model object from RandomFields package
    ## get the 'model generator'
    modgen <- getRandomFieldsModelGen(model)
    ## now create a RandomFields 'model' object
    rfmodel <- do.call(modgen, append(as.list(param), list(...)))
    if(!inherits(rfmodel, "RMmodel"))
      stop("Unable to create RandomFields model object", call.=FALSE)
    
    ## secret exit
    if(modelonly)
      return(rfmodel)
    
    ## simulation window
    win.given <- !is.null(win)
    mu.image <- is.im(mu)
    win <- if(win.given) as.owin(win) else if(mu.image) as.owin(mu) else owin()
    
    if(win.given && mu.image && !is.subset.owin(win, as.owin(mu)))
      stop(paste("The spatial domain of the pixel image", sQuote("mu"),
                 "does not cover the simulation window", sQuote("win")))
    
    ## convert win to a mask
    w <- as.mask(w=win, eps=eps, dimyx=dimyx, xy=xy)
    xcol <- w$xcol
    yrow <- w$yrow
    dim <- w$dim
    xy <- expand.grid(x=xcol, y=yrow)
    xx <- xy$x
    yy <- xy$y
    
    muxy <- if(is.numeric(mu)) mu else
      if (is.function(mu)) mu(xx,yy) else
        lookup.im(mu, xx, yy, naok=TRUE, strict=TRUE)
    muxy[is.na(muxy)] <- -Inf
    
    logitPthin <- if(is.numeric(thin)) thin else
      if (is.function(thin)) thin(xx,yy) else
        lookup.im(thin, xx, yy, naok=TRUE, strict=TRUE)
    logitPthin[is.na(muxy)] <- -Inf
    
    stopifnot(nsim >= 1)
    result <- vector(mode="list", length=nsim)
    for(i in 1:nsim) {
      ## generate zero-mean Gaussian random field
      spc <- RandomFields::RFoptions()$general$spConform
      if(spc) RandomFields::RFoptions(spConform=FALSE)
      z <- RandomFields::RFsimulate(rfmodel, xcol, yrow, grid = TRUE)
      if(spc) RandomFields::RFoptions(spConform=TRUE)
      
      ## convert to log-Gaussian image
      logLambda <- muxy + z
      Lambda <- matrix(exp(logLambda), nrow=dim[1], ncol=dim[2], byrow=TRUE)
      Pthin <- matrix(1/(1+exp(-logitPthin)), nrow=dim[1], ncol=dim[2], byrow=TRUE)
      Lambdathin <- Lambda * Pthin
      Lambdathin <- as.im(Lambdathin, W=w)
      ## generate Poisson points
      X <- rpoispp(Lambdathin)[win]
      ## 
      if(saveLambda)
        attr(X, "Lambda") <- Lambdathin
      result[[i]] <- X
      ##Visualize
      par(mar=c(0, 0, 0, 0))
      plot(Lambdathin, main = "")
      points(X, pch = 19)
    }
    if(drop && nsim == 1)
      return(result[[1]])
    names(result) <- paste("Simulation", 1:nsim)
    return(as.solist(result))
  }
  
  rthinLGCP
})