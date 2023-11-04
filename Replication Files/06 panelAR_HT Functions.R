## Necessary functions for panelAR


prais.correct <- function(method,env.base){	
  ### pull objects from base environment
  lm.out <- get("lm.out",envir=env.base)
  yX <- get("yX",envir=env.base)
  units <- get("units",envir=env.base)
  panel.vec <- get("panel.vec",envir=env.base)
  e.mat <- get("e.mat",envir=env.base)
  panel.weight <- get("panel.weight",envir=env.base)
  N.times <- get("N.times",envir=env.base)
  N.units <- get("N.units",envir=env.base)
  bound.rho <- get("bound.rho",envir=env.base)
  rho.na.rm <- get("rho.na.rm",envir=env.base)
  rhotype <- get("rhotype",envir=env.base) 
  obs.mat <- get("obs.mat",envir=env.base)
  rank <- get("rank",envir=env.base)
  singular.ok <- get("singular.ok",envir=env.base)
  p <- ncol(yX)
  
  if(method=="none"){
    pw.output <- list(pw.lm = lm.out, pw.rho = NULL)
  } else{			
    ### calculate rho by panel
    rhos <- apply(e.mat, MARGIN=2, function(e) est.rho(e, k=rank,rhotype=rhotype))
    
    ### deal with missing rhos
    if(any(is.na(rhos))){
      if(!rho.na.rm){
        stop("Cannot estimate at least one panel-specific autocorrelation. Consider setting rho.na.rm to 'TRUE'.",call.=FALSE)
      } else {
        if(method=="ar1"){
          rhos <- na.omit(rhos)
        } else{
          rhos[is.na(rhos)] <- 0
          message("Setting panel-specific correlation to 0 for at least one panel because unable to estimate autocorrelation.")
        }
      }
    }
    
    ### set bounds
    if(bound.rho & any(abs(rhos)>1)){
      rhos[rhos > 1] <- 1
      rhos[rhos < -1] <- -1
      message("Panel-specific correlations bounded to [-1,1]")
    }
    
    ### create average rho if ar1
    if(method=="ar1"){
      # calculate total length of runs, not counting the first observation in each run
      weights <- apply(obs.mat,MARGIN=1,function(x) sum(rle(x)$length[rle(x)$values==TRUE]-1))
      # remove weights for rhos that were NA (and so that weight vectors matches rho vector)
      if(!is.null(attr(rhos,"na.action"))){
        weights <- weights[-attr(rhos,"na.action")]
      }
      
      if(panel.weight=="t"){
        # if panel.weight is "t", adjust by 1
        weights <- weights + 1
      }
      
      # calculate average rho
      rho.avg <- sum(rhos*weights)/sum(weights)
      rhos <- rep(rho.avg,N.units)
    }
    
    ### do Prais-Winsten transformation here
    transformed.data <- lapply(1:N.units,function(el) prais.transform(el,rhos=rhos,p=p,N.times=N.times,units=units,panel.vec=panel.vec,yX=yX,obs.mat=obs.mat))
    transformed.df <- do.call(rbind,transformed.data)
    transformed.X <- transformed.df[,2:p]
    transformed.y <- transformed.df[,1]
    
    # 2nd stage regression
    pw.lm <- lm(transformed.y ~ transformed.X - 1, singular.ok= singular.ok)
    
    if(method=="ar1"){
      rho.out <- rho.avg
    } else{
      rho.out <- rhos
    }
    
    pw.output <- list(pw.lm = pw.lm, pw.rho = rho.out)
  } # end of autocorrelation else statement
  
  return(pw.output)
}

panelAR_HT <- function(formula, data, panelVar, timeVar, autoCorr = c("ar1", 
                                                                      "none", "psar1"), panelCorrMethod = c("none","phet","pcse","pwls","parks"), rhotype ="breg", bound.rho = FALSE, rho.na.rm = FALSE, panel.weight = c("t-1", "t"), dof.correction = FALSE, complete.case = FALSE, seq.times = FALSE,  singular.ok=TRUE) {
  # save environment
  env.base <- environment()
  
  # save call
  call <- match.call()
  
  # check for presence of a formula and data object
  ind <- match(c("formula","data","timeVar"),names(call),nomatch=0)
  if(ind[1]==0){
    stop("A formula argument is required.",call.=FALSE)
  }
  if(ind[2]==0){
    stop("A data argument is required.",call.=FALSE)
  }
  if(ind[3]==0){
    stop("You must specify time ID using timeVar.",call.=FALSE)
  }
  
  # identify autocorrelation and panel structure
  autoCorr.method <- match.arg(autoCorr)
  panelCorr.method <- match.arg(panelCorrMethod)
  
  # identify 2nd stage method
  pMethod <- switch(panelCorr.method, "none"="OLS", "phet"="OLS", "pcse"="OLS", "pwls"="GLS", "parks"="GLS")
  
  # check rhotype
  rhotype <- match.arg(rhotype,c("breg","freg","dw","theil-nagar","scorr","theil"))
  
  # check panel.weight & rho.na.rm arguments
  panel.weight <- match.arg(panel.weight) 
  if(!is.logical(rho.na.rm)){
    stop("rho.na.rm should be logical argument.")
  }
  
  # make sure that timevar is in data object; check for NA values
  if (!(timeVar %in% colnames(data))) {
    stop("Please make sure that timeVar is in the data object.",call.=FALSE)
  }
  # extract time vector and check for NA values
  time.vec <- data[, timeVar]
  if (any(is.na(time.vec))) {
    stop("You cannot have NA values for the time variable.",call.=FALSE)
  }
  if (!is.integer(time.vec) & !is.numeric(time.vec)) {
    stop("The time variable must be defined as an integer.",call.=FALSE)
  }
  
  # extract panel variable and check if is.null
  if (!is.null(panelVar)){
    if (!(panelVar %in% colnames(data))) {
      stop("Please make sure that panelVar is in the data object.",call.=FALSE)
    } else{
      panel.vec <- as.character(as.vector(data[, panelVar]))
      if (any(is.na(panel.vec))) {
        stop("You cannot have NA values for the panel ID variable.",call.=FALSE)
      }
    }
  } else{   
    panel.vec <- rep("1",nrow(data))
    if(panelCorr.method!="none") warning("Without specifying a panel ID variable, data is assumed to come from a single panel and variance is assumed to be homoskedastic within that panel. Panel heteroskedasticity and/or correlation is ignored.",call.=FALSE)
  }
  
  # check for duplicate times in each panel
  if (!is.null(timeVar)) {
    if (any(by(time.vec, panel.vec, function(x) any(table(x) > 1)))) {
      stop("You must specify unique times for each observation in a given panel.",call.=FALSE)
    }
  }
  
  # sort dataframe
  order.index <- order(panel.vec, time.vec)
  data <- data[order.index, ]
  panel.vec <- panel.vec[order.index]
  time.vec <- time.vec[order.index]
  
  # run OLS and extract results
  lm.out <- lm(formula = formula, data = data,singular.ok=singular.ok)
  
  # extract terms
  mterms <- lm.out$terms
  # aliased coefficients
  aliased <- is.na(coef(lm.out))
  # get model matrix, frame, and response from OLS regression
  X <- model.matrix(lm.out)[,!aliased]
  mf <- model.frame(lm.out)
  y <- model.response(mf)
  # create yX matrix
  yX <- cbind(y, X)
  # extract initial residuals from OLS regression
  original.e <- residuals(lm.out)
  # save variable names
  var.names <- colnames(X)
  # number of observations in regression
  N <- length(y)
  
  # get rank and residual dof (rdf)
  rank <- lm.out$rank
  rdf <- N-rank
  
  # lm automatically row-deletes missing data 
  # so we need to do same for panelVar and timeVar
  # store na.action
  obs.dropped <- lm.out$na.action
  if (!is.null(obs.dropped)) {
    data <- data[-obs.dropped, ]
    panel.vec <- panel.vec[-obs.dropped]
    time.vec <- time.vec[-obs.dropped]
  }
  
  # generate time sequentially if seq.times
  if (seq.times) {
    time.vec <- as.vector(unlist(by(data, panel.vec, function(x) 1:nrow(x))))
    data[, timeVar] <- time.vec # not sure if we need this
  }
  
  # sort units and times
  units <- sort(unique(panel.vec))
  times <- sort(unique(time.vec))
  # number of units and times
  N.units <- length(units)
  N.times <- length(times)
  # average number of units per panel
  N.avgperpanel <- N/N.units
  
  # stop if just 1 time period
  if(N.times<2){
    stop("More than two time periods required.",call.=FALSE)
  }
  
  # check on Parks method
  if(panelCorr.method=="parks" & (N.units > N.times)){
    stop("Cannot estimate Parks-Kmenta method because of singularity.")
  }
  
  # expected number of observations if panels balanced
  NT <- N.units * N.times
  # check if balanced
  balanced <- ifelse(N == NT, TRUE, FALSE)
  
  ### create observations matrix: rows are units and columns are times
  # dimensions: N_p x T
  # if cell is TRUE, panel i at time t is observed
  obs.mat <- reshape(cbind(data[, c(panelVar, timeVar)], TRUE), timevar = timeVar, 
                     idvar = panelVar, direction = "wide", new.row.names = units)[,-1]
  col.order <- order(as.integer(gsub("TRUE.", "", colnames(obs.mat))))
  obs.mat <- obs.mat[, col.order]
  colnames(obs.mat) <- times
  obs.mat[is.na(obs.mat)] <- FALSE
  obs.mat <- as.matrix(obs.mat)
  
  ### runs analysis: generate number of runs by panel
  # output is vector of counts of number of runs per panel
  # length of vector is N_p
  N.runs.panel <- apply(obs.mat, MARGIN = 1, function(x) sum(rle(x)$values == TRUE))
  if (any(N.runs.panel > 1)) {
    message(paste("The following units have non-consecutive observations. Use runs.analysis() on output for additional details: ", paste(names(N.runs.panel[N.runs.panel>1]),collapse=", "),".",sep=""))
  }
  
  ### reshape residuals into T x N_p matrix
  # matrix will have NAs if unbalanced design
  e.mat <- matrix(NA, nrow = ncol(obs.mat), ncol = nrow(obs.mat))
  e.mat[t(obs.mat) == TRUE] <- original.e
  
  ### Run Prais-Winsten correction
  pw.output <- prais.correct(method = autoCorr.method,env.base=env.base)
  
  ### Extract residuals from Prais-Winsten regression
  transformed.resids <- pw.output$pw.lm$residuals
  ### Extract model matrix and response from Prais-Winsten regression
  model.mat.pw <- model.matrix(pw.output$pw.lm)
  
  ### In cases where rhos not bounded to [-1,1], we lose first observation
  ### so need to remake observation matrix
  # observation matrix is now T x N_p
  obs.mat.pw <- t(obs.mat)
  
  # create new observation matrix and update panel and time ID vectors 
  if (!is.null(pw.output$pw.lm$na.action)) {
    panel.vec.pw <- panel.vec[-pw.output$pw.lm$na.action]
    time.vec.pw <- time.vec[-pw.output$pw.lm$na.action]
    obs.mat.pw[obs.mat.pw == TRUE][pw.output$pw.lm$na.action] <- FALSE
  } else {
    panel.vec.pw <- panel.vec
    time.vec.pw <- time.vec
  }
  
  ### Panel homoskedasticity
  if (panelCorr.method == "none") {
    # calculate residual variance
    sigma <- mean(transformed.resids^2)
    # N_p x N_p matrix of panel covariances
    Sigma <- diag(sigma, N.units)
    # matrix of residual covariances
    Omega <- diag(sigma, nrow(model.mat.pw))
    # number of panel covariances calculated
    N.cov <- 1
    # run estimation
    res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
    ### Panel heteroskedasticity
  } else if (panelCorr.method == "phet"|panelCorr.method == "pwls") {
    # N_p x N_p matrix of panel covariances
    sigma.vec <- as.vector(by(transformed.resids, panel.vec.pw, function(x) mean(x^2)))
    if(length(sigma.vec)>1){
      Sigma <- diag(sigma.vec)
    } else{
      Sigma <- sigma.vec
    }
    # matrix of residual covariances
    Omega <- diag(rep(sigma.vec, times = as.integer(table(panel.vec.pw))))
    # number of panel covariances calculated
    N.cov <- length(unique(panel.vec.pw))
    # run estimation
    res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
  } else {
    
    ### Correlated panels / PCSE
    # if balanced panel design
    if (balanced) {
      # set up new T x N_p matrix for residuals
      E <- matrix(transformed.resids, nrow = N.times, ncol = N.units, byrow = FALSE)
      # E'E
      E.E <- crossprod(E)
      # N_p x N_p matrix which 
      # gives number of overlapping time observations between each panel
      weight.mat <- crossprod(obs.mat.pw)
    } else {
      
      # set up new T x N_p matrix for residuals
      E <- obs.mat.pw
      # set missing residuals to 0
      E[E == TRUE] <- transformed.resids
      E[E == FALSE] <- 0
      
      # complete case calculation of E'E and weight.mat
      if (complete.case) {
        # average times per panel
        # select only those time periods with complete cases
        I.com.case <- apply(obs.mat.pw, MARGIN = 1, function(x) prod(x) == 1)
        if (!any(I.com.case)) {
          stop("Unable to compute correlated SEs / PCSEs because there are no time periods in common across all units. Instead, consider setting complete.case=FALSE.",call.=FALSE)
        } else {
          if(sum(I.com.case)<(0.5*N.avgperpanel)){
            warning(paste("The number of time periods used for the calculation of correlated SEs / PCSEs (",as.character(sum(I.com.case)),") is less than half the average number of time periods per panel (",as.character(round(N.avgperpanel,digits=2)),"). Consider setting complete.case=FALSE.",sep=""),call.=FALSE)
          }
          E[!I.com.case, ] <- 0
          E.E <- crossprod(E)
          weight.mat <- matrix(data = sum(I.com.case), nrow = nrow(E.E), ncol = nrow(E.E))
        }
      } else {
        # pairwise option
        E.E <- crossprod(E)
        weight.mat <- crossprod(obs.mat.pw)
      }
    }
    
    # N_p x N_p matrix of panel covariances
    Sigma <- E.E/weight.mat
    # Number of panel covariances calculated
    N.cov <- length(Sigma[lower.tri(Sigma, diag = TRUE)])
    # assume covariance is 0 between panels that do not overlap in time periods
    Sigma <- replace(Sigma, is.na(Sigma), 0)
    
    # matrix of residual covariances
    Omega <- kronecker(Sigma, diag(1, N.times))
    if (!balanced) {
      Omega <- Omega[as.vector(obs.mat.pw), as.vector(obs.mat.pw)]
    }
    
    # run estimation
    res <- switch(pMethod,OLS=ols(env.base), GLS=gls(env.base))
  }
  # end of correlated SE calculations    
  
  # clean up coef vector and var-covar matrix
  coef <- as.vector(res$coef)
  vcov <- res$vcov
  colnames(vcov) <- rownames(vcov) <- names(coef) <- var.names
  
  # calculate yhat
  yhat <- as.vector(X %*% coef)
  names(yhat) <- row.names(X)
  
  # calculate residuals as y-yhat
  resids <- y - yhat
  
  # dof correction to vcov
  if (dof.correction) {
    vcov <- vcov * (N/(N - rank))
  }
  
  # create panelStructure list
  panelStructure <- list(obs.mat = obs.mat, rho = pw.output$pw.rho, Sigma = Sigma, N.cov=N.cov)
  
  if(autoCorr.method=="psar1"){
    names(panelStructure$rho) <- rownames(obs.mat)
  }
  
  # R2 (if ols)
  if(pMethod=="OLS"){
    transform.y.vec <- model.response(model.frame(pw.output$pw.lm))
    r2 <- 1 - sum(transformed.resids^2)/sum((transform.y.vec-mean(transform.y.vec))^2)
  } else{
    r2 <- NULL
  }
  
  # insert panelVar and timeVar into model frame for output
  mf[,panelVar] <- panel.vec
  mf[,timeVar] <- time.vec
  
  # create output list
  fit <- list(coefficients = coef, residuals = resids,  fitted.values = yhat, rank = rank,
              df.residual = rdf,call = call,terms = mterms, model = mf,
              aliased=aliased, na.action = obs.dropped, vcov=vcov, 
              r2=r2, panelStructure = panelStructure)
  class(fit) <- "panelAR"
  return(fit)
}  # end of panelAR fxn

# EST RHO

### Function to estimate rho
### Author: Konstantin Kashin
### August 1, 2013

# input: vector of residuals for a given panel
# k is rank
# output: correlation coefficient
# function is robust to panels with 1 time period
est.rho <- function(e,k,rhotype){
  mat <- embed(e,2)
  
  if(rhotype %in% c("breg","freg")){
    if(rhotype=="breg"){
      # backward regression
      # col 1 is resid, col 2 is lag
      num <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
      denom <- sum((mat[!is.na(mat[,1]),2])^2,na.rm=TRUE)
      rho <- num/denom
    } else {
      #forward regression
      # col 2 is resid, col 1 is forward
      num <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
      denom <- sum((mat[!is.na(mat[,2]),1])^2,na.rm=TRUE)
      rho <- num/denom
    }
  } else if(rhotype %in% c("dw","theil-nagar")){
    sse <- sum(e^2,na.rm=TRUE)
    sseN <- length(na.omit(e))
    dwal <- sum((mat[,1]-mat[,2])^2, na.rm=TRUE)/sse
    
    if(rhotype=="dw"){
      # Durbin-Watson calculation 
      rho <- 1-dwal/2
    } else{
      # Theil-Nagar
      rho <- (sseN^2 * (1-dwal/2) + k^2)/(sseN^2 - k^2)
    }	
  } else{
    sse <- sum(e^2,na.rm=TRUE)
    sseN <- length(na.omit(e))
    cov <- sum(apply(mat,MARGIN=1,prod),na.rm=TRUE)
    # scorr
    rho <- cov/sse
    if(rhotype=="theil"){
      # Theil rho
      # scale scorr by (T-k)/(T-1)
      # Ncov is T
      # Ncov is sseN-1 = T-1
      Ncov <- length(na.omit(apply(mat,MARGIN=1,prod)))
      rho <- (sseN - k)/Ncov*rho
    }
  }
  rho
}

### estim 
### Functions for OLS and GLS Estimation
### Author: Konstantin Kashin
### August 1, 2013

ols <- function(env.base){
  coef <- coef(get("pw.output",envir=env.base)$pw.lm)
  model.mat.pw <- get("model.mat.pw",envir=env.base)
  Omega <- get("Omega",envir=env.base)
  sandwich <- t(model.mat.pw) %*% Omega %*% model.mat.pw
  vcov <- solve(t(model.mat.pw) %*% model.mat.pw) %*% sandwich %*% solve(t(model.mat.pw) %*% model.mat.pw)
  out <- list(coef=coef,vcov=vcov)
  return(out)	
}

gls <- function(env.base){
  model.mat.pw <- get("model.mat.pw",envir=env.base)
  Omega <- get("Omega",envir=env.base)
  response <- model.response(model.frame(get("pw.output",envir=env.base)$pw.lm), "numeric")    
  coef <- solve(t(model.mat.pw) %*% solve(Omega) %*% model.mat.pw) %*% t(model.mat.pw) %*% solve(Omega) %*% response
  vcov <- solve(t(model.mat.pw) %*% solve(Omega) %*% model.mat.pw)
  out <- list(coef=coef,vcov=vcov)
  return(out)	
}

## plot_panelAR

plot.panelAR <- function(x,legend=TRUE,rot.axis=c(0,0),...){
  if(class(x)!="panelAR"){
    stop("x must be of class 'panelAR'.")
  }
  obs.mat <- x$panelStructure$obs.mat
  N.time <- ncol(obs.mat)
  N.panel <- nrow(obs.mat)
  times <- colnames(obs.mat)
  units <- rownames(obs.mat)
  image(z=t(obs.mat), ylab="", xlab="Time", axes=FALSE, col=c("darkred","wheat"),...)
  axis(1, at = seq(from = 0, to = 1, 
                   length = N.time), tck = 0, labels=FALSE, lwd = 0, las = 2)
  text(seq(from = 0, to = 1, 
           length = N.time), par("usr")[3], labels = times, srt = rot.axis[1], pos = 1, 
       xpd =TRUE, cex=par("cex.axis"))
  axis(2, labels = FALSE, at = seq(from = 0, 
                                   to = 1, length = N.panel), tck = 0, lwd = 0, 
       las = 1)
  text(par("usr")[1],seq(from = 0, to = 1, 
                         length = N.panel), labels =units, pos = 2, xpd = TRUE, cex=par("cex.axis"), 
       srt=rot.axis[2])
  if (legend) {
    par(xpd = TRUE)
    legend(x = 0.95, y = 1.1, col = c("darkred","wheat"), bty = "n", 
           xjust = 1, legend = c("Missing", "Observed"), 
           fill = c("darkred","wheat"), horiz = TRUE)
  }	
}

## prais transform

### Function to do prais-winsten transformation
### Author: Konstantin Kashin
### August 1, 2013
prais.transform <- function(i,rhos,p,N.times,units,panel.vec,yX,obs.mat){
  mat <- matrix(NA, nrow=N.times, ncol=p)
  unit.i <- units[i]
  rho.i <- rhos[i]
  mat[obs.mat[i,],] <- yX[panel.vec %in% unit.i,]
  
  mat.L.mat <- rbind(c(mat[1,],rep(NA,p)),embed(mat,2))
  mat.diff <- mat.L.mat[,1:p]-rho.i*mat.L.mat[,(p+1):(2*p)]
  
  # Prais correction for start of runs
  begin.run <- which(!is.na(mat.L.mat[,1]) & is.na(mat.L.mat[,(p+1)]))
  mat.diff[begin.run,] <- mat[begin.run,]*(1-rho.i^2)^0.5
  
  # remove missing values
  keep.index <- sort(union(begin.run,which(complete.cases(mat.diff))))
  mat.diff <- mat.diff[keep.index,]
}

## predict panel AR

predict.panelAR <- function(object,newdata=NULL,se.fit = FALSE,
                            conf.interval = FALSE, conf.level = 0.95, na.action=na.pass,...){
  beta <- object$coefficients
  tt <- terms(object)
  terms <- delete.response(tt)
  df <- object$df
  
  if (missing(newdata) || is.null(newdata)) {
    yhat <- object$fitted.values
    X <- model.matrix(terms,out$model)
  } else{
    mf <- model.frame(terms, newdata, na.action = na.action)
    # pull class descriptions of variables and see if they match the model frame
    .checkMFClasses(attr(terms, "dataClasses"), mf)
    X <- model.matrix(terms, mf)
    yhat <- X %*% beta
  }	
  
  if(se.fit || conf.interval){
    fit.out <- data.frame(fit=yhat)
    d <- diag(X%*%object$vcov%*%t(X))
    d[d<0] <- NA
    se <- sqrt(d)
    if(se.fit){
      fit.out$se <- se		
    }
    if(conf.interval){
      fit.out$lb <- yhat + se*qt((1-conf.level)/2,df=df)
      fit.out$ub <- yhat + se*qt(1-(1-conf.level)/2,df=df)
    }
  } else{
    fit.out <- yhat
  }
  
  out <- list(fit=fit.out,df=df)
  out
}


### print panel ar

print.panelAR <- function(x, digits=max(3,getOption("digits")-3),...){
  if(x$call$autoCorr=="none"){
    autoCorr.Method <- "no autocorrelation"
  } else{
    autoCorr.Method <- "AR(1) Prais-Winsten correction"
  }
  panelCorr.Method  <- switch(x$call$panelCorrMethod,none="homoskedastic variance",phet="panel heteroskedasticity-robust standard errors",pwls="panel weighted least squares",pcse="panel-corrected standard errors",parks="Parks-Kmenta FGLS")
  
  cat(paste("\nPanel Regression with ",autoCorr.Method, " and ",panelCorr.Method,"\n", sep = "")) 
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  if(any(x$aliased)){
    coef <- rep(NA,length(x$aliased))
    names(coef) <- names(x$aliased)
    coef[!x$aliased] <- coef(x)
  } else{
    coef <- coef(x)
  }
  cat("Coefficients:\n")
  print.default(format(coef, digits = digits), print.gap = 2, 
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#print summary panel ar

print.summary.panelAR <- function(x,digits = max(3, getOption("digits") - 3), 
                                  signif.stars = getOption("show.signif.stars"),...){
  if(x$call$autoCorr=="none"){
    autoCorr.Method <- "no autocorrelation"
  } else{
    autoCorr.Method <- "AR(1) Prais-Winsten correction"
  }
  panelCorr.Method  <- switch(x$call$panelCorrMethod,none="homoskedastic variance",phet="panel heteroskedasticity-robust standard errors",pwls="panel weighted least squares",pcse="panel-corrected standard errors",parks="Parks-Kmenta FGLS")
  
  table.structure <- cbind(c("Total obs.:","Number of panels:","Number of times:"),
                           c(x$panelStructure$N,x$panelStructure$N.panel,x$panelStructure$N.time),
                           c("Avg obs. per panel","Max obs. per panel","Min obs. per panel"),
                           c(round(x$panelStructure$N.avg,digits),x$panelStructure$N.max,x$panelStructure$N.min))
  dimnames(table.structure) <- list(rep("",nrow(table.structure)),rep("",ncol(table.structure)))
  
  cat(paste("\nPanel Regression with ",autoCorr.Method, " and ",panelCorr.Method,"\n", sep = ""))
  cat(paste("\n",ifelse(x$panelStructure$balanced,"Balanced","Unbalanced")," Panel Design:",sep=""))
  print.default(table.structure, quote=F,print.gap=1)
  
  if (any(x$aliased)) {
    cnames <- names(x$aliased)
    coefs <- matrix(NA, length(x$aliased), 4, dimnames = list(cnames, 
                                                              colnames(x$coefficients)))
    coefs[!x$aliased, ] <- x$coefficients
  } else{
    coefs <- x$coefficients
  }
  cat("\nCoefficients:\n")
  printCoefmat(coefs,digits = digits, signif.stars = signif.stars, 
               na.print = "NA", ...)
  if(!is.null(x$r2)){cat(paste("\nR-squared: ",round(x$r2,4),sep=""))}
  cat(paste("\nWald statistic: ",round(x$wald["value"],4),", Pr(>Chisq(",x$wald["df"],")): ",round(x$wald["Pr(>Chisq)"],4),"\n",sep=""))
}

## run summary

getRuns <- function(x,times){
  runs <- rle(x)
  runs.length <- runs$lengths[runs$values==TRUE]
  time.end.run <- names(runs$values[runs$values==TRUE])
  time.start.run <- times[which(times %in% time.end.run)-runs.length+1]
  run.table <- data.frame(Start=as.integer(time.start.run),End=as.integer(time.end.run),stringsAsFactors=FALSE)
  return(run.table)
}

run.analysis <- function(object){
  if(class(object)!="panelAR"){
    stop("object must be of class 'panelAR'.",call.=FALSE)
  }
  obs.mat <- object$panelStructure$obs.mat
  times <- colnames(obs.mat)
  units <- rownames(obs.mat)
  
  run.list <- apply(obs.mat,MARGIN=1,function(x) getRuns(x,times=times))
  run.count <- sapply(run.list,function(x) nrow(x))
  run.table <- do.call(rbind,run.list)
  run.table$Length <- run.table$End-run.table$Start+1
  run.table <- as.matrix(run.table)
  rownames(run.table) <- rep(units,run.count)
  
  
  out <- list(run.count=run.count, runs=run.table,rho=object$panelStructure$rho)
  class(out) <- "panelAR.runs"
  out
}

print.panelAR.runs <- function(x,...){
  count.table <- data.frame(Unit=names(x$run.count),Runs=x$run.count)
  
  if(any(x$run.count>1) & !is.null(x)){
    names <- paste(names(x$run.count[which(x$run.count>1)]), collapse=", ")
    cat(paste("Calculation of autocorrelation coefficient restarted for each run for the following panels: ",names,"\n\n" ,sep=""))
  }
  cat("Run Counts:\n")
  print(count.table, row.names=FALSE)
}

# summary.panelAR

summary.panelAR <- function(object,...){
  rdf <- object$df.residual
  rank <- object$rank
  N <- length(object$residuals)
  k <- length(object$aliased)
  df <- c(rank,rdf,k)
  
  # SE
  se <- sqrt(diag(object$vcov))
  
  # test statistics
  coef <- object$coefficients
  t.stat <- (coef)/se
  p.val <- 2*pt(abs(t.stat), rdf, lower.tail=FALSE)
  tab <- cbind(coef,se,t.stat,p.val)
  dimnames(tab) <- list(names(coef), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  
  # set up hypothesis for wald test
  hyp <- names(coef)
  hyp <- hyp[hyp!="(Intercept)"]
  
  # wald test
  lh <- linearHypothesis(object, hyp, test=c("Chisq", "F"), vcov.=object$vcov, singular.ok=FALSE)
  wald <- c(lh$Chisq[2],lh$Df[2],lh[["Pr(>Chisq)"]][2])
  names(wald) <- c("value","df","Pr(>Chisq)") # wald stat, model dof, p stat
  
  # check if balanced
  N.panel <- nrow(object$panelStructure$obs.mat)
  N.time <- ncol(object$panelStructure$obs.mat)
  balanced <- ifelse(N.panel*N.time==N,T,F)
  
  # calculate number of obs per panel, number of time observations, balanced vs. unbalanced 
  N.per.panel <- rowSums(object$panelStructure$obs.mat)
  N.min <- min(N.per.panel)
  N.max <- max(N.per.panel)
  N.avg <- N/N.panel # average units per panel # 
  
  # create list with variables that describe panel structure
  panelStruct <- list(N=N,N.panel=N.panel,N.time=N.time,balanced=balanced,N.min=N.min,N.max=N.max,N.avg=N.avg,N.per.panel=N.per.panel)
  
  out <- list(call=object$call,terms=object$terms,coefficients=tab,residuals=object$residuals, aliased=object$aliased, df=df, rho = object$panelStructure$rho, Sigma=object$panelStructure$Sigma, r2=object$r2, wald=wald, vcov=object$vcov, na.action=object$na.action,panelStructure=panelStruct)
  
  class(out) <- "summary.panelAR"
  out
}

#vcov panelar

# vcov() method for panelAR
vcov.panelAR <- function(object,...){
  object$vcov
}