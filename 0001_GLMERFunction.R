vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

kappa.mer <- function (fit,
                       scale = TRUE, center = FALSE,
                       add.intercept = TRUE,
                       exact = FALSE) {
  X <- fit@pp$X
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  nrp <- sum(1 * (nam == "(Intercept)"))
  if (nrp > 0) {
    X <- X[, -(1:nrp), drop = FALSE]
    nam <- nam[-(1:nrp)]
  }
  
  if (add.intercept) {
    X <- cbind(rep(1), scale(X, scale = scale, center = center))
    kappa(X, exact = exact)
  } else {
    kappa(scale(X, scale = scale, center = scale), exact = exact)
  }
}

colldiag.mer <- function (fit,
                          scale = TRUE, center = FALSE,
                          add.intercept = TRUE) {
  ## adapted from perturb::colldiag, method in Belsley, Kuh, and
  ## Welsch (1980).  look for a high condition index (> 30) with
  ## more than one high variance propotion.  see ?colldiag for more
  ## tips.
  result <- NULL
  if (center) 
    add.intercept <- FALSE
  if (is.matrix(fit) || is.data.frame(fit)) {
    X <- as.matrix(fit)
    nms <- colnames(fit)
  }
  else if (class(fit) == "mer" | class(fit) == "glmerMod" ) {
    nms <- names(fixef(fit))
    X <- fit@pp$X
    if (any(grepl("(Intercept)", nms))) {
      add.intercept <- FALSE
    }
  }
  X <- X[!is.na(apply(X, 1, all)), ]
  
  if (add.intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "(Intercept)"
  }
  X <- scale(X, scale = scale, center = center)
  
  svdX <- svd(X)
  svdX$d
  condindx <- max(svdX$d)/svdX$d
  dim(condindx) <- c(length(condindx), 1)
  
  Phi = svdX$v %*% diag(1/svdX$d)
  Phi <- t(Phi^2)
  pi <- prop.table(Phi, 2)
  colnames(condindx) <- "cond.index"
  if (!is.null(nms)) {
    rownames(condindx) <- nms
    colnames(pi) <- nms
    rownames(pi) <- nms
  } else {
    rownames(condindx) <- 1:length(condindx)
    colnames(pi) <- 1:ncol(pi)
    rownames(pi) <- 1:nrow(pi)
  }         
  
  result <- data.frame(cbind(condindx, pi))
  zapsmall(result)
}

maxcorr.mer <- function (fit,
                         exclude.intercept = TRUE) {
  so <- summary(fit)
  corF <- so$vcov
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0 & exclude.intercept) {
    corF <- corF[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  corF[!lower.tri(corF)] <- 0
  maxCor <- max(corF)
  minCor <- min(corF)
  if (abs(maxCor) > abs(minCor)) {
    zapsmall(maxCor)
  } else {
    zapsmall(minCor)
  }
}

fit_glmer_mc = function(out, models, co) {
  library(doParallel)
  cl <- makeCluster(co)
  registerDoParallel(cl) # Do a parallel cluster of co cores
  #if (is.data.frame(data)) 
  #  data = replicate(length(calls), data, simplify = F)
  r <- foreach(i= 1:length(models), .packages=c("lme4"),.errorhandling = "pass") %do% {
    mod <- try( eval(parse(text = models[i])),silent = T) 
    if(class(mod)=="try-error") mod <- NA
    return(mod)
  }
  stopCluster(cl); stopImplicitCluster() # Stop the cluster
  return(r)
}

fit_gamm4_mc = function(out, models, co) {
  library(doParallel)
  cl <- makeCluster(co)
  registerDoParallel(cl) # Do a parallel cluster of co cores
  #if (is.data.frame(data)) 
  #  data = replicate(length(calls), data, simplify = F)
  r <- foreach(i= 1:length(models), .packages=c("gamm4"),.errorhandling = "pass") %do% {
    mod <- try( eval(parse(text = models[i])),silent = T) 
    if(class(mod)=="try-error") mod <- NA else
      b <- summary(mod$gam)
    mod <- b$r.sq
    return(mod)
  }
  stopCluster(cl); stopImplicitCluster() # Stop the cluster
  return(r)
}

# Extract coefficient function
extractCoef <- function(mod,group){
  # Extracts coefficients of with conditional variances from an lme4 object
  # How to from Ben bolker
  # http://stackoverflow.com/questions/26198958/extracting-coefficients-and-their-standard-error-from-lme
  stopifnot(class(mod) =="lmerMod")
  fixed.vars <- diag(vcov(mod))
  ## extract variances of conditional modes
  r1 <- ranef(mod,condVar=TRUE)
  cmode.vars <- t(apply(cv <- attr(r1[[1]],"postVar"),3,diag))
  seVals <- sqrt(sweep(cmode.vars,2,fixed.vars,"+"))
  res <- cbind(coef(mod)[[group]],seVals)
  res2 <- setNames(res[,c(1,3,2,4)],
                   c("int","int_se","slope","slope_se"))
  return(res2)
}