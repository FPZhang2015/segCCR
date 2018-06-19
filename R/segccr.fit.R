
segccr.fit <- function(y, x, par.ini, tm, NB=5, sig.level)
{
  if (is.null(n <- nrow(x))) stop("'x' must be a matrix")
  if(n == 0L) stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L) {
    ## ops, null model
    return(list(coefficients = numeric() ))
  }
  
  
  ny <- NCOL(y)
  ## treat one-col matrix as vector
  if(is.matrix(y) && ny != 2)
    stop("y need have 2 columns")
  if (NROW(y) != n)
    stop("incompatible dimensions")
  # chkDots(...)
  
  y1 <- y[,1]
  y2 <- y[,2]
  dat <- cbind(y1, y2, x)
  
  xf <- factor(x[,-1])
  px <- nlevels(xf)
  
  ## some preliminary functions
  gfun <- function(u){log(u+1e-10)}
  gfun.inv <- function(u){exp(u)}
  gfun.inv.dev1 <- function(u){exp(u)} 
  gfun.inv.dev2 <- function(u){exp(u)} 
  
  hfun <- function(u){log(u+1e-10)}
  hfun.dev1 <- function(u) {1/(u+1e-10)}  
  hfun.dev2 <- function(u) {-1/(u^2+1e-10)}  
  
  ################################################
  
  # initial values
  if(is.null(par.ini)){
    par.ini = c(0.5, 2, 1, rep(0.1, 2*px))
  }
  ## initial values for regression coefficients
  bets.ini = par.ini[-1]
  
  ## initial for the change point
  pini = par.ini[1]
  if(length(tm)<=30){
    ptlow <- max(pini - 0.2, 1e-2)
    ptup <- min(pini + 0.2, 1-1e-2)
  }else{
    ptlow <- max(pini - 0.1, 1e-2)
    ptup <- min(pini + 0.1, 1-1e-2)
  }

  taus <- tm[ptlow <= tm & tm <= ptup]   # grid search range for change point
  
  ## tm and psi.tm are declared before
  prepare.Uim <- function(a, tm){ ### Uim
    
    n <- nrow(a)
    m <- length(tm)
    
    # a: (y_i1, y_i2) scores
    a.cdf1 <- ecdf(a[,1])
    a.cdf2 <- ecdf(a[,2])
    cdf1.value <- a.cdf1(a[,1])*n/(n+1)
    cdf2.value <- a.cdf2(a[,2])*n/(n+1)
    
    
    # get counts in each category
    wt.temp <- 1*(1*outer(cdf1.value, tm, "<=") + 1*outer(cdf2.value, tm, "<=")==2) 
    wt <- cbind(wt.temp[, 1], wt.temp[, -1] - wt.temp[, -m])     # n*m
    
    invisible(wt)
  }
  
  ### estimated by grid search method 
  likGrid <- function(dat, tm){
    
    y1 <- dat[,1]
    y2 <- dat[,2]
    x <- as.matrix(dat[,-c(1,2)])
    
    xf <- factor(x[,-1])
   
    likFun <- function(bets, tau){
      
      if(nlevels(xf)==0){
        
        wt <- prepare.Uim(cbind(y1, y2), tm)
        ht <- hfun(tm)
        htau <- hfun(tau) 
        zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
        zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
        coef <-  drop(bets[1]*zz1 + bets[2]*zz2)
        
        ght <- gfun.inv(coef)
        dght <- c(ght[1], diff(ght))
        lpr <- ifelse(dght > 0, log(dght), Inf) 
        lik <- sum(wt %*% lpr)
        
      }else{
        (modmat <- model.matrix(~ xf))
        lables <- apply(modmat, 1, paste, collapse ="|")
        lables_uqi <- unique(lables)
        
        lik <- 0
        for (i in 1:length(lables_uqi)){
          is.in <- which(lables==lables_uqi[i])
          y1.x <- y1[is.in]
          y2.x <- y2[is.in]
          
          x.i <- unlist(strsplit(lables_uqi[i],split = "|"))
          x.i <- as.numeric(x.i[which(!(x.i=="|"))])

          wt <- prepare.Uim(cbind(y1.x, y2.x), tm)
          
          ht <- hfun(tm)
          htau <- hfun(tau) 
          zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
          zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
          zz <- rbind(zz1, zz2)
          
          coef <- drop(x.i%*% matrix(bets, ncol=2, byrow = TRUE) %*% zz)
          ght <- gfun.inv(coef)
          dght <- c(ght[1], diff(ght))
          
          lpr <- ifelse(dght > 0, log(dght), Inf) 
          lik <- lik + sum(wt %*% lpr)
        }
      }
        
      
      
      return(lik)
    }
    
    
    lik <- NULL
    n.tau <- length(taus)
    betcoef <- matrix(0, nrow = n.tau, ncol = 2*px)
    for (kk in 1:n.tau){   
      
      fit <- try(
        optim(
          par = bets.ini, 
          fn = function(bets){-likFun(bets, taus[kk])}, 
          method = "Nelder-Mead", 
          hessian = FALSE),
        silent = TRUE)
      
      if(isTRUE(class(fit)=="try-error")) { NA } 
      else {
        # optimize minimization
        betcoef[kk, ] <- fit$par
        #betse[ , ,kk] <- fit$hessian
        lik[kk] <- (-1)*fit$value
      }
      
    }  
    
    
    ## find the change-point 
    bp <- min(taus[lik == max(lik, na.rm = TRUE)] )  
    
    ## estimate the regression coefficients
    bethat <- betcoef[taus == bp, ]
    #betse <- sqrt(diag(ginv(betse[ , , taus == bp])))
    
    thethat <- c(bp, bethat)
    
    invisible(return(list(thethat= thethat)))
  }
  
  ## by bootstrap with parallel computation
  likGrid.boot <- function(dat, boot.indx){
    
    # bootstrap the data
    #boot.indx <- sample(nrow(dat), replace = T)
    dat <- dat[boot.indx, ]
    
    fit <- likGrid(dat, tm)
    
    thethat <- fit$thethat
    
    invisible(return(thethat))
    
  }
  
  # bootstrap variance estimate in bulk
  eseGrid.boot <- function(dat, tm){
    #require(boot)
    
    thet.boot <- boot(dat, likGrid.boot, R = NB, parallel="multicore",
                      ncpus = getOption("boot.ncpus", 4L))
    thet.ese <- apply(thet.boot$t, 2, sd)
    
    invisible(return(thet.ese))
  }
  
  ## estimate procedure
  fit <- try(likGrid(dat,  tm), silent = TRUE)
  if(isTRUE(class(fit)=="try-error")){
    print("error")
  }else{ 
    thethat <- fit$thethat
  }
  
  
  # the bootstrap standard errors estimate, which is time consuming
  thetese <- eseGrid.boot(dat, tm)
  
  
  ## names 
  est <- thethat
  ese <- thetese
 
  if(px==0){
    dn1 <- c(expression(tau),
             paste0("baseline",":-", seq = ""),
             paste0("baseline", ":+", seq = ""))
  }else{

    dn0 <- levels(xf)[-1]
    dn1 <- c(expression(tau),
             paste0("baseline",":-", seq = ""),
             paste0("baseline", ":+", seq = ""),
             paste0("level", rep(dn0, each=2), 
                    rep(c(":-", ":+"), times=length(dn0)), seq = ""))
  }


  names(est) <- dn1
  
  
  CI.low <- est - qnorm(1-sig.level/2, 0, 1) * ese
  CI.up <- est + qnorm(1-sig.level/2, 0, 1) * ese
  
  z = list()
  z$coefficients <- est
  z$std.error <- ese
  z$CI <- cbind(CI.low, CI.up)
  
  return(c(z[c("coefficients", "std.error", "CI")]))
  
}
