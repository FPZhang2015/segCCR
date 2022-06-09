

#' @title test the existence of a change point in segCCR
#'
#' @rdname segCCR_test
#'
#'
#' @param dat
#' an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model. 
#' @param tm
#' The vector of tm.
#' @param NB
#' The replicate times to simulate the p-values.
#' 
#' @description Test the existence of a change in segCCR.
#'
#' @details Please refer to \strong{Zhang, F. and Li, Q. (2022)}.
#'
#' @references  Zhang, F. and Li, Q. (2022).
#' Segmented correspondence curve regression for quantifying 
#' covariate effects on the reproducibility of high-throughput experiments.
#'
#' @return A list with the elements:
#'
#' \describe{
#'   \item{pv}{the p-value.}
#'   \item{Tn}{the Tn test statistic.}
#'   \item{Tn.rw}{the Tn test statistic by random weighted method.}
#'   }
#'
#' @author Feipeng Zhang and Qunhua Li
#'
#' @keywords segCCR_test

#' @import stats
#' @importFrom maxLik maxLik
#' @export
#'
#' @examples
#'## The example of ChIP-seq data
#'\dontrun{
#' data(ChIPseq)
#' ## test
#' m = 100
#' tm <- seq(0.01, 0.999, length.out = m)
#' fit = segCCR_test(ChIPseq,
#'            tm=tm,
#'            NB = 10)
#' }


segCCR_test <- function(dat, tm, NB = 100){
  
  y1 <- dat[,1]
  y2 <- dat[,2]
  x <- dat[,3]
  x.unique <- levels(factor(x)) 
  nx <- nlevels(factor(x))
  
  pv = c()
  
  for (k in 1:nx ){
    y1.k = y1[x==as.numeric(x.unique)[k]]
    y2.k = y2[x==as.numeric(x.unique)[k]]
    QLR = segCCR.test(y1.k, y2.k,  tm,  NB)
    pv[k] = QLR$pv
  }
  
  return(pv)
}


segCCR.test <- function(y1, y2,  tm,  NB){
  
  
  ## some preliminary functions
  gfun <- function(x){log(x+1e-10)}
  gfun.inv <- function(x){exp(x)}
  gfun.inv.dev1 <- function(x){exp(x)}
  gfun.inv.dev2 <- function(x){exp(x)}
  
  hfun <- function(x){log(x+1e-10)}
  hfun.dev1 <- function(x) {1/(x+1e-10)}
  hfun.dev2 <- function(x) {-1/(x+1e-10)^2}
  
  ### calculation Uim and correpondence curve
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
  
  ###############################################
  par.ini = c(2, 1)
  par.ini.test = c( 0.5, par.ini)  # initial value
  
  n = length(y1)

  wt <- prepare.Uim(cbind(y1, y2), tm)
  ht <- hfun(tm)
  
  ## H1 
  likFun <- function(thets){
    tau = thets[1]
    bets = thets[-1]
    
    htau <- hfun(tau) 
    zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
    zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
    zz <- cbind(zz1, zz2)
    coef <- drop(zz %*% bets)
    
    ght <- gfun.inv(coef)
    dght <- c(ght[1], diff(ght))
    lpr <- ifelse(dght > 0, log(dght), Inf) 
    lik <- sum(drop(wt %*% lpr))
    
    
    return(lik)
  }
  
  ## H0 
  likFun0 <- function(bet){
    
    zz <- ht
    ht <-   bet*zz
    
    ght <- gfun.inv(ht)
    dght <- c(ght[1], diff(ght))
    lpr <- ifelse(dght > 0, log(dght), Inf) 
    lik <- sum(drop(wt %*% lpr))
    
    return(lik)
  }
  
  
  
  ## original data
  lr <- function(){
    
    ### estimated by MLE 
    ## optimize the objective function by maxLik
    fit <- try(
      maxLik(
        start = par.ini.test,
        logLik = function(thets){likFun(thets)},
        method = "NR", 
        control = list(tol=1e-4, iterlim=50)
      ),
      silent = TRUE)
    
    if(isTRUE(class(fit)=="try-error")) { 
      print("error") 
      est <- rep(NA, length(par.ini))
      llik <- NA
    } else {
      est <- coef(fit)
      llik <- fit$maximum
    }
    
    
    
    ## optimize the objective function by maxLik
    fit0 <- try(
      maxLik(
        start = 1,
        logLik = function(bet){likFun0(bet)},
        method = "NR", 
        control = list(tol=1e-4, iterlim=50)
      ),
      silent = TRUE)
    
    est0 = fit0$est
    llik0 = fit0$maximum
    
    Tn <- llik - llik0
    
    return(list(Tn=Tn, est0=est0, est=est) )
    
  }
  
  fit.org = lr()
  Tn = fit.org$Tn
  est0 = fit.org$est0
  est = fit.org$est
  
  
  ## H1 with random weights
  likFun.rw <- function(thets){
    tau = thets[1]
    bets = thets[-1]
    
    htau <- hfun(tau) 
    zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
    zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
    zz <- cbind(zz1, zz2)
    coef <- drop(zz %*% bets)
    
    ght <- gfun.inv(coef)
    dght <- c(ght[1], diff(ght))
    lpr <- ifelse(dght > 0, log(dght), Inf) 
    lik <- sum(ww*drop(wt %*% lpr))
    
    return(lik)
  }
  
  ## H0 with random weights
  likFun0.rw <- function(bet){
    
    zz <- ht
    ht <-   bet*zz
    
    ght <- gfun.inv(ht)
    dght <- c(ght[1], diff(ght))
    lpr <- ifelse(dght > 0, log(dght), Inf) 
    lik <- sum(ww * drop(wt %*% lpr))
    
    return(lik)
  }
  
  ## random weighted
  lr.rw <- function(ww){
    
    
    
    ### estimated by MLE 
    ## optimize the objective function by maxLik
    fit.rw <- try(
      maxLik(
        start = par.ini.test,
        logLik = function(thets){likFun.rw(thets)},
        method = "NR", 
        control = list(tol=1e-4, iterlim=50)
      ),
      silent = TRUE)
    
    if(isTRUE(class(fit.rw)=="try-error")) { 
      print("error") 
      est <- rep(NA, length(par.ini.test))
      llik.rw <- NA
    } else {
      est.rw <- coef(fit.rw)
      llik.rw <- fit.rw$maximum
    }
    
    
    
    ## optimize the objective function by maxLik
    fit0.rw <- try(
      maxLik(
        start = 1,
        logLik = function(bet){likFun0.rw(bet)},
        method = "NR", 
        control = list(tol=1e-4, iterlim=50)
      ),
      silent = TRUE)
    
    est0.rw = fit0.rw$est
    llik0.rw = fit0.rw$maximum
    
    llik.org = likFun.rw(est)
    llik0.org = likFun0.rw(est0)
    Tn.rw <- llik.rw - llik0.rw - (llik.org - llik0.org)
    
    return(Tn.rw)
    
  }
  
  ## by for
  Tn.rw = c()
  for(b in 1:NB){
    ww = rexp(n, rate = 1)
    Tn.rw[b] = lr.rw(ww)
  }
  

  
  pv = mean(Tn.rw > Tn)
  
  
  rm(list=ls()[! ls() %in% c("pv", "Tn", "Tn.rw")])
  
  return(list(pv = pv, Tn = Tn, Tn.rw = Tn.rw))
}

