
#' @title Fit the segmented correspondence curve regression
#'
#' @rdname segCCR
#'
#' @param formula
#' object of class 'formula'
#' (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' The details of model specification are given under 'Details'.
#'
#' @param data
#' an optional data frame, list or environment
#' (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model. If not found in data,
#' the variables are taken from environment(formula),
#' typically the environment from which segCCR is called.
#' @param model
#' logical. If TRUE the corresponding components of the fit  are returned.
#' @param par.ini
#' the initial values for the estimate parameters. The first component is
#' the change point. If is.null(par.ini) == TRUE, par.ini is set in the
#' the details.
#' @param tm
#' The vector of tm.
#' @param NB
#' The bootstrap times to obtain the estimated standard errors.
#' @param method
#' The method to be used for fitting, currently only method = 'segCCR_fit'.
#' @param sig.level
#' The significant level. Default is 0.05.
#'
#' @description Fit the correspondence curve regression
#'
#' @details Please refer to \strong{Zhang, F. and Li, Q. (2022)}.
#'
#' @references Zhang, F. and Li, Q. (2022).
#' Segmented correspondence curve regression for quantifying 
#' covariate effects on the reproducibility of high-throughput experiments.
#'
#' @return A list with the elements:
#'
#' \describe{
#'   \item{coefficients}{a named vector of coefficients.}
#'   \item{model}{if requested (the default), the model frame used.}
#'   \item{call}{the matched call.}
#'   \item{std.error}{the estimated standard errors.}
#'   \item{CI}{the confidence intervals.}
#'   \item{pv}{the p-values for individual test.}
#'   \item{jointpv}{the p-values for joint test.}
#'   }
#'
#' @author Feipeng Zhang and Qunhua Li
#'
#' @keywords segCCR

#' @import stats
#' @importFrom boot boot
#' @export
#'
#' @examples
#'## The example of ChIP-seq data
#'\dontrun{
#' data(ChIPseq)
#' ## estimate
#' m = 100
#' tm <- seq(0.01, 0.999, length.out = m)
#' nx = nlevels(factor(ChIPseq$x))
#' par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
#' fit = segCCR(cbind(y1, y2)~x, data = ChIPseq,
#'            par.ini = par.ini,
#'            tm=tm,
#'            NB = 5)
#' }



segCCR_old <- function (formula, data, model = TRUE,
                        par.ini,
                        tm = seq(0.01, 0.999, length.out = 50),
                        NB = 5, sig.level = 0.05,
                        method = 'segCCR_fit')
{
  
  #method <- match.arg(method)
  cl <- match.call()
  
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c('formula', 'data'), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(identical(method, 'model.frame')) return(mf)
  
  # if (!is.character(method) && !is.function(method))
  # stop('invalid 'method' argument')
  
  mt <- attr(mf, 'terms') # allow model.frame to have updated it
  y <- model.response(mf, 'numeric')  # response
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients =  numeric() )
  }
  else {
    ## We want to set the name on this call and the one below for the
    ## sake of messages from the fitter function
    x <- model.matrix(mt, mf, contrasts)
    
    #z <- segccr.fit(y, x, par.ini, tm, NB, sig.level) # which equals to ...
    z <- eval(call(if(is.function(method)) 'method' else method,
                   y, x, par.ini, tm, NB, sig.level))
  }
  
  z$contrasts <- attr(x, 'contrasts')
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)   z$model <- mf
  
  z$x <- x
  z$y <- y
  
  z
}





segCCR_fit <- function(y, x, par.ini, tm, NB=10, sig.level)
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
  
  m <- length(tm)
  
  y1 <- y[,1]
  y2 <- y[,2]
  dat <- cbind(y1, y2, x)
  
  xf <- factor(x[,-1])
  px <- nlevels(xf)   ## without baseline
  
  
  x.unique <- levels(xf)
  nx <- px+1
  
  
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
  
  
  ##### fit segCCR
  segCCR.fit <- function(dat, tm, taus, par.ini){
    
    #### likelihood function
    likFun <- function(bets, tau){
      
      if(px==0){
        
        wt <- prepare.Uim(cbind(y1, y2), tm)
        ht <- hfun(tm)
        htau <- hfun(tau)
        zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
        zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
        zz = cbind(zz1, zz2)
        coef <-  drop(zz %*% bets)
        
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
    betcoef <- matrix(0, nrow = n.tau, ncol = 2*nx)
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
        lik[kk] <- (-1)*fit$value
      }
      
    }
    
    
    ## find the change-point
    bp <- min(taus[lik == max(lik, na.rm = TRUE)] )
    
    ## estimate the regression coefficients
    bethat <- betcoef[taus == bp, ]
    
    thethat <- c(bp, bethat)
    
    return(list(thethat=thethat))
  }
  
  ## estimated coefficients
  fit =  segCCR.fit(dat, tm, taus, par.ini)
  thethat = fit$thethat
  
  ## bootstrap variance estimate
  segccr_cov <- function(dat, tm, taus, par.ini, NB){
    
    ## by bootstrap with parallel computation
    segccr_boot <- function(dat, boot.indx){
      
      # bootstrap the data
      #boot.indx <- sample(nrow(dat), replace = T)
      dat <- dat[boot.indx, ]
      
      fit <- segCCR.fit(dat, tm, taus, par.ini)
      
      thethat <- fit$thethat
      
      invisible(return(thethat))
      
    }
    
    # thet.boot <- boot(dat, segccr_boot, R = NB, parallel="multicore",
    #                   ncpus = getOption("boot.ncpus", n.cls))
    thet.boot <- boot(dat, segccr_boot, R = NB)
    
    thet.cov <- cov(thet.boot$t)
    #thet.ese <- apply(thet.boot$t, 2, sd, na.rm = TRUE)
    
    invisible(return(thet.cov))
  }
  
  
  ## save the bootstrap variance estimate, which is time consuming
  thetcov = segccr_cov(dat, tm, taus, par.ini, NB)
  thetse = sqrt(diag(thetcov))
  
  ## H0: beta_sj=0
  thetpv <- 2*(1-pnorm(abs(thethat/thetse), 0, 1))
  
  
  ## joint test: H0: beta_s1=beta_s2=0
  est.beta = thethat[-1]
  cov.beta = thetcov[-1, -1]
  jtest.pvs = c()
  for (k in 1:nx){
    Rmat.k = cbind(matrix(0, nrow = 2, ncol = 2*(k-1)), diag(2), matrix(0, nrow = 2, ncol = 2*(nx-k) ) )
    Wn.k = drop(t(Rmat.k%*%est.beta) %*% solve(Rmat.k%*%cov.beta%*%t(Rmat.k)) %*% (Rmat.k%*%est.beta))
    # jtest.pvs[k] = 2*min(pchisq(q=Wn.k, df=2, lower.tail=FALSE), pchisq(q=Wn.k, df=2, lower.tail=TRUE))
    jtest.pvs[k] = pchisq(q=Wn.k, df=2, lower.tail=FALSE)
  }
  
  
  
  
  
  ## names
  est <- thethat
  ese <- thetse
  
  if(px==0){
    dn1 <- c(expression(tau),
             expression(beta["01"]),
             expression(beta["02"])
    )
  }else{
    
    dn0 <- levels(xf)[-1]
    dn1 <- c(expression(tau),
             expression(beta["01"]),
             expression(beta["02"]),
             paste0(rep(dn0, each=2),
                    rep(c("1", "2"), times=length(dn0)), seq = ""))
  }
  
  # plot(1,1, main=expression('title'^2))  #superscript
  # plot(1,1, main=expression('title'[2])) #subscript
  
  
  
  names(est) <- dn1
  
  
  CI.low <- est - qnorm(1-sig.level/2, 0, 1) * ese
  CI.up <- est + qnorm(1-sig.level/2, 0, 1) * ese
  
  z = list()
  z$coefficients <- est
  z$std.error <- ese
  z$CI <- cbind(CI.low, CI.up)
  z$pv <- thetpv
  z$jointpv <- jtest.pvs
  
  return(c(z[c("coefficients", "std.error", "CI", 'pv', "jointpv")]))
  
}
