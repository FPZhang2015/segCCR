
#' @title Test the existence of change point in
#' the segmented correspondence curve regression model.
#'
#' @rdname segccrTest
#' 
#' @param y
#' an matrix with two columns. 
#' 
#' @param tm
#' The vector of tm.
#' 
#' @param NB
#' The bootstrap times to obtain p-values.
#' 
#' @param sig.level
#' The significant level. Default is 0.05.
#'
#' @description Test the existence of change point in
#'  the correspondence curve regression. 
#' 
#' @details Please refer to \strong{Zhang, F. and Li, Q. (2018)}.
#' 
#' @references  Zhang, F. and Li, Q. (2018). 
#' Segmented correspondence curve regression model for quantifying 
#' reproducibility of high-throughput experiments.  Submitted.
#' 
#' @return A list with the elements:
#' 
#' \describe{
#'   \item{pv}{p-values.}
#'   }
#'   
#' @author Feipeng Zhang and Qunhua Li
#' 
#' @keywords ccrTest

#' @import stats 
#' @importFrom boot boot
#' @export
#'
#' @examples
#' ## The example of salmon data
#' \dontrun{
#' data(data_salmon)
#' ## estimate
#' segccrTest(cbind(data_salmon$spawners, data_salmon$recruits), 
#'            tm=seq(0.001, 0.999, length.out = 28), 
#'            NB = 2)
#'  }   


segccrTest <- function (y, 
                        tm = seq(0.01, 0.999, length.out = 50), 
                        NB = 5, sig.level = 0.05) 
{
  ny <- NCOL(y)
  ## treat one-col matrix as vector
  if(is.matrix(y) && ny != 2)
    stop("y need have 2 columns")
  
  ##################  some transformation functions  ################
  ## g=h=log()  
  
  gfun <- function(x){log(x+1e-10)}
  gfun.inv <- function(x){exp(x)}
  gfun.inv.dev <- function(x){exp(x)} 
  hfun <- function(x){log(x+1e-10)}
  hfun.dev <- function(x) {1/(x+1e-10)}  

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
  
  
  ### estimated by grid search method under H1: beta1 != beta2
  likGrid <- function(dat, tm, taus){
    
    y1 <- dat[, 1]
    y2 <- dat[, 2]
    
    n <- length(y1)
    m <- length(tm)
    
    likFun <- function(bets, tau){
      
      alp1 = bets[1]
      alp2 = bets[2]
      
      wt <- prepare.Uim(cbind(y1, y2), tm)  # n*m
      coef1 <- alp1 
      coef2 <- alp2 
      
      ht <- hfun(tm)
      htau <- hfun(tau) 
      zz1 <- (ht-htau) * ifelse(tm <= tau, 1, 0)
      zz2 <- htau + (ht-htau) * ifelse(tm > tau, 1, 0)
      
      ght <- gfun.inv(coef1*zz1 + coef2*zz2)
      dght <- c(ght[1], ght[-1] - ght[-length(ght)])
      
      lpr <- ifelse(dght > 0, log(dght), Inf) 
      lik <- sum(wt %*% lpr)
      
      return(lik)
    }
    
    lik <- NULL
    n.tau <- length(taus)
    betcoef <- matrix(0, nrow = n.tau, ncol=2)
    #betse <- array(0, dim = c(4, 4, n.tau))  # no need for ese
    for (kk in 1:n.tau){   
      
      fit <- try(
        optim(
          par = c(1, 2), 
          fn = function(bets){-likFun(bets, taus[kk])}, 
          method = "Nelder-Mead", hessian = FALSE),
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
    likhat <- lik[taus == bp]
    
    thethat <- c(bp, bethat)
    
    invisible(return(list(thethat= thethat, 
                          bet.prof = betcoef, likhat = likhat)))
  }
  
  ### estimate by naive estimate under H0: beta1 = beta2
  naivelikest <- function(dat, tm){
    
    likNaiveFun<- function(dat, tm, bet){
      
      bet1 = bet[1]
      y1 <- dat[, 1]
      y2 <- dat[, 2]
      
      wt <- prepare.Uim(cbind(y1, y2), tm)    # n * m
      
      zz<- hfun(tm)
      ht<-   bet1*zz
      
      
      ght <- gfun.inv(ht)
      dght <- c(ght[1], diff(ght))
      lpr <- ifelse(dght > 0, log(dght), Inf) 
      lik <- sum(wt %*% lpr)
      
      return(lik)
    }
    
    fit <- optim(par=c(1), 
                 fn=function(bet){-likNaiveFun(dat, tm, bet)}, 
                 method="Nelder-Mead", hessian = TRUE ) 
    
    bet.naive <- fit$par
    bet.se <- sqrt(diag(solve(fit$hessian)))
    lik <- (-1) * fit$value
    return(list(est=bet.naive, ese = bet.se, lik = lik))
  }
  
  #####################################################
  ## test the existence of change point
  
  ## QTR: quasi-likelihood ratio
  QLR <- function(dat, tm, taus){
    
    n <- nrow(dat) 
    
    ## H0: beta1 = beta2
    fit.null <- naivelikest(dat, tm)
    bet.null <- fit.null$est
    lik.null <- fit.null$lik
    
    ## under H1: beta1 != beta2
    fit.grid <- likGrid(dat, tm, taus)
    thet.grid <- fit.grid$thethat
    bet.prof <- fit.grid$bet.prof
    lik.grid <- fit.grid$likhat
    
    ## find the change-point 
    Tn <- n*(lik.grid - lik.null)
    
    return(list(Tn = Tn, bet.null = bet.null, 
                bet.prof = bet.prof, thet.grid = thet.grid))
  }
  
  QLR.boot <- function(dat, tm, taus, bet.null, bet.prof){
    
    y1 <- dat[, 1]
    y2 <- dat[, 2]
    n <- length(y1)
    m <- length(tm)
    
    wt <- prepare.Uim(cbind(y1, y2), tm)  # n*m 
    
    ht <- hfun(tm)            # m*1
    vv <- rnorm(n, 0, 1)    # N(0,1), n*1
    
    ## H0 
    betZt1 <- bet.null * ht                # m*1
    ginv1 <- gfun.inv(betZt1)            # m*1
    ginvdev1 <- gfun.inv.dev(betZt1) * ht       # m*1
    ght1 <- c(ginvdev1[1]/ginv1[1]*(ginv1[1] != 0), 
              diff(ginvdev1)/(diff(ginv1)+1e-8) * (diff(ginv1) != 0))  # m*1
    
    sco1 <- wt %*% ght1            # n*1
    sco.vv1 <- n^(-1/2) * sum(sco1*vv)   # 1*1
    
    Imat1 <- mean(t(sco1) %*% sco1)      # 1*1
    TT1 <- sco.vv1 * 1/Imat1 * sco.vv1
    
    ## H1
    TT <- NULL
    
    for (kk in 1:length(taus)){
      
      tau <- taus[kk]
      htau <- hfun(tau)
      
      Zt <- cbind((ht-htau) * ifelse(tm <= tau, 1, 0), 
                  htau + (ht-htau) * ifelse(tm > tau, 1, 0) )  # m*2
      betZt <- Zt %*% bet.prof[kk,]                    # m*1
      ginv <- gfun.inv(betZt) %*% matrix(1, 1, 2)          # m*2
      ginvdev <- (gfun.inv.dev(betZt) %*% matrix(1, 1, 2)) * Zt      # m*2
      ght <- rbind(ginvdev[1,]/(ginv[1,]+1e-8) * (ginv[1, ] != 0), 
                   apply(ginvdev, 2, diff)/(apply(ginv, 2, diff)+1e-8) * 
                     (apply(ginv, 2, diff) != 0) )    # m*2
      sco <- wt %*% ght        # n*2
      sco.vv <-  n^(-1/2) * as.vector(vv %*% sco)     # 1*2
      Imat <- n^(-1)* t(sco) %*% sco      # 2*2
      TT[kk] <- sco.vv %*% ginv(Imat+1e-8) %*% sco.vv
      
      
    }
    
    Tn <- 1/2 * (max(TT[-c(1, length(TT))]) - TT1)
    return(Tn)
    
  }
  
  ###############################################################
  ## calculate p-value
  # grid search range for change point
  taus <- tm[quantile(tm, 0.1)<tm & tm < quantile(tm, 0.9)] 
  
  fit.QLR <-  QLR(y, tm, taus)
  bet.null <- fit.QLR$bet.null
  bet.prof <- fit.QLR$bet.prof
  
  Tn.QLR <- fit.QLR$Tn
  Tn.QLR.NB <- replicate(NB, QLR.boot(y, tm, taus, bet.null, bet.prof))
  
  pv <- mean(Tn.QLR.NB > Tn.QLR,  na.omit=T) 
  
  return(pv)
}
  
  

  
  