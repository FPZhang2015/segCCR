
#' @title Simulated data
#' 
#' @rdname datSim
#' 
#' @description Generate dara from mixture bivariate Gaussian distribution
#' 
#' @param n    sample size
#' @param bp0  change point
#' @param mu01  mean 
#' @param mu02  mean 
#' @param mu11  mean 
#' @param mu12  mean 
#' @param var01  variance 
#' @param var02  variance
#' @param var11  variance
#' @param var12  variance
#' @param rho01  correltation coefficient
#' @param rho02  correltation coefficient 
#' @param rho11  correltation coefficient
#' @param rho12  correltation coefficient 
#' 
#' @return data.frame as (y1, y2, x)
#' @author Feipeng Zhang
#' @keywords myadd
#' @references
#' \itemize{
#' \item Feipeng Zhang (2017).
#' A toy example for creating R package.
#' Qunhua Li's group metting, 12, xxx--xxx.
#' }
#' 
#' @importFrom MASS mvrnorm
#' @export
#' 
#' @examples
#' \dontrun{
#' #' n = 10000
#' bp0 = 0.8
#' ## for x=0
#' rho01 <- 0.0   ## noise,  irreproducible
#' rho02 <- 0.6   ## signal, reproducible
#' mu01 <-  c(0, 0)
#' mu02 <- c(3, 3)
#' var01 <- 1
#' var02 <- 1
#' # for x=1
#' rho11 <- rho01  ## noise,  irreproducible
#' rho12 <- 0.9  ## signal, reproducible
#' mu11 <-  mu01
#' mu12 <-  mu02
#' var11 <- var01
#' var12 <- var02
#' 
#' ## generate data 
#' dat <- datSim(n, bp0,
#' mu01, mu02, var01, var02, rho01, rho02, 
#' mu11, mu12, var11, var12, rho11, rho12)
#' }



### data generation from mixture bivariate Gaussian distribution  
datSim <- function(n, bp0, 
                   mu01, mu02, var01, var02, rho01, rho02, 
                   mu11, mu12, var11, var12, rho11, rho12){
  
  ##  for protocol x=0
  pi01 <- bp0   ## break point 
  pi02 <- 1-pi01
  n0 <- round(n*pi01)
  n1 <- round(n*pi02)
  Sigma01 <- matrix(c(var01, rho01*var01, rho01*var01, var01), 2, 2)
  Sigma02 <- matrix(c(var02, rho02*var02, rho02*var02, var02), 2, 2) 
  dat0 <- rbind(mvrnorm(n0, mu01, Sigma01), mvrnorm(n1, mu02, Sigma02)) 
  x0 <- rep(0, n)
  
  
  ##  for protocol x=1
  pi11 <- pi01   ## break point
  pi12 <- pi02 
  Sigma11 <- matrix(c(var11, rho11*var11, rho11*var11, var11), 2, 2)
  Sigma12 <- matrix(c(var12, rho12*var12, rho12*var12, var12), 2, 2) 
  dat1 <- rbind(mvrnorm(n0, mu11, Sigma11), mvrnorm(n1, mu12, Sigma12)) 
  x1 <- rep(1, n)
  
  
  ## combine dat0 and dat1 
  dat <- rbind(cbind(dat0, x0), cbind(dat1, x1)) 
  dat <- as.data.frame(dat)
  colnames(dat) <- c("y1", "y2", "x")
  invisible(return(dat))
}
