
#' @title Fit the segmented correspondence curve regression
#'
#' @rdname segccrFit
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
#' typically the environment from which ccrFit is called.
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
#' The method to be used for fitting, currently only method = 'ccr.fit'.
#' @param sig.level
#' The significant level. Default is 0.05.
#'
#' @description Fit the correspondence curve regression
#' 
#' @details Please refer to \strong{Zhang, F. and Li, Q. (2018)}.
#' 
#' @references  Zhang, F. and Li, Q. (2018). 
#' Segmented correspondence curve regression model for quantifying 
#' reproducibility of high-throughput experiments. Submitted.
#' 
#' @return A list with the elements:
#' 
#' \describe{
#'   \item{coefficients}{a named vector of coefficients.}
#'   \item{model}{if requested (the default), the model frame used.}
#'   \item{call}{the matched call.}
#'   \item{std.error}{the estimated standard errors.}
#'   \item{CI}{the confidence intervals.}
#'   }
#'   
#' @author Feipeng Zhang and Qunhua Li
#' 
#' @keywords ccrFit

#' @import stats 
#' @importFrom boot boot
#' @export
#'
#' @examples
#' ## simualted example
#' \dontrun{
#' n = 1000
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
#' y1 = dat[,1]
#' y2 = dat[,2]
#' x = dat[,3]
#' 
#' #####  initial point 
#' par.ini <- c(0.8, 2, 1, 0.1, 0.1)
#' fit = segccrFit(cbind(y1, y2)~x, data = dat,
#'            par.ini = par.ini,
#'            tm=seq(0.01, 0.999, length.out = 50), 
#'            NB = 2)
#' }
#' 
#' ## The example of salmon data
#' \dontrun{
#' data(data_salmon)
#' ## estimate
#' par.ini = c(0.5, 1, 2)  # initial value
#' fit = segccrFit(cbind(spawners, recruits)~1, data = data_salmon,
#'            par.ini = par.ini,
#'            tm=seq(0.001, 0.999, length.out = 28), 
#'            NB = 2)
#'  }   
#'  
#'## The example of ChIP-seq data
#'\dontrun{
#' data(data_ChIP)
#' ## estimate
#' par.ini = c(0.5, 2, 1, rep(0.1, 8))  # initial value
#' fit = segccrFit(cbind(y1, y2)~x, data = data_ChIP,
#'            par.ini = par.ini,
#'            tm=seq(0.01, 0.999, length.out = 100),
#'            NB = 2)
#' }



segccrFit <- function (formula, data, model = TRUE,  
                    par.ini, 
                    tm = seq(0.01, 0.999, length.out = 50), 
                    NB = 5, sig.level = 0.05, 
                    method = 'segccr.fit') 
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