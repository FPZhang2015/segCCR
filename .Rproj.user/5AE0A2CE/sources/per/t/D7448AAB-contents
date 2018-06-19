
#' @title Salmon dataset
#' @name data_salmon
#' 
#' @description A dataset of two different populations is from Simonoff (1990, p 161). 
#' It concerns the size of the annual spawning stock and its
#' production of new catchable-sized fish for 1940 through 1967 
#' for the Skeena river sockeye salmon stock (in thousands of fish). 
#' It has three columns, year, spawners and recruits. 
#' 
#' The variables are as follows:
#'
#' @format A data frame with 28 rows and 3 variables
#' 
#' \describe{
#'   \item{year}{year}
#'   \item{spawners}{the annual spawning stock}
#'   \item{recruits}{the production of new catchable-sized fish}
#' }

#' @source Data is salmon.dat in Simonoff (1990), 
#' which can be downloaded from the book’s website.
#' 
#' @docType data
#' @keywords data_salmon
#' @usage data(data_salmon)
#' 
#' 
#' @references
#' \itemize{
#' \item Simonoff, J. S. (1990).  
#' Smoothing Methods in Statistics, New York: Spinger-Verlag
#'   }
#'   
#'
#' @examples 
#' \dontrun{
#' data(data_salmon)
#' summary(data_salmon)
#' ## estimate
#' par.ini = c(0.5, 1, 2)  # initial value
#' fit = segccrFit(cbind(spawners, recruits)~1, data = data_salmon,
#'            par.ini = par.ini,
#'            tm=seq(0.01, 0.9999, length.out = 28),
#'            NB = 10)
#' }
#' 
  'data_salmon'


