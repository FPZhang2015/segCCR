
#' @title ChIP-seq dataset
#' @name ChIPseq
#' 
#' @description ChIP-seq data from Chen et al. (2012) 
#' to characterize the effect of sequencing depth on the reproducibility 
#' of binding site identification. 
#' At the sequencing depth of 0.45, 0.9, 2.7, 5.4 and 16.2 million reads, 
#' 1335, 2198, 3813, 4631, and 5499 binding sites are identified 
#' on both replicates, respectively.
#' It has three columns, y1, y2 and x. 
#' 
#' The variables are as follows:
#'
#' @format A data frame with 17476 rows and 3 variables
#' 
#' \describe{
#'   \item{y1}{the scores of replicate 1}
#'   \item{y2}{the scores of replicate 2}
#'   \item{x}{the factor variabel of 5 depths. 0 for baseline at depth 0.45M, 
#'   and 1,2,3,4 are for the for the depths of 0.9M, 2.7M, 5.4M and 16.2M, respectively.}
#' }

#' @source Data is from Chen et al. (2012). 
#' 
#' @docType data
#' @keywords ChIPseq
#' @usage data(ChIPseq)
#' 
#' 
#' @references
#' \itemize{
#' \item Chen, Y., Negre, N., Li, Q., Mieczkowska, J. O., Slattery, 
#' M., Liu, T., Zhang, Y., Kim, T.-K., He, H. H., Zieba, J., et al. (2012). 
#' Systematic evaluation of factors influencing ChIP-seq fidelity. 
#' Nature Methods, 9:609鈥?14.
#'   }
#'   
#'
#' @examples 
#' \dontrun{
#' data(ChIPseq)
#' ## estimate
#' m = 100
#' tm <- seq(0.01, 0.999, length.out = m)
#' nx = nlevels(factor(ChIPseq$x))
#' par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
#' nx = nlevels(factor(ChIPseq$x))
#' par.ini = c(0.5, 2, 1, rep(0.1, 2*(nx-1)))  # initial value
#' fit = segCCR(data = ChIPseq,
#'            par.ini = par.ini,
#'            tm=tm,
#'            NB = 5)
#' }
#' 

'ChIPseq'


