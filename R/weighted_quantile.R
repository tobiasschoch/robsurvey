#' Weighted sample quantiles
#'
#' \code{weighted_quantile} computes the weighted sample quantiles
#'
#' \describe{
#'     \item{Overview.}{\code{weighted_quantile} computes the weighted sample quantiles; argument \code{probs} allows vector inputs.}
#'     \item{Implementation.}{The function is based on a weighted version of the quickselect algorithm with the Bentley and McIlroy (1993) 3-way partitioning scheme. For very small arrays, we use insertion sort.}
#'    \item{Compatibility.}{For equal weighting, i.e. when all elements in \code{w} are equal, \code{weighted_quantile} computes quantiles that are identical with \code{type = 2} in \code{stats::quantile}; see also Hyndman and Fan (1996).}
#' }
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param probs \code{[numeric vector]} vector of probabilities with values in \code{[0,1]}.
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return Weighted estimate of the population quantiles
#'
#' @examples
#' data(workplace) 
#'
#' # Weighted 25% quantile
#' weighted_quantile(workplace$employment, workplace$weight, 0.25) 
#' @references  
#' Bentley, J.L. and D.M. McIlroy (1993). Engineering a Sort Function, \emph{Software - Practice and Experience} 23, pp. 1249-1265.
#'
#' Hyndman, R.J. and Y. Fan (1996). Sample Quantiles in Statistical Packages, \emph{The American Statistician} 50, pp. 361-365.
#' @seealso \code{\link{weighted_median}}
#' @export 
#' @useDynLib robsurvey wquantile
weighted_quantile <- function(x, w, probs, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   if (any(probs < 0) | any(probs > 1)){
      stop("Argument 'probs' not in [0, 1]\n", call. = FALSE)
   }
   res <- NULL
   for (i in 1:length(probs)){
      tmp <- .C("wquantile", x = as.double(dat$x), w = as.double(dat$w), 
	 n = as.integer(dat$n), probs = as.double(probs[i]), 
	 q = as.double(numeric(1)), PACKAGE = "robsurvey") 
      res <- c(res, tmp$q)
   }
   names(res) <- paste0(probs * 100, "%")
   return(res)
} 
