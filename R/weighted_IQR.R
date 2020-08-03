#' Weighted interquartile range (IQR) 
#'
#' \code{weighted_IQR} computes weighted (normalized) interquartile range 
#'
#' By default, the weighted IQR is normalized to be an unbiased estimate of scale at the Gaussian core model. If normalization is not wanted, put \code{constant = 1}.
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param constant \code{[double]} constant scaling factor to make the weighted IQR a consistent estimator of the scale (default: \code{0.7413}).
#' @return Weighted IQR 
#'
#' @examples
#' data(workplace) 
#'
#' # normalized weighted IQR (default) 
#' weighted_IQR(workplace$employment, workplace$weight)
#'
#' # weighted IQR (without normalization) 
#' weighted_IQR(workplace$employment, workplace$weight, constant = 1)
#' @export 
weighted_IQR <- function(x, w, na.rm = FALSE, constant = 0.7413)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   qs <- weighted_quantile(dat$x, dat$w, probs = c(0.25, 0.75), na.rm)
   unname((qs[2] - qs[1]) * constant)
} 
