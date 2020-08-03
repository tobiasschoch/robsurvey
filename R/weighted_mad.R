#' Weighted median absolute deviation from the median (MAD)
#'
#' \code{weighted_mad} computes weighted median of the absolute deviations from the weighted median
#'
#' The weighted MAD is computed as the (normalized) weighted median of the absolute deviation from the weighted median; see \code{\link{weighted_median}}. The weighted MAD is normalized to be an unbiased estimate of scale at the Gaussian core model. If normalization is not wanted, put \code{constant = 1}.
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param constant \code{[double]} constant scaling factor to make the MAD a consistent estimator of the scale (default: \code{1.4826}).
#' @return Weighted median absolute deviation from the (weighted) median
#'
#' @examples
#' data(workplace) 
#'
#' # normalized weighted MAD (default) 
#' weighted_mad(workplace$employment, workplace$weight)
#'
#' # weighted MAD (without normalization) 
#' weighted_mad(workplace$employment, workplace$weight, constant = 1)
#' @export 
weighted_mad <- function(x, w, na.rm = FALSE, constant = 1.4826)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   med <- weighted_quantile(dat$x, dat$w, probs = 0.5, na.rm)
   mad <- weighted_quantile(abs(dat$x - med), dat$w, probs = 0.5, na.rm)
   unname(mad * constant)
} 
