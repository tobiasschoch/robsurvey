#' Weighted median
#'
#' \code{weighted_median} computes the weighted sample quantile
#'
#' Weighted sample median; see \code{\link{weighted_quantile}} for more information. 
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return Weighted estimate of the population median 
#'
#' @examples
#' data(workplace) 
#'
#' weighted_median(workplace$employment, workplace$weight) 
#' @seealso \code{\link{weighted_quantile}}
#' @export 
weighted_median <- function(x, w, na.rm = FALSE)
{
   weighted_quantile(x, w, 0.5, na.rm)
} 
