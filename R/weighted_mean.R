#' Weighted total and mean (Horvitz-Thompson and Hajek estimators)
#'
#' Weighted total and mean (Horvitz-Thompson and Hajek estimators)
#'
#' Computation of the Horvitz-Thompson and the Hajek estimator of, 
#' respectively, the total and the mean
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return estimated population mean or total 
#'
#' @examples
#' data(workplace) 
#'
#' # Horvitz-Thompson estimator of the total
#' weighted_total(workplace$employment, workplace$weight)
#'
#' # Hajek estimator of the mean
#' weighted_mean(workplace$employment, workplace$weight)
#' @export 
weighted_mean <- function(x, w, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   return(sum(dat$x * dat$w) / sum(dat$w))
} 

#' @rdname weighted_mean
#' @export 
weighted_total <- function(x, w, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   return(sum(dat$w * dat$x))
} 
