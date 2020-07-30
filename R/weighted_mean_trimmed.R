#' Weighted trimmed mean and total (bare-bone functions)  
#'
#' Weighted trimmed mean and total (bare-bone functions with limited functionality; see \code{\link{svymean_trimmed}} and \code{\link{svytotal_trimmed}} for more capable methods)
#'
#' \describe{
#'    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote the estimated trimmed population mean; then, the estimated trimmed total is given by \eqn{N \mu} with \eqn{N =\sum w_i}, where summation is over all observations in the sample.}
#'    \item{Trimming.}{The methods trims the \code{LB}\eqn{~\cdot 100\%} percentage of the smallest observations and the (1 - \code{UB})\eqn{~\cdot 100\%} percentage of the largest observations from the data.} 
#'    \item{Variance estimation.}{See the related but more capable functions: 
#'    \itemize{
#'	 \item \code{\link{svymean_trimmed}}, 
#'	 \item \code{\link{svytotal_trimmed}}.
#'	 } 
#'    } 
#' } 
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param LB \code{[double]} lower bound of trimming, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param UB \code{[double]} upper bound of trimming, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param info \code{[logical]} indicating whether additional information should be returned (default: \code{FALSE}). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return The return value depends on \code{info}: \describe{
#'   \item{\code{info = FALSE}:}{estimate of mean or total \code{[double]}}
#'   \item{\code{info = TRUE}:}{a \code{[list]} with items: \itemize{
#'	 \item \code{characteristic} \code{[character]},
#'	 \item \code{estimator} \code{[character]},
#'	 \item \code{estimate} \code{[double]},
#'	 \item \code{variance} (default: \code{NA}),
#'	 \item \code{robust} \code{[list]},
#'	 \item \code{residuals} \code{[numeric vector]},
#'	 \item \code{model} \code{[list]},
#'	 \item \code{design} (default: \code{NA}),
#'	 \item \code{[call]}
#'	 }
#'    }
#' }
#' @examples
#' data(workplace) 
#'
#' # Estimated trimmed population total (5% symmetric trimming) 
#' weighted_total_trimmed(workplace$employment, workplace$weight, LB = 0.05, 
#'                        UB = 0.95) 
#' 
#' # Estimated trimmed population mean (5% trimming at the top of the distr.) 
#' weighted_mean_trimmed(workplace$employment, workplace$weight, UB = 0.95) 
#' @seealso \code{\link{svymean_trimmed}} and \code{\link{svytotal_trimmed}}
#' @export 
#' @useDynLib robsurvey wtrimmedmean
weighted_mean_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE, 
   na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   if (LB >= UB) stop("Argument 'LB' must be smaller than 'UB'!", call. = FALSE)
   if (LB < 0) stop("Argument 'LB' must not be < 0!", call. = FALSE)
   if (UB > 1) stop("Argument 'UB' must not be > 1!", call. = FALSE)
   tmp <- .C("wtrimmedmean", x = as.double(dat$x), w = as.double(dat$w),
      lb = as.double(LB), ub = as.double(UB), loc = as.double(numeric(1)),
      n = as.integer(dat$n), PACKAGE = "robsurvey")
   if (info){
      resid <- dat$x - tmp$loc
      res <- list(
	 characteristic = "mean", 
	 estimator = paste0("Weighted trimmed estimator (", LB, ", ", UB, ")"), 
	 estimate = tmp$loc, 
	 variance = NA,
	 robust = list(UB = UB, LB = LB),  
	 residuals = resid, 
	 model = list(y = dat$x, w = dat$w),
	 design = NA, 
	 call = match.call())
      return(res)
   }else{
      return(tmp$loc)
   }
}

#' @rdname weighted_mean_trimmed
#' @export 
weighted_total_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE, 
   na.rm = FALSE)
{
   res <- weighted_mean_trimmed(x, w, LB, UB, info, na.rm)
   if (length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
      res$call <- match.call()
   }
   return(res)
} 
