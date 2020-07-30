#' Weighted winsorized mean and total (bare-bone functions)  
#'
#' Weighted winsorized mean and total (bare-bone functions with limited functionality; see \code{\link{svymean_winsorized}} and \code{\link{svytotal_winsorized}} for more capable methods)
#'
#' \describe{
#'    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote the estimated winsorized population mean; then, the estimated winsorized total is given by \eqn{N \mu} with \eqn{N =\sum w_i}, where summation is over all observations in the sample.}
#'    \item{Modes of winsorization.}{The amount of winsorization can be specified in relative or absolute terms:
#'	 \itemize{
#'	    \item \emph{relative:} By specifying \code{LB} and \code{UB}, the methods winsorizes the \code{LB}\eqn{~\cdot 100\%} percentage of the smallest observations and the (1 - \code{UB})\eqn{~\cdot 100\%} percentage of the largest observations from the data. 
#'	    \item \emph{absolute:} By specifying argument \code{k} in the functions with the "infix" \code{_k_} in their name, the largest \eqn{k} observations are winsorized, \eqn{0<k<n}, where \eqn{n} denotes the sample size.
#'	 }
#'    } 
#'    \item{Variance estimation.}{See the related but more capable functions:
#'    \itemize{
#'	 \item \code{\link{svymean_winsorized}}, 
#'	 \item \code{\link{svytotal_winsorized}}, 
#'	 \item \code{\link{svymean_k_winsorized}}, 
#'	 \item \code{\link{svytotal_k_winsorized}}.
#'	 } 
#'    }
#' }
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param LB \code{[double]} lower bound of winsorization, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}. 
#' @param UB \code{[double]} upper bound of winsorization, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param k \code{[integer]} number of observations to be winsorized at the top of the distribution. 
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
#' # Estimated winsorized population mean (5% symmetric winsorization)
#' weighted_mean_winsorized(workplace$employment, workplace$weight, LB = 0.05)
#'
#' # Estimated one-sided k winsorized population total (2 observations are 
#' # winsorized at the top of the distribution)
#' weighted_total_k_winsorized(workplace$employment, workplace$weight, k = 2)
#' @seealso \code{\link{svymean_winsorized}}, \code{\link{svymean_k_winsorized}}, \code{\link{svytotal_winsorized}}, and \code{\link{svytotal_k_winsorized}}
#' @export 
#' @useDynLib robsurvey 
weighted_mean_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE, 
   na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   if (LB >= UB) stop("Argument 'LB' must be smaller than 'UB'!", call. = FALSE)
   if (LB < 0) stop("Argument 'LB' must not be < 0!", call. = FALSE)
   if (UB > 1) stop("Argument 'UB' must not be > 1!", call. = FALSE)
   tmp <- .C("wwinsorizedmean", x = as.double(dat$x), w = as.double(dat$w),
      lb = as.double(LB), ub = as.double(UB), loc = as.double(numeric(1)),
      n = as.integer(dat$n), PACKAGE = "robsurvey")
   if (info){
      resid <- dat$x - tmp$loc
      res <- list(
	 characteristic = "mean", 
	 estimator = paste0("Weighted winsorized estimator (", LB, ", ", UB, 
	    ")"), 
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

#' @rdname weighted_mean_winsorized
#' @export 
#' @useDynLib robsurvey 
weighted_mean_k_winsorized <- function(x, w, k, info = FALSE, na.rm = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   if (k %% 1 > 0){
      k <- as.integer(k)
      cat(paste0("Argument 'k' is casted to integer: k = ", k,"\n"))
   }
   n <- dat$n
   if (k >= n) stop("k must be smaller than n\n", call. = FALSE)
   if (k < 1) stop("k must larger than 1\n", call. = FALSE)
   tmp <- .C("wkwinsorizedmean", x = as.double(dat$x), w = as.double(dat$w),
      k = as.integer(k - 1), loc = as.double(numeric(1)), n = as.integer(n), 
      prob = as.double(numeric(1)), PACKAGE = "robsurvey")
   if (info){
      res <- list(
	 characteristic = "mean",
	 estimator = paste0("weighted k winsorized estimator (k = ", k, ")"),
	 estimate = tmp$loc,
	 variance = NA,
	 robust = list(k = k, UB = tmp$prob), 
	 residuals = dat$x - tmp$loc,
	 model = list(y = x, w = w),
	 design = NA,
	 call = match.call())
      return(res)
   }else{
      return(tmp$loc)
   }
}

#' @rdname weighted_mean_winsorized
#' @export 
weighted_total_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB, 
   info = FALSE, na.rm = FALSE)
{
   res <- weighted_mean_winsorized(x, w, LB, UB, info, na.rm)
   if (length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
      res$call <- match.call()
   }
   return(res)
}

#' @rdname weighted_mean_winsorized
#' @export 
weighted_total_k_winsorized <- function(x, w, k, info = FALSE, na.rm = FALSE)
{ 
   res <- weighted_mean_k_winsorized(x, w, k, info, na.rm)
   if (length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
      res$call <- match.call()
   }
   return(res)
} 
