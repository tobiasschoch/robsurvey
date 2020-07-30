#' Weighted Huber Proposal 2 estimator 
#'
#' Weighted Huber Proposal 2 estimator 
#'
#' The function \code{huber2} computes the weighted Huber (1964) Proposal 2 estimates of location and scale. 
#' 
#' The method is initialized by the weighted median (location) and the weighted interquartile range (scale). 
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param maxit \code{[integer]} maximum number of iterations to use (default: \code{50}). 
#' @param tol \code{[double]} numerical tolerance criterion to stop the iterations (default: \code{1e-04}). 
#' @param info \code{[logical]} observations (default: \code{FALSE}).
#' @return The return value depends on \code{info}: \describe{
#'   \item{\code{info = FALSE}:}{estimate of mean or total \code{[double]}}
#'   \item{\code{info = TRUE}:}{a \code{[list]} with items: 
#'	 \itemize{
#'	    \item \code{characteristic} \code{[character]},
#'	    \item \code{estimator} \code{[character]},
#'	    \item \code{estimate} \code{[double]},
#'	    \item \code{variance} (default: \code{NA}),
#'	    \item \code{robust} \code{[list]},
#'	    \item \code{residuals} \code{[numeric vector]},
#'	    \item \code{model} \code{[list]},
#'	    \item \code{design} (default: \code{NA}),
#'	    \item \code{[call]}
#'	 }
#'    }
#' }
#' @examples
#' data(workplace) 
#'
#' # Weighted Proposal 2 M-estimator of the mean
#' huber2(workplace$employment, workplace$weight, k = 8) 
#' @references Huber, P. J. (1964). Robust Estimation of a Location Parameter, \emph{Ann. Math. Statist.} 35, pp. 73--101.   
#' @export 
#' @useDynLib robsurvey huberm 
huber2 <- function(x, w, k, na.rm = FALSE, maxit = 50, tol = 1e-4, info = FALSE)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   if (is.finite(k)) 
      kk <- k
   else     
      kk <- 1e5 

   tmp <- .C("huberm", x = as.double(x), w = as.double(w), 
      robwgt = as.double(numeric(dat$n)), k = as.double(kk), 
      loc = as.double(numeric(1)), scale = as.double(numeric(1)), 
      n = as.integer(dat$n), maxit = as.integer(maxit), tol = as.double(tol), 
      PACKAGE = "robsurvey")
   if (info){
      resid <- dat$x - tmp$loc
      res <- list(
	 characteristic = "mean", 
	 estimator = paste0("Weighted Huber proposal 2 estimator (k=", k,")"), 
	 estimate = tmp$loc, 
	 variance = NA,
	 robust = list(k = k, robweights = tmp$robwgt),  
	 residuals = tmp$resid, 
	 model = list(y = dat$x, w = dat$w),
	 design = NA, 
	 call = match.call())
      return(res)
   }else{
      return(tmp$loc)
   }
} 
