#' Dalen estimator of the mean and total 
#'
#' Dalén's estimator of the mean and total. 
#'
#' Weight reduction  
#'
#' @param x \code{[numeric vector]} observations.
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param censored \code{[double]} threshold above which all observations are censored (i.e., set equal to the threshold). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param verbose \code{[logical]} indicating whether additional information should be \emph{printed} to the console (default: \code{FALSE}).
#' @param info \code{[logical]} indicating whether additional information should be returned (default: \code{FALSE}). 
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
#'
#' @examples
#' data(workplace) 
#'
#' # Dalen's estimator of the total (with censoring threshold: 100000)  
#' weighted_total_dalen(workplace$employment, workplace$weight, 100000)
#' @references Dalén, J. (1987). Practical Estimators of a Population Total Which Reduce the Impact of Large Observations, Research Report, Statistics Sweden.
#' @export 
weighted_mean_dalen <- function(x, w, censored, na.rm = FALSE, verbose = TRUE, 
   info = FALSE)
{ 
   res <- robsurvey::weighted_total_dalen(x, w, censored, na.rm, verbose, 
      info = TRUE)
   res$characteristic <- "mean"
   res$estimate <- res$estimate / sum(res$robust$weightsmod)
   res$call <- match.call()
   if (info){
      return(res)
   }else{
      return(res$estimate)
   }
}

#' @rdname weighted_mean_dalen
#' @export 
weighted_total_dalen <- function(x, w, censored, na.rm = FALSE, verbose = TRUE, 
   info = FALSE)
{ 
   stopifnot(censored > 0)
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   xw <- dat$x * dat$w
   if (verbose){
      cat(paste0(sum(xw > censored), " of ", length(x), 
	 " observations censored\n"))
   }
   at <- xw > censored
   if (sum(at) > 0){
      xw[at] <- censored + (xw[at] - censored) / dat$w[at]
      weightsmod <- xw / dat$x
   }
   if (info){
      res <- list(
	 characteristic <- "total",
	 estimator = paste0("Dalen estimator (censored at ", censored, ")"),
	 estimate = sum(xw),
	 variance = NA,
	 robust = list(censored = censored, weightsmod = weightsmod),
	 residuals = NA,
	 model = list(y = dat$x, w = dat$w),
	 design = NA,
	 call = match.call()
      )
   }else{
      res <- sum(xw)
   }
   return(res)
} 
