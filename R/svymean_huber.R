#' Weighted Huber mean and total - Robust Horvitz-Thompson estimator  
#'
#' Weighted Huber M-estimator of the population mean and total (robust Horvitz-Thompson estimator) 
#'
#' \describe{
#'    \item{Overview.}{}
#'    \item{Methods.}{\code{type = "rht"} anad \code{type = "rwm"}}
#'    \item{Variance estimation.}{Taylor linearization (residual variance estimator).}
#'    \item{Utility functions.}{\code{\link[=svystat.rob]{summary}}, \code{\link[=svystat.rob]{coef}}, \code{\link[=svystat.rob]{SE}}, \code{\link[=svystat.rob]{vcov}}, \code{\link[=svystat.rob]{residuals}}, \code{\link[=svystat.rob]{fitted}}, and \code{\link[=svystat.rob]{robweights}}.}
#'    \item{Bare-bone functions.}{See \code{\link{weighted_mean_huber}} and \code{\link{weighted_total_huber}}.}
#' } 
#'
#' @section Failure of convergence: 
#' By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}. 
#'
#' @param x a one-sided \code{[formula]}, e.g., \code{~myVariable}. 
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}; default: \code{k = 1.5}). 
#' @param type \code{[character]} type of method: \code{"rwm"} or \code{"rht"}. 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param ... additional arguments passed to the method (e.g., \code{maxit}: maxit number of iterations, etc.). 
#' @return object of class \code{\link{svystat.rob}} 
#' @seealso \code{\link{weighted_mean_huber}} and \code{\link{weighted_total_huber}}
#' @references Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators, \emph{Surv. Methodol.} 21, pp. 79-87.
#' @examples
#' data(workplace) 
#'
#' library(survey)
#' # Survey design for simple random sampling without replacement
#' dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight, 
#'                 data = workplace)
#'
#' # Robust Horvitz-Thompson M-estimator of the population total 
#' svytotal_huber(~employment, dn, k = 9, type = "rht")
#' 
#' # Robust weighted M-estimator of the population mean 
#' svymean_huber(~employment, dn, k = 12, type = "rwm")
#' @export
svymean_huber <- function(x, design, k = 1.5, type = "rwm", na.rm = FALSE, ...)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_huber(dat$x, dat$w, k, type, info = TRUE, na.rm, ...)
   # modify residuals for type 'rht' (only for variance estimation)
   if (type == "rht"){
      r <- sqrt(res$model$var) * res$model$y - res$estimate 
   } else {
      r <- res$residuals
   }
   # compute variance

# FIXME: take y and from res
   infl <- res$robust$robweights * r * dat$w / sum(dat$w) 
   res$variance <- survey::svyrecvar(infl, design$cluster, design$strata, 
        design$fpc, postStrata = design$postStrata)
   names(res$estimate) <- dat$xname
   res$call <- match.call()
   res$design <- design
   class(res) <- "svystat.rob"
   res
}

#' @rdname svymean_huber 
#' @export 
svytotal_huber <- function(x, design, k = 1.5, type = "rwm", na.rm = FALSE, ...)
{
   res <- svymean_huber(x, design, k, type, na.rm, ...)  
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
}

# #' @describeIn svymean_huber 
# #' @export
# fixed_downweighting <- function(object, at = 0.95){
#    if (class(object) != "svystat.rob") stop("method not supported for class: ", 
#       class(object), "\n") 
#    if (object$call[[1]] == "svymean_huber" && 
#       object$call[[1]] == "svytotal_huber"){
#       stop("method: '", object$call[[1]], "' is not supported\n") 
#    }
#    stopifnot(at > 0.5, at < 1)
#    ctrl <- svyreg_control()
#
#    cl <- object$call 
#    foo <- function(k, at){
#       cl$k <- k
#       tmp <- eval(cl)
#       mean(tmp$robust$robweights) - at
#    }
#    res <- uniroot(foo, interval = c(0.0001, ctrl$k_Inf), at = at)
#    res$root
# }
# 
