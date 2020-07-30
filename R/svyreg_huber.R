#' Huber robust survey regression M-estimator 
#' 
#' Huber robust survey regression M-estimator
#' 
#' Details
#' 
#' @section Failure of convergence: 
#' By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}. 
#'
#' @param formula a \code{[formula]} object (i.e., symbolic description of the model)
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param var \code{[numeric vector]} heteroscedastic variance (default: \code{NULL}, i.e., homoscedastic variance).
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param ... additional arguments passed to the method (e.g., \code{maxit}: maxit number of iterations, etc.). 
#' @return object of class \code{svyreg.rob} 
#' @export 
svyreg_huber <- function(formula, design, var = NULL, k, na.rm = FALSE, ...)
{
#FIXME: better check function (check for NA and NaN and infinite)
   mf <- stats::model.frame(formula, design$variables, na.action = 
      stats::na.fail())
   n <- nrow(mf)
   # heteroscedasticity
   if (!is.null(var)){
      tmp <- stats::model.frame(var, design$variables, na.action = 
	 stats::na.fail())
      var <- tmp[[1]]
   } 
   res <- robsvyreg(x = stats::model.matrix(formula, mf), y = mf[[1]], 
      w = as.numeric(1 / design$prob), k = k,
      intercept = attr(attributes(mf)$terms, "intercept"), var = var, 
      na.rm = na.rm, ...)
   res$design <- design
   res$call <- match.call()
   res$model$yname <- names(mf)[[1]] 
   class(res) <- c("svyreg", "svyregrob")
   res
} 
