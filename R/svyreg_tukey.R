#' Huber robust survey regression M- and GM-estimator 
#' 
#' Huber robust survey regression M- and GM-estimator (Mallows and Schweppe type) 
#' 
#' Details
#' 
#' @section Failure of convergence: 
#' By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}. 
#'
#' @param formula a \code{[formula]} object (i.e., symbolic description of the model)
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}).
#' @param type \code{[character]} \code{"Mallows"} or \code{"Schweppe"}. 
#' @param xwgt \code{[numerical vector]} of weights in the design space (default: \code{NULL}); \code{xwgt} is only relevant if \code{type = "Mallows"} or \code{type = "Schweppe"}. 
#' @param var \code{[numeric vector]} heteroscedastic variance (default: \code{NULL}, i.e., homoscedastic variance).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param ... additional arguments passed to the method (e.g., \code{maxit}: maxit number of iterations, etc.). 
#' @return object of class \code{svyreg.rob} 
#' @export 
svyreg_tukey <- function(formula, design, k, var = NULL, na.rm = FALSE, 
   ...)
{
   dat <- .checkreg(formula, design, var, na.rm)
   res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, 0, NULL, dat$var, ...)
   res$design <- design
   res$call <- match.call()
   res$model$intercept <- dat$intercept
   res$model$yname <- dat$yname 
   class(res) <- "svyreg_rob"
   res
}

#' @rdname svyreg_tukey
#' @export 
svyreg_tukeyGM <- function(formula, design, k, type, xwgt, var = NULL, 
   na.rm = FALSE, ...)
{
   type_int <- switch(toupper(type), "MALLOWS" = 1, "SCHWEPPE" = 2)
   if (is.null(type_int))
      stop("Type '", type,"' is not defined\n")

   dat <- .checkreg(formula, design, var, na.rm)

   if (NCOL(xwgt) > 1) {
      xwgt <- as.numeric(xwgt[, 1])
      warning("Only first column of argument 'xwgt' is used\n")
   }
   if (length(xwgt) != length(dat$y))
      stop("Argument 'xwgt' is not of length n\n")

   res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, type_int, xwgt, dat$var, ...)
   res$design <- design
   res$call <- match.call()
   res$model$intercept <- dat$intercept
   res$model$yname <- dat$yname 
   res$model$xwgt <- xwgt 
   class(res) <- "svyreg_rob"
   res
} 
