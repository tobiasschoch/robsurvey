#' Utility functions 
#'
#' Methods and utility functions for objects of class \code{svystat.rob} 
#'
#' Utility functions:
#' \itemize{
#'    \item \code{summary} gives a summary of the estimation properties
#'    \item \code{robweights} extracts the robustness weights
#'    \item \code{coef} extracts the estimates 
#'    \item \code{SE} extracts the (estimated) standard error
#'    \item \code{vcov} extracts the (estimated) covariance matrix
#'    \item \code{residuals} extracts the residuals
#'    \item \code{fitted} extracts the fitted values
#' }
#'
#' @param object object of class \code{svystat.rob}.
#' @param x object of class \code{svystat.rob}.
#' @param digits \code{[integer]} minimal number of significant digits.
#' @param ... additional arguments passed to the method. 
#' @name class_svystat.rob 
#' @aliases svystat.rob 

NULL

#' @rdname class_svystat.rob 
#' @export 
summary.svystat.rob <- function(object, digits = max(3L, getOption("digits") - 
   3L), ...)
{
   cat(paste0(object$estimator, " of the sample ",
      object$characteristic, "\n"))
   cat("\n")
   est <- cbind(object$estimate, sqrt(object$variance))
   colnames(est) <- c(object$characteristic, "SE")
   print(est, digits)
   cat("\n")
   if(!is.null(object$optim)){
      cat("Robustness:\n")
      cat(paste0("  Psi-function: ", object$robust$psifunction, " with k = ",
	 object$robust$k, "\n"))
      cat(paste0("  mean of robustness weights: ", 
	 round(mean(object$robust$robweights), digits), "\n"))
      cat("\n")
      cat("Algorithm performance:\n")
      if (object$optim$converged){
	 cat(paste0("  converged in ", object$optim$niter, " iterations \n"))
	 cat(paste0("  with residual scale (weighted MAD): ", 
	    format(object$robust$scale, digits = digits), "\n"))
      }else{
	 cat(paste0("  FAILURE of convergence in ", object$optim$niter, 
	    " iterations \n"))
	 cat(paste0("  with residual scale (weighted MAD): ", 
	    round(object$robust$scale, digits), "\n"))
      }
   }
#FIXME:
   #cat("Sampling design:\n")
   #print(object$design)
}

#' @rdname class_svystat.rob 
#' @export 
coef.svystat.rob <- function(object, ...)
{
   object$estimate
}

#' @rdname class_svystat.rob 
#' @export 
SE.svystat.rob <- function(object, ...)
{
   sqrt(object$variance)
}

#' @rdname class_svystat.rob 
#' @export 
vcov.svystat.rob <- function(object, ...)
{
   v <- as.matrix(object$variance)
   rownames(v) <- names(object$estimate)
   colnames(v) <- "Variance"
   v
}

#' @rdname class_svystat.rob 
#' @export 
residuals.svystat.rob <- function(object, ...)
{
   object$residuals
}

#' @rdname class_svystat.rob 
#' @export 
fitted.svystat.rob <- function(object, ...)
{
   object$model$y - object$residuals
}

#' @rdname class_svystat.rob 
#' @export 
robweights <- function(object)
{
   UseMethod("robweights", object)
}

#' @rdname class_svystat.rob 
#' @export 
robweights.svystat.rob <- function(object)
{
   tmp <- object$robust$robweights
   if (is.null(tmp)){ 
      stop("Robustness weights are not available\n")
   } else {
      tmp
   }
}

#' @rdname class_svystat.rob 
#' @export 
print.svystat.rob <- function(x, digits = max(3L, getOption("digits") - 3L), 
   ...)
{
   conv <- TRUE
   if(!is.null(x$optim)){
      conv <- x$optim$converged
   }
   if(conv){
      m <- cbind(x$estimate, sqrt(x$variance))
      colnames(m) <- c(x$characteristic, "SE")
      print(m, digits)
   }else{
      cat(paste0(x$call[[1]], ": failure of convergence in ", x$optim$niter,
	 " steps\n"))
      cat("(you may use the 'summary' method to see more details)\n")
   }
} 
