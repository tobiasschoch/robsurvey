#' Utility functions 
#'
#' Methods and utility functions for objects of classs code{svyreg_rob} 
#'
#' Utility functions
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
#' @param object, x object of class \code{svyreg_rob}.
#' @param x object of class \code{svyreg_rob}.
#' @param digits \code{[integer]} minimal number of significant digits.
#' @param which indicating which plots to be drawn; if a subset of the plots is required, you can specify a subset of the numbers \code{1:5}. 
#' @param ... additional arguments passed to the method. 
#' @name class_svyreg_rob
#' @aliases svyreg_rob 

NULL

#' @rdname class_svyreg_rob
#' @method print svyreg_rob 
print.svyreg_rob <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
   converged <- x$optim$converged 
   if (is.null(converged) || converged) { 
      cat(paste0("\n", x$estimator, "\n"))
      cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	 "\n\n", sep = "")
      if (length(x$estimate)) {
	 cat("Coefficients:\n")
	 print.default(format(x$estimate, digits = digits), print.gap = 2L, 
	    quote = FALSE)
      }
   } else { 
      cat(paste0("\nIRWLS not converged in ", x$optim$niter, 
	 " iterations (with tol = ", x$optim$tol, ")\n"))
   }
   cat("\n")
   invisible(x)
} 

#' @rdname class_svyreg_rob
#' @method summary svyreg_rob 
summary.svyreg_rob <- function(object, digits = max(3L, getOption("digits") 
   - 3L), ...)
{
   converged <- object$optim$converged 
   if (is.null(converged) || converged) { 
      w <- object$model$w
      r <- object$residuals
      yhat <- object$model$y - r 
      n <- object$model$n; p <- object$model$p
      intercept <- object$model$intercept

      cat("\nCall:\n")
      print(object$call)
      cat("\nResiduals:\n")
      print(summary(r))
      cat("\nCoefficients:\n")

   #FIXME:
      # residual variance
      res_var <- object$variance * n / (n - p)

      # inverse of X^T %*% X
      xtx_inv <- chol2inv(qr.R(qr(object$model$x)))
      se <- sqrt(diag(xtx_inv) * res_var)
      est <- object$estimate
      tval <- est / se
      stats::printCoefmat(cbind(Estimate = est, 'Std. Error' = se, 
	 't value' = tval, 'Pr(>|t|)' = 2 * stats::pt(abs(tval), n - p, 
	 lower.tail = FALSE)), digits = digits)
      cat("---\n")
      cat("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n\n") 
   #FIXME:
   #   cat("Residual standard error: 150 on 198 degrees of freedom\n")

      # coefficient of determination
      rss <- sum(w * r^2)
      if (intercept) {
	 m <- sum(w * yhat) / sum(w)
	 mss <- sum(w * (yhat - m)^2) 
      } else
	 mss <- sum(w * yhat^2)

      R2 <- mss / (mss + rss)   

#FIXME: N in place of n?
      R2_adj <- 1 - (1 - R2) * (n - intercept) / (n - p)
      cat("Multiple R-squared: ", round(R2, digits), "   Adjusted R-squared: ", 
	 round(R2_adj, digits), "\n\n")

      if (!is.null(object$robust)) {
	 cat("\nRobustness:\n")
      }
   } else {
      cat(paste0("\nIRWLS not converged in ", object$optim$niter, 
	 " iterations (with tol = ", object$optim$tol, ")\n"))
   }
   cat("\n")
}

#' @rdname class_svyreg_rob 
#' @method coef svyreg_rob
coef.svyreg_rob <- function(object, ...)
{
   object$estimate
}

#' @rdname class_svyreg_rob
#' @method vcov svyreg_rob
vcov.svyreg_rob <- function(object, ...)
{
#FIXME:
   v <- as.matrix(object$variance)
   rownames(v) <- names(object$estimate)
   colnames(v) <- "Variance"
   v
}

#' @rdname class_svyreg_rob
#' @method residual svyreg_rob
residuals.svyreg_rob <- function(object, ...)
{
   object$residuals
}

#' @rdname class_svyreg_rob
#' @method fitted svyreg_rob
fitted.svyreg_rob <- function(object, ...)
{
   object$model$y - object$residuals
}

#' @rdname class_svyreg_rob
#' @method robweights svyreg_rob 
robweights.svyreg_rob <- function(object)
{
   tmp <- object$robust$robweights
   if (is.null(tmp)){ 
      warning("Robustness weights are not available\n")
      NA
   } else {
      tmp
   }
}

#' @rdname class_svyreg_rob
#' @method plot svyreg_rob 
#' @importFrom graphics plot 
plot.svyreg_rob <- function (x, which = 1:5, ...) 
{
}


# b <- lm(Petal.Width ~ Sepal.Width + Sepal.Length, iris, x = TRUE, y = TRUE)
# m <- list(coefficients = b$coefficients,
#    model = list(x = b$x, y = b$y, w = rep(1, 150), var = rep(1, 150)),
#    robust = list(robweights = rep(1, 150), scale = 1),
#    fitted.values = fitted.values(b), 
#    residuals = residuals(b),
#    qr = b$qr,
#    call = "a"
# )
#
# plot.svyreg(m)
#
#
