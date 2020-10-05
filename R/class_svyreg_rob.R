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

coef.svyreg_rob <- function(object, ...)
{
   object$estimate
}

vcov.svyreg_rob <- function(object, ...)
{
#FIXME:
   v <- as.matrix(object$variance)
   rownames(v) <- names(object$estimate)
   colnames(v) <- "Variance"
   v
}

residuals.svyreg_rob <- function(object, ...)
{
   object$residuals
}

fitted.svyreg_rob <- function(object, ...)
{
   object$model$y - object$residuals
}

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
