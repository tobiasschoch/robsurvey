#' Utility functions 
#'
#' Methods and utility functions for objects of class \code{svyreg.rob} 
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
#' @param object, x object of class \code{svyreg.rob}.
#' @param x object of class \code{svyreg.rob}.
#' @param digits \code{[integer]} minimal number of significant digits.
#' @param which indicating which plots to be drawn; if a subset of the plots is required, you can specify a subset of the numbers \code{1:5}. 
#' @param ... additional arguments passed to the method. 
#' @name class_svyreg.rob

NULL

#' @rdname class_svyreg.rob 
#' @method summary svyreg 
summary.svyreg <- function(object, digits = max(3L, getOption("digits") - 3L), 
   ...)
{
   cat("\nCall:\n")
   print(object$call)
   cat("\nResiduals:\n")
   print(summary(object$residuals))
   cat("\nCoefficients:\n")
   n <- object$model$n; p <- object$model$p
   # residual variance
   res_var <- object$variance * n / (n - p)
#   res_var <- sum(object$model$w * object$residual^2) / (n-p)

   if (!is.null(object$robust)){ 
      res_var <- res_var * object$robust$Epsi2 / object$robust$Epsiprime^2 
   }
   # inverse of X^T %*% X
   xtx_inv <- chol2inv(qr.R(qr(object$model$x)))
   se <- sqrt(diag(xtx_inv) * res_var)
   est <- object$estimate
   tval <- est / se
   stats::printCoefmat(cbind(Estimate = est, 'Std. Error' = se, 
      't value' = tval, 'Pr(>|t|)' = 2 * stats::pt(abs(tval), n - p, 
      lower.tail = FALSE)), digits = digits)
   y <- object$model$y 
   if (object$model$intercept){
      r2 <- 1 - sum(object$residuals^2) / sum((y - mean(y))^2)
   }else{
      r2 <- 1 - sum(object$residuals^2) / sum(y^2)
   }
   cat("Multiple R-squared: ", round(r2, digits), "   Adjusted R-squared: ", 
      round(1 - (1 - r2) * (n - 1) / (n - p - 1), digits), "\n\n")
}

#' @rdname class_svyreg.rob 
#' @method plot svyreg
#' @importFrom graphics par mtext plot symbols abline grid panel.smooth text 
plot.svyreg <- function(x, which = 1:5)
{
   fitted <- x$model$y - x$residuals
   show <- rep(FALSE, 5)
   show[which] <- TRUE
   par(ask = TRUE)
   n <- length(x$model$y)
   v <- sqrt(x$variance)
   stdres <- x$residuals / v
   if (show[1L]){
      r <- range(c(x$model$y, fitted))
      plot(x$model$y, fitted, xlab = x$model$yname, ylab = "fitted", 
	 main = "Response variable vs. fitted", type = "n", ylim = r,
	 xlim = r)
      mtext(text = "(circe size is proportional to sampling weights)", 
	 side = 3, line = 0.5)
      scale <- (diff(range(x$model$y)) / n) / mean(x$model$w) 
      bg <- rep(NA, n)
      if (!is.null(x$robust)){
	 bg[x$robust$robweights < 1] <- "black"
      }
      symbols(x$model$y, fitted, circles = scale * x$model$w, inches = FALSE, 
	 add = TRUE, bg = bg)
      abline(0, 1)
      grid(col = "grey65")
   }
   if (show[2L]){
      plot(fitted, stdres, main = "Standarized Residuals vs. Fitted", 
	 xlab = "Fitted value", ylab = "Standarized residual", type = "n")
      panel.smooth(fitted, stdres)
      grid(col = "grey65")
      abline(h = 0, lty = 2, col = "gray65")
      at <- which(abs(stdres) > 1.96)
      if (length(at) > 0){
	 atboundary <- min(fitted) == fitted[at]
	 if (sum(atboundary) > 0 & sum(atboundary) < length(atboundary)){
	    text(fitted[at[atboundary]], stdres[at[atboundary]], 
	       labels = at[atboundary], cex = 0.8, pos = 4)
	    text(fitted[at[!atboundary]], stdres[at[!atboundary]], 
	       labels = at[!atboundary], cex = 0.8, pos = 2)
	 }else{
	    text(fitted[at], stdres[at], labels = at, cex = 0.8, pos = 2)
	 }
      }
   }
   if (show[3L]){
      u <- stats::qqnorm(stdres, main = "Normal QQ-Plot", ylab = "Standardized 
	 residual", xlab = "Theoretical quantile")
      stats::qqline(stdres, lty = 2, col = "gray65") 
      grid(col = "grey65")
      at <- which(u$y < -1.96)
      if (length(at) > 0){
	 text(u$x[at], u$y[at], labels = at, cex = 0.8, pos = 4)
      }
      at <- which(u$y > 1.96)
      if (length(at) > 0){
	 text(u$x[at], u$y[at], labels = at, cex = 0.8, pos = 2)
      }
   }   
   if (show[4L]){
      hii <- rowSums(qr.Q(qr(x$model$x / v))^2)
      plot(hii, stdres, main = "Standardized Residuals vs. Leverage", 
	 xlab = "Leverage", ylab = "Standardized residual", 
	 xlim = c(0, max(hii)), type = "n")
      panel.smooth(hii, stdres)
      grid(col = "grey65")
      at <- which(abs(stdres) > 1.96)
      if (length(at) > 0){
	 text(hii[at], stdres[at], labels = at, cex = 0.8, pos = 2)
      }
      abline(h = 0, lty = 2, col = "gray65")
      abline(v = 0, lty = 2, col = "gray65")
   }
   if (show[5L] & exists("wBACON")){
#FIXME:
#       xdat <- x$model$x[, (1 + x$model$intercept):model$p]
# #      tmp <- wBACON(xdat, x$model$Di, verbose = FALSE)
#       tmp <- 1
#
#       plot(tmp$dist, stdres, xlab = "Robust distance", type = "n",
# 	 ylab = "Standardized residual", xlim = c(0, 1.05 * max(tmp$dist)), 
# 	 xaxs = "i")
#       outer <- par("usr")
#       rect(outer[1], outer[3], outer[2], outer[4], col = "grey90", 
# 	 border = NA)
#       inner <- c(0, -1.96, tmp$cutoff, 1.96)
#       rect(inner[1], inner[2], inner[3], inner[4], border = NA, 
# 	 col = "white")
#       box()
#       grid(col = "grey65")
#       points(tmp$dist, stdres)
#       at <- c(which(abs(stdres) > 1.96), which(tmp$dist > tmp$cutoff))
#       flag <- at[duplicated(at)]
#       if (!is.null(flag) > 0){
# 	 text(tmp$dist[flag], stdres[flag], labels = flag, cex = 0.8, 
# 	    pos = 2)
#       }
#       at <- unique(at)
#       if (!is.null(at)){
# 	 points(tmp$dist[at], stdres[at])
#       }
#       abline(h = 0, lty = 2)
#       abline(v = tmp$cutoff, lty = 3)
#       abline(h = c(-1.96, 1.96), lty = 3)
   }
}

#' @rdname class_svyreg.rob 
#' @method print svyreg
print.svyreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
   cat(x$estimator, "\n\nCall:\n")
   print(x$call)
   cat("\nCoefficients:\n")
   print.default(format(coef.svystat.rob(x), digits = digits), print.gap = 2L, 
      quote = FALSE)
} 
