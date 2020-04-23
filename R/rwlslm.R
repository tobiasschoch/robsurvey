#==============================================================================
# robust regression 
# ----------------------------------------------------------------------------- 
# PROJECT  sctbase 							     
# SUBEJCT  R function and test cases  
# AUTHORS  Tobias Schoch (tobias.schoch@fhnw.ch), February 10, 2020	     
# LICENSE  GPL >= 2							     
# COMMENT  [none]							     
#==============================================================================
setwd("C:/My/C/sctbase")
dyn.load("sctbase.dll")

#' Huber robust survey regression M-estimator 
#' 
#' Huber regression
#'  
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted.  
#' @param design survey design (object of class \code{svydesign})
#' @param k robustness tuning constant (Huber psi-function) 
#' @param ... further arguments (see \code{svyreg.control})  
#' @return object of class \code{svyreg} and \code{svyregrob} 
svyreg.huber.survey.design <- function(formula, design, k, ...){
   stopifnot(k > 0)
   call <- match.call()
   call[[1]] <- substitute(svyreg.huber)
   ctrl <- svyreg.control(...)
   mf <- model.frame(formula, design$variables, na.action = na.fail)
   xmat <- model.matrix(formula, mf)
   intercept <- attr(attributes(mf)$terms, "intercept")
   n <- nrow(xmat); p <- ncol(xmat)
   y <- mf[[1]]
   w <- as.numeric(weights(design))

   tmp <- .C("rwlslm", x = as.double(xmat), y = as.double(y), 
      w = as.double(w), resid = as.double(numeric(n)), 
      robwgt = as.double(numeric(n)), n = as.integer(n), p = as.integer(p), 
      k = as.double(k), beta = as.double(numeric(p)), 
      scale = as.double(numeric(1)), maxit = as.integer(ctrl$maxit), 
      tol = as.double(ctrl$tol), psitype = as.integer(ctrl$psi), 
      Epsi2 = as.double(numeric(1)), Epsiprime = as.double(numeric(1)))

   psi_fun <- ifelse(ctrl$psi == 0, "Huber", "asymHuber")
   names(tmp$beta) <- colnames(xmat) 
   design$variables <- NULL # drop the variables in the svydesign object
   res <- list(characteristic = "regression", 
      estimator = paste0("Survey regression M-estimator (", psi_fun,", k = ", 
	 k, ")"), 
      estimate = tmp$beta, 
      variance = tmp$scale^2,
      robust = list(psifunction = ifelse(ctrl$psi == 0, "Huber", "asymHuber"), 
	 k = k, robweights = tmp$robwgt, outliers = sum(tmp$robwgt < 1), 
	 Epsi2 = tmp$Epsi2, Epsiprime = tmp$Epsiprime), 
      optim = list(converged = (tmp$maxit != 0), niter = ifelse(tmp$maxit == 0, 
	 ctrl$maxit, tmp$maxit)), 
      residuals = tmp$resid, 
      model = list(y = tmp$y, x = xmat, w = tmp$w, n = n, p = p, 
	 intercept = intercept, yname = names(mf)[1]), 
      design = design, 
      call = call)

   class(res) <- c("svyreg", "svyregrob")
   res
}

#' @describedIn svyreg.survey.design
svyreg.huber <- function(formula, design, ...){
   UseMethod("svyreg.huber", design)
}

#' Survey regression estimator 
#' 
#' Weighted regression estimator
#'  
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted.  
#' @param design survey design (object of class \code{svydesign})
#' @param object (see \code{summary}, \code{coef}, and \code{fitted} method)
#' @param x (see \code{print} and \code{plot} method)
#' @param which (see  \code{plot} method)
#' @return object of class \code{svyreg} 
svyreg.survey.design <- function(formula, design){
   call <- match.call()
   call[[1]] <- substitute(svyreg)
   res <- svyreg.huber.survey.design(formula, design, k = 10000) 
   res$estimator <- "Survey regression estimator"
   res$call <- call
   res$robust <- NULL 
   res$optim <- NULL
   class(res) <- "svyreg"
   res
}

#' @describedIn svyreg.survey.design
svyreg <- function(formula, design, ...){
   UseMethod("svyreg", design)
}

#' @describedIn svyreg.survey.design
coef.svyreg <- function(object, ...){
   object$estimate
}

#' @describedIn svyreg.survey.design
fitted.svyreg <- function(object, ...){
   object$model$y - object$residuals 
}

#' @describedIn svyreg.survey.design
print.svyreg <- function(x, digits = max(3L, getOption("digits") - 3L),
      ...){
   cat(x$estimator, "\n\nCall:\n")
   print(x$call)
   cat("\nCoefficients:\n")
   print.default(format(coef(x), digits = digits), print.gap = 2L, 
      quote = FALSE)
}

#' @describedIn svyreg.survey.design
summary.svyreg <- function(object, digits = max(3L, getOption("digits") - 3L), 
      ...){
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
   tval <- coef(object) / se
   printCoefmat(cbind(Estimate = coef(object), 'Std. Error' = se, 
      't value' = tval, 'Pr(>|t|)' = 2 * pt(abs(tval), n - p, 
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

#' @describedIn svyreg.huber.survey.design
svyreg.control <- function(tol = 1e-5, maxit = 100, psi = "Huber", ...){
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be 
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi, 
      "Huber" = 0L,
      "asymHuber" = 1L)
   list(tol= unname(tol), maxit = unname(maxit), psi = unname(psi0))
}

#' Regression predictor (GREG) 
#' 
#' regression estimator
#'  
#' @param object fitted regression model (object of class \code{svyreg}) 
#' @param mean_auxiliary vector of population means of the auxiliary variables 
#' @return object of class \code{svystatrob} 
svymean.greg.svyreg <- function(object, mean_auxiliary,...){
   call <- match.call()
   call[[1]] <- substitute(svymean.greg)
   # check dimensions
   if (length(mean_auxiliary) != object$model$p){
      stop("Dimension of argument 'mean_auxiliary' is not correct\n")
   }
   # greg estimate
   w <- object$model$w
   est <- sum(mean_auxiliary * coef(object)) 
   if (object$model$intercept == 0){
      est <- est + sum(w * object$residuals) / sum(w)
   }
   names(est) <- object$model$yname
   # greg variance (using pkg survey)
   design <- object$design
   design$variables <- data.frame(z = object$residuals)
   v <- as.numeric(attr(svymean(~z, design), "var"))
   design$variables <- NULL
   res <- list(characteristic = "mean", estimator = "GREG", estimate = est, 
      variance = v, design = design, call = call)
   class(res) <- "svystat.rob"
   res
}

#' @describedIn svyreg.greg.svyreg
svymean.greg <- function(object, ...){
   UseMethod("svymean.greg", object)
}

#' Robust regression predictor (GREG) 
#' 
#' regression estimator
#'  
#' @param object fitted regression model (object of class \code{svyreg}) 
#' @param mean_auxiliary vector of population means of the auxiliary variables 
#' @param k robustness prediction tuning constant (Huber psi-function) 
#' @return object of class \code{svystatrob} 
svymean.greg.huber.svyregrob <- function(object, mean_auxiliary, k){
   stopifnot(k > 0)
   call <- match.call()
   call[[1]] <- substitute(svymean.greg.huber)
   # check dimensions
   if (length(mean_auxiliary) != object$model$p){
      stop("Dimension of argument 'mean_auxiliary' is not correct\n")
   }
   # robust greg estimate
   w <- object$model$w
   resid_winsorized <- pmin.int(k, pmax.int(-k, object$residuals))
   est <- sum(mean_auxiliary * coef(object)) + sum(w * resid_winsorized) / sum(w)
   names(est) <- object$model$yname
   # greg variance (using pkg survey)
   design <- object$design
   design$variables <- data.frame(z = object$residuals)
   v <- as.numeric(attr(svymean(~z, design), "var"))
   res <- list(characteristic = "mean", estimator = "GREG M-estiamtor", 
      estimate = est, variance = v, design = design, call = call)
   class(res) <- "svystat.rob"
   res
}

svymean.greg.huber <- function(object, ...){
   UseMethod("svymean.greg.huber", object)
}

#' @describedIn svyreg.survey.design
plot.svyreg <- function(x, which = 1:5){
   show <- rep(FALSE, 5)
   show[which] <- TRUE
   par(ask = TRUE)
   n <- length(x$model$y)
   v <- sqrt(x$variance)
   stdres <- x$residuals / v
   if (show[1L]){
      r <- range(c(x$model$y, fitted(x)))
      plot(x$model$y, fitted(x), xlab = x$model$yname, ylab = "fitted", 
	 main = "Response variable vs. fitted", type = "n", ylim = r,
	 xlim = r)
      mtext(text = "(circe size is proportional to sampling weights)", 
	 side = 3, line = 0.5)
      scale <- (diff(range(x$model$y)) / n) / mean(x$model$w) 
      bg <- rep(NA, n)
      if (!is.null(x$robust)){
	 bg[x$robust$robweights < 1] <- "black"
      }
      symbols(x$model$y, fitted(x), circles = scale * x$model$w, 
	 inches = FALSE, add = TRUE, bg = bg)
      abline(0, 1)
      grid(col = "grey65")
   }
   if (show[2L]){
      f <- fitted(x)
      plot(f, stdres, main = "Standarized Residuals vs. Fitted", 
	 xlab = "Fitted value", ylab = "Standarized residual", type = "n")
      panel.smooth(f, stdres)
      grid(col = "grey65")
      abline(h = 0, lty = 2, col = "gray65")
      at <- which(abs(stdres) > 1.96)
      if (length(at) > 0){
	 atboundary <- min(f) == f[at]
	 if (sum(atboundary) > 0 & sum(atboundary) < length(atboundary)){
	    text(f[at[atboundary]], stdres[at[atboundary]], 
	       labels = at[atboundary], cex = 0.8, pos = 4)
	    text(f[at[!atboundary]], stdres[at[!atboundary]], 
	       labels = at[!atboundary], cex = 0.8, pos = 2)
	 }else{
	    text(f[at], stdres[at], labels = at, cex = 0.8, pos = 2)
	 }
      }
   }
   if (show[3L]){
      u <- qqnorm(stdres, main = "Normal QQ-Plot", ylab = "Standardized 
	 residual", xlab = "Theoretical quantile")
      qqline(stdres, lty = 2, col = "gray65") 
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
      xdat <- x$model$x[, (1 + x$model$intercept):model$p]
      tmp <- wBACON(xdat, x$model$Di, verbose = FALSE)
      plot(tmp$dist, stdres, xlab = "Robust distance", ylab = "Standardized 
	 residual", xlim = c(0, 1.05 * max(tmp$dist)), type = "n", xaxs = "i")
      outer <- par("usr")
      rect(outer[1], outer[3], outer[2], outer[4], col = "grey90", border = NA)
      inner <- c(0, -1.96, tmp$cutoff, 1.96)
      rect(inner[1], inner[2], inner[3], inner[4], border = NA, col = "white")
      box()
      grid(col = "grey65")
      points(tmp$dist, stdres)
      at <- c(which(abs(stdres) > 1.96), which(tmp$dist > tmp$cutoff))
      flag <- at[duplicated(at)]
      if (!is.null(flag) > 0){
	 text(tmp$dist[flag], stdres[flag], labels = flag, cex = 0.8, pos = 2)
      }
      at <- unique(at)
      if (!is.null(at)){
	 points(tmp$dist[at], stdres[at])
      }
      abline(h = 0, lty = 2)
      abline(v = tmp$cutoff, lty = 3)
      abline(h = c(-1.96, 1.96), lty = 3)
   }
}




library(survey)
# data(api)
# design <- svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
#
# g <- svyreg(api00 ~ api99, design)
# svymean.greg(g, 1:2)
#
#
# r <- svyreg.huber(api00 ~ api99, design, k = 1.3)
# svymean.greg(r, 1:2)
#
#



library(data.table)
# MU284 population
library(sampling)
data(MU284)
MU284 <- as.data.table(MU284)
# add region size
setkey(MU284, REG)
MU284 <- MU284[MU284[, .N, keyby = REG]]
setkey(MU284, LABEL)
# stratified srswor by region (sample size per stratum = 9)
dat <- as.data.table(sampling:::strata(MU284, "REG", rep(9, 8), method = "srswor"))
dat[, REG := NULL] 
setkey(dat, ID_unit)
dat <- MU284[dat]
dat[, weight := 1 / Prob]
d_strat_MU284 <- svydesign(id = ~1, strata = ~ Stratum, weights = ~weight, 
   fpc = ~N, data = dat)

m <- svyreg(REV84 ~ P85 + P75, d_strat_MU284)
x_pop_mean <- c(1, MU284[, c(mean(P85), mean(P75))])
svymean.greg(m, x_pop_mean)

svymean(~REV84, d_strat_MU284)

r <- svyreg.huber(REV84 ~ P85 + P75, d_strat_MU284, k = 1.34)
svymean.greg.huber(r, x_pop_mean, k = 4)

l <- svyreg(log(REV84) ~ P85 + P75, d_strat_MU284)
plot(l)




d <- d_strat_MU284
d
m <- svyreg(REV84 ~ P85 + P75, d)


summary(m)

l0 <- lm(REV84 ~ P85 + P75, weights = weights(d_strat_MU284), 
   data = d_strat_MU284$variables)

l0 <- lm(REV84 ~ P85 + P75, data = d_strat_MU284$variables)





r <- residuals(l0)
w <- weights(d_strat_MU284)
p <- 3
n <- 72
rdf <- n - p
rss <- sum(w * r^2)
resvar <- rss/rdf
p1 <- 1L:p
Qr <- l0$qr
R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
se <- sqrt(diag(R) * resvar)
est <- z$coefficients[Qr$pivot[p1]]
tval <- est/se
ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
ans$residuals <- r
ans$coefficients <- cbind(Estimate = est, `Std. Error` = se, 
        `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
            lower.tail = FALSE))
