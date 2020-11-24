GM <- function(formula, data, xwgt, zero = FALSE, Mallows = TRUE, tol = 1e-5, 
   low = FALSE, high = FALSE, k = 1.345) 
{
   mf <- model.frame(formula, data)
   y <- as.numeric(model.response(mf))
   x <- model.matrix(terms(formula), mf)  
   n <- NROW(y); p <- NCOL(x)

   w <- rep(1, n); sum_w <- sum(w)
   tmp <- lm.fit(x, y)
   beta0 <- tmp$coefficients 
   resid <- tmp$resid 
   scale0 <- mad(resid)   
 
   if (Mallows) { # determine consistency correction for Mallows type estimator
      foo <- function(x, wgt) sum(pnorm(x, 0, sqrt(wgt))) / length(wgt) - 0.75 
      const <- uniroot(foo, c(0.1, 5), xwgt)$root
   } else
      const <- 0.6744898

   while (1) {
      robwgt <- if (Mallows)
	    huberWgt(resid / scale0, k)
	 else
	    huberWgt(resid / (xwgt * scale0), k)

      tmp <- lm.wfit(x, y, robwgt)
      beta <- tmp$coefficients 
      resid <- tmp$residuals

      r <- if (Mallows) 
	 resid * sqrt(xwgt)
      else 
	 resid

      scale <- if (zero)
	 mad(r, center = 0, constant = 1 / const, high = TRUE)
      else
	 mad(r, constant = 1 / const)

      if (norm(as.matrix(beta - beta0)) < tol * scale) 
	 break
      else {
	 beta0 <- beta; scale0 <- scale
      }
   }

   # covariance
   if (Mallows) {
      psiprime <- huberPsi@Dpsi(resid / scale, k)
      Epsi <- sum(w * psiprime) / sum_w
      Epsi2 <- sum(w * psiprime^2) / sum_w
      S <- sum(w * (resid * robwgt)^2) / (sum_w - p) / Epsi^2 
      qr <- qr(sqrt(w * xwgt) * x)
      R <- qr.R(qr)
      Qw <- sqrt(xwgt) * qr.Q(qr)
      Rinv <- backsolve(R, diag(p))
      Sigma <- crossprod(Qw %*% t(Rinv))
   } else { 
      s1 <- 
      s2 <-   
   }
   list(coef = beta, scale = scale, sigma = Sigma * S)
}



H <- function(formula, data, zero = FALSE, tol = 1e-5, low = FALSE, 
   high = FALSE) 
{
   mf <- model.frame(formula, data)
   y <- as.numeric(model.response(mf))
   x <- model.matrix(terms(formula), mf)  
   n <- NROW(y); p <- NCOL(x)
   w <- rep(1, n); sum_w <- sum(w)

   tmp <- lm.fit(x, y)
   beta0 <- tmp$coefficients 
   resid <- tmp$resid 
   scale0 <- mad(resid)   

   while (1) {
      robwgt <- huberWgt(resid / scale0)
      tmp <- lm.wfit(x, y, robwgt)
      beta <- tmp$coefficients 
      resid <- tmp$residuals

      scale <- if (zero)
	 median(abs(resid))/0.6745 #mad(resid, center = 0, high = TRUE)
      else
	 mad(resid)

      if (norm(as.matrix(beta - beta0)) < tol * scale) 
	 break
      else {
	 beta0 <- beta; scale0 <- scale
      }
   }
  
   # covariance
   R <- qr(x)$qr[1L:p, 1L:p, drop = FALSE]
   R[lower.tri(R)] <- 0
   Rinv <- solve(R, diag(p))
   Sigma <- tcrossprod(Rinv, Rinv)

   psiprime <- huberPsi@Dpsi(resid / scale)
   Epsi <- sum(w * psiprime) / sum_w
   Epsi2 <- sum(w * psiprime^2) / sum_w

   S <- sum(w * (resid * robwgt)^2) / (sum_w - p) 
   kappa <- 1 + p / sum(w) * (Epsi2 / Epsi^2  - 1) * n / (n - 1) 
   stddev <- kappa * sqrt(S) / Epsi

   list(coef = beta, scale = scale, kappa = kappa, stddev = stddev, sigma = Sigma * stddev^2)
}


library(robsurvey); library(survey); library(robustbase)
data(education, package = "robustbase")
design <- svydesign(id = ~1, weights = rep(1, nrow(education)), 
      data = education)
f <- Y ~ X1 + X2 + X3

H(f, education)

svyreg_huber(f, design, k = 1.345)

svyreg_huberGM(f, design, k = 1.345, type = "Mallows", xwgt = rep(1, 50))

GM(f, education, Mallows = TRUE, xwgt = rep(1, 50))


tmp <- svyreg_huberGM(f, design, k = 1.345, type = "Mallows", xwgt = rep(1, 50))
dat <- list(x = tmp$model$x, y = tmp$model$y, w = tmp$model$w, var = NULL)
res <- test(dat$x, dat$y, dat$w, k = 1.345, 0, 1, xwgt = rep(1, 50), dat$var)

res$variance$covmat

R <- qr.R(qr(dat$x)); Q <- qr.Q(qr(dat$x)); Rinv <- solve(R)

test <- function (x, y, w, k, psi, type, xwgt = NULL, var = NULL, ...) 
{
    ctrl <- svyreg_control(...)
    if (k <= 0) 
        stop("Argument 'k' must be > 0\n", call. = FALSE)
    if (k == Inf) 
        k <- ctrl$k_Inf
    n <- length(y)
    p <- NCOL(x)
    if (is.null(xwgt)) 
        xwgt <- rep(1, n)
    if (!is.null(var)) {
        x <- x/sqrt(var)
        y <- y/sqrt(var)
    }
    tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), w = as.double(w), 
        resid = as.double(numeric(n)), robwgt = as.double(numeric(n)), 
        xwgt = as.double(xwgt), n = as.integer(n), p = as.integer(p), 
        k = as.double(k), beta = as.double(numeric(p)), scale = as.double(numeric(1)), 
        scale2 = as.double(numeric(1)), tol = as.double(ctrl$tol), 
        maxit = as.integer(ctrl$maxit), psi = as.integer(psi), 
        type = as.integer(type), PACKAGE = "robsurvey")
    psi_fun <- switch(psi + 1, "Huber", "asymHuber", "Tukey")
    names(tmp$beta) <- colnames(x)
    list(characteristic = "regression", estimator = paste0("Survey regression ", 
        switch(type + 1, "", "Mallows G", "Schweppe G"), "M-estimator (", 
        psi_fun, " psi, k = ", k, ")"), estimate = tmp$beta, 
        variance = list(
#	       covmat = matrix(tmp$x[1:(p * p)], ncol = p), 
	       covmat = matrix(tmp$x, ncol = p), 
            scale = tmp$scale2), robust = list(psifunction = psi_fun, 
            k = k, robweights = tmp$robwgt, outliers = 1 * (tmp$robwgt < 
                1), Epsi2 = tmp$Epsi2, Epsiprime = tmp$Epsiprime, 
            scale = tmp$scale), optim = list(converged = (tmp$maxit != 
            0), niter = ifelse(tmp$maxit == 0, ctrl$maxit, tmp$maxit), 
            tol = ctrl$tol), residuals = tmp$resid, model = list(x = x, 
            y = y, w = w, var = var, n = n, p = p), design = NA, 
        call = NA)
}

