sapply(c("robsurvey", "survey", "MASS"), require, character.only = TRUE)

# reference function
test_svyreg_huber <- function(formula, design, k, tol = 1e-5, mad0 = FALSE) 
{
   dat <- robsurvey:::.checkreg(formula, design, NULL)
   x <- dat$x; y <- dat$y; w <- dat$w
   n <- NROW(y); p <- NCOL(x)

   # initialize the regression estimate
   tmp <- lm.wfit(x, y, w)
   beta0 <- tmp$coefficients; resid <- tmp$resid 

   # initialize the scale
   scale0 <- if (mad0)
	 robsurvey::weighted_median(abs(resid), w) / 0.6745 
      else
	 robsurvey::weighted_mad(resid, w)   

   while (1) {
      # update the regression estimate
      robwgt <- huberWgt(resid / scale0, k)
      tmp <- lm.wfit(x, y, robwgt * w)
      beta <- tmp$coefficients; resid <- tmp$residuals

      # update the scale
      scale <- if (mad0)
	 robsurvey::weighted_median(abs(resid), w) / 0.6745 
      else
	 robsurvey::weighted_mad(resid, w)   

      # check for convergence
      if (norm(as.matrix(beta - beta0)) < tol * scale) 
	 break
      else 
	 beta0 <- beta; scale0 <- scale
   }
  
   # covariance
   R <- qr(sqrt(w) * x)$qr[1L:p, 1L:p, drop = FALSE]
   R[lower.tri(R)] <- 0
   Rinv <- solve(R, diag(p))
   Sigma <- tcrossprod(Rinv, Rinv)

   psiprime <- abs(resid / scale) <= k
   sum_w <- sum(w)
   Epsi <- sum(w * psiprime) / sum_w
   Epsi2 <- sum(w * psiprime^2) / sum_w

   S <- sum(w * (resid * robwgt)^2) / (sum_w - p) 
   kappa <- 1 + p / sum(w) * (Epsi2 / Epsi^2  - 1) * n / (n - 1) 
   stddev <- kappa * sqrt(S) / Epsi

   list(coef = beta, scale = as.numeric(scale), kappa = kappa, stddev = stddev, 
      sigma = Sigma * stddev^2, Epsi = Epsi, Epsi2 = Epsi2, resid = resid, 
      robwgt = robwgt)
}

#===============================================================================
data(workplace)
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight, 
   data = workplace)
f <- payroll ~ employment

test_svyreg_huber(f, dn, k = 1.345, tol = 1e-8)
summary(svyreg_huber(f, dn, k = 1.345, tol = 1e-8))

# test_svyreg_huber(f, dn, k = 1.345, tol = 1e-4, mad0 = TRUE)
# r <- rlm(formula = f, data = workplace, weights = as.numeric(weights(dn)), 
#     k = 1.345, scale.est = "MAD", wt.method = "case")

