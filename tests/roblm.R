setwd("C:/My/code/robsurvey/src")
dyn.load("robsurvey.dll")

huber_wgt <- function(x, k) pmin.int(1, k / abs(x))
huber_psi <- function(x, k) pmin.int(k, pmax.int(-k, x))
huber_psi_prime <- function(x, k) abs(x) <= k

w_median <- function(x, w = NULL){
   n <- length(x)
   if (is.null(w)) w <- rep(1, n)
   .C("wquantile", x = as.double(x), w = as.double(w), n = as.integer(n), 
      p = as.double(0.5), med = as.double(numeric(1)))$med
}

w_mad <- function(x, w = NULL){
   if (is.null(w)) w <- rep(1, length(x))
   w_median(abs(x - w_median(x, w)), w) * 1.4826
}

# NOTE: 
# MASS:rlm(x, y, weights = w, k = 2, init = "ls", scale.est = "MAD", 
#	   method = "M", wt.method = "case", test.vec = "coef")
# MASS: mad about zero (not about median of abs. residuals => can be a problem
#       when the model does not have an intercept, see A. Welsh, Ann. Stat, 
#       1986)
# MASS: different termination rule

roblm <- function(x, y, w, k){
   tol <- 1e-4; iter <- 1 
   n <- length(y)

   # inital estiamtes
   tmp <- lm.wfit(x, y, w)
   beta0 <- tmp$coefficients
   res <- tmp$residuals
   s0 <- w_mad(res, w)
   
   while (iter < 50){
      iter <- iter + 1
      ui <- huber_wgt(res / s0, k)

      tmp <- lm.wfit(x, y, w * ui)
      beta <- tmp$coefficients
      res <- tmp$residuals

      s <- w_mad(res, w)

      if (norm(as.matrix(beta - beta0), type = "F") < s * tol){
	 break
      }

      beta0 <- beta
      s0 <- s
   }

   # u_i's with (normalized) weighted mad as estimate of scale 
   ui <- huber_wgt(res / s, k)

   # compute "Proposal 2" variance estimate 
   kappa <- ifelse(k < 10, 1 - 2 * (k * dnorm(k) + (1 - k * k) * pnorm(k, 
      lower = FALSE)), 1)
   scale <- sqrt(sum(w * (ui * res)^2) / (n * kappa))

   # empirical estimates of Epsi2 and Epsi_prime
   Epsi2 <- sum(huber_psi(res / scale, k)^2) / n
   Epsi_prime <- sum(huber_psi_prime(res / scale, k)) / n

   list(beta = beta, scale = scale, Epsi2 = Epsi2, Epsi_prime = Epsi_prime)
}





roblm2 <- function(x, y, w, k){
   n <- length(y)
   p <- ncol(x)
   maxit = 50
   tol = 1e-4


   list(beta = tmp$beta, scale = tmp$scale)
}



data(iris)
y <- iris[, 1]
n <- length(y)
x <- as.matrix(cbind(rep(1, n), iris[, 2:3]))
colnames(x)[1] <- "(Intercept)" 
w <- rep(1, n)
lm.wfit(x, y, w)
#m1 <- roblm(x, y, w, k = 2)
#m2 <- roblm2(x, y, w, k = 2)

# variance of MLE variance estimator (i.e. without accounting for the loss of
# degrees of freedom 
# sqrt(sum(lm.fit(x, y)$residuals^2) / length(y))

svyreg_control <- function(tol = 1e-5, maxit = 100, psi = "Huber", k_Inf = 1e5,
   ...)
{
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be 
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi, "Huber" = 0L, "asymHuber" = 1L)
   list(tol= unname(tol), maxit = unname(maxit), psi = unname(psi0), 
      k_Inf = unname(k_Inf))
}




x <- cbind(x, x[,3])

   var <- NULL
   k <- Inf#1.345 #k <- 4.68
   ctrl <- svyreg_control()
   if (k <= 0) stop("Argument k must be > 0\n", call. = FALSE)
   if (k == Inf) k <- ctrl$k_Inf 

   n <- length(y); p <- ifelse(is.null(ncol(x)), 1, ncol(x))
   # account for heteroscedasticity
   if (!is.null(var)){
      x <- x / sqrt(var); y <- y / sqrt(var)
   }
   #
   tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), 
      w = as.double(w), resid = as.double(numeric(n)), 
      robwgt = as.double(numeric(n)), n = as.integer(n), p = as.integer(p), 
      k = as.double(k), beta = as.double(numeric(p)), 
      scale = as.double(numeric(1)), maxit = as.integer(ctrl$maxit), 
      tol = as.double(ctrl$tol), psitype = as.integer(0), 
      Epsi2 = as.double(numeric(1)), Epsiprime = as.double(numeric(1)))
   tmp$beta




