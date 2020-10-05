library(robsurvey)
library(survey)
library(MASS)
library(robustbase)
library(robeth)
#===============================================================================
# wrapper function for the coefficients of the Huber regression M-estimator 
robeth_reg <- function(formula, data, W = TRUE, coef = TRUE) 
{
   # NOTE: robustness tuning constant is fixed at k = 1.345 (default)
   mf <- model.frame(formula, data)
   y <- as.numeric(model.response(mf))
   x <- model.matrix(terms(formula), mf)
   n <- length(y); p <- ncol(x) 
   
   dfvals()			       # set default paramters  
   dfrpar(cbind(x, y), "huber")	       # set parameters for 'huber' M-est.

   wgt <- vector("numeric", n)
   ribet0(wgt)			       # consistency correction: MAD 

   # initial estimate of regression: LS
   tmp <- riclls(x, y)	

   s <- liepsh()
   epsi2 <- s$epsi2; epsip <- s$epsip

   if (W) {
      # initial estimate of covariance 
      cv <- ktaskv(x, f = epsi2 / epsip^2)
      # Huber regression M-estimator (using W-algorithm)
      rob <- rywalg(x, y, tmp$theta, wgt, 
	 cov = cv$cov,		 # initival cov [computed above]
	 psp0 = 1,			 # value of \psi'(0) 
	 expsi = psi,		 # [default] 
	 exchi = chi,		 # [default]
	 exrho = rho,		 # [default]
	 sigmai = tmp$sigma,	 # initial sigma [computed above] 
	 tol = 0.0001,		 # relative precision of convergence crit.
	 gam = 1,			 # relaxation factor 
	 tau = 1e-5,		 # tolerance for determinaton of pseudo rank 
	 itype = 1,		 # 1: Huber, 2: Mallows, 3: Schweppe
	 isigma = 2,		 # estimator of sigma (see RYSIGM)
	 icnv = 1,			 # convergence criterion 
	 maxit = 50,		 # max. iterations of IRWLS steps 
	 maxis = 1,		 # max. iterations for scale step 
	 nitmon = 0)		 # toggle convergence monitoring (nitmon > 0) 
   } else {
      # initial estimate of covariance 
      cv <- kiascv(tmp$xt, fu = epsi2 / epsip^2, fb = 0)
      # Huber regression M-estimator (using H-algorithm)
      rob <- ryhalg(x, y, tmp$theta, wgt, 
	 cov = cv$cov,
	 isigma = 2,		 # MAD
	 sigmai = tmp$sigma,
	 ic = 0)
   }
   # extract coefficients
   if (coef)
      rob$theta[1:p]		
   else
      rob$sigmaf
}

#===============================================================================
data(education, package="robustbase")
design <- svydesign(id = ~1, weights = rep(1, nrow(education)), 
   data = education)
formula <- Y ~ X1 + X2 + X3

robsurvey <- svyreg_huber(formula, design, k = 1.345)
mass <- rlm(formula, education, k = 1.345, method = "M", scale.est = "MAD", 
   acc = 1e-5, test.vec = "coef")
robeth_w <- robeth_reg(formula, education) 
robeth_h <- robeth_reg(formula, education, FALSE) 

coeff <- rbind(coef(robsurvey), robeth_w, robeth_h, coef(mass))
rownames(coeff) <- c("robsurvey", "ROBETH (W-alg.)", "ROBETH (H-alg.)", 
   "MASS (rlm)")
round(coeff, 4)

#===============================================================================
data(stackloss, package = "datasets")
design <- svydesign(id = ~1, weights = rep(1, nrow(stackloss)), 
   data = stackloss)
formula <- stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.

robsurvey <- svyreg_huber(formula, design, k = 1.345)
mass <- rlm(formula, stackloss, k = 1.345, method = "M", scale.est = "MAD", 
   acc = 1e-5, test.vec = "coef")
robeth_w <- robeth_reg(formula, stackloss) 
robeth_h <- robeth_reg(formula, stackloss, FALSE) 

coeff <- rbind(coef(robsurvey), robeth_w, robeth_h, coef(mass))
rownames(coeff) <- c("robsurvey", "ROBETH (W-alg.)", "ROBETH (H-alg.)", 
   "MASS (rlm)")
round(coeff, 4)


#===============================================================================
data(delivery, package = "robustbase")



