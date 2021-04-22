# wrapper function for the coefficients of the Huber regression M-estimator
# NOTE: robustness tuning constant is fixed at k = 1.345 (default)
robeth_reg <- function(formula, data, W = TRUE)
{
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y); p <- ncol(x)

    dfvals()                        # set default paramters
    dfrpar(cbind(x, y), "huber")    # set parameters for 'huber' M-est.
    wgt <- vector("numeric", n)     # weights
    ribet0(wgt)			            # consistency correction: MAD

    # initial estimate of regression: LS
    tmp <- riclls(x, y)

    s <- liepsh()
    epsi2 <- s$epsi2
    epsip <- s$epsip

    if (W) {
        # initial estimate of covariance
        cv <- ktaskv(x, f = epsi2 / epsip^2)
        # Huber regression M-estimator (using W-algorithm)
        rob <- rywalg(x, y, tmp$theta, wgt,
	        cov = cv$cov,		# initival cov [computed above]
	        psp0 = 1,			# value of \psi'(0)
	        expsi = psi,		# [default]
	        exchi = chi,		# [default]
	        exrho = rho,		# [default]
	        sigmai = tmp$sigma,	# initial sigma [computed above]
	        tol = 0.0001,		# relative precision of convergence crit.
	        gam = 1,			# relaxation factor
	        tau = 1e-5,		    # tolerance for determinaton of pseudo rank
	        itype = 1,		    # 1: Huber, 2: Mallows, 3: Schweppe
	        isigma = 2,		    # estimator of sigma (see RYSIGM)
	        icnv = 1,			# convergence criterion
	        maxit = 50,		    # max. iterations of IRWLS steps
	        maxis = 1,		    # max. iterations for scale step
	        nitmon = 0)		    # toggle convergence monitoring (nitmon > 0)
    } else {
        # initial estimate of covariance
        cv <- kiascv(tmp$xt, fu = epsi2 / epsip^2, fb = 0)
        # Huber regression M-estimator (using H-algorithm)
        rob <- ryhalg(x, y, tmp$theta, wgt,
	        cov = cv$cov,
	        isigma = 2,         # MAD
	        sigmai = tmp$sigma, # sigma
  	        ic = 0)
    }
    # extract coefficients
    c(rob$theta[1:p], rob$sigmaf)
}

# ROBETH: this function is used to compute initial estimates
rbmost <- function(x, y, cc, usext = userfd) {
    n      <- nrow(x); np <- ncol(x); dfcomn(xk=np)
    .dFvPut(1,"itw")
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y, exu = usext); nitw <- z$nit
    wgt <- 1 / z$dist; wgt[wgt>1.e6] <- 1.e6
    z <- comval()
    bto <- z$bt0;    ipso   <- z$ipsi; co <- z$c
    z <- ribet0(wgt, itype = 2, isqw = 0)
    xt <- x * wgt;    yt     <- y * wgt
    z <- rilars(xt, yt)
    theta0 <- z$theta;  sigma0 <- z$sigma
    rs <- z$rs / wgt; r1     <- rs/sigma0
    dfcomn(ipsi = 1,c = cc)
    z <- liepsh(cc)
    den <- z$epsip
    g <- Psp(r1) / den          # (see Psi in Chpt. 14)
    dfcomn(ipsi = ipso, c = co, bet0 = bto)
    list(theta = theta0, sigma = sigma0, rs = rs, g = g, nitw = nitw)
}

# ROBETH: Mallows-type standard estimator (with wyfalg and rywalg)
robeth_mallows <- function(formula, data){
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y); np <- ncol(x)
    b2 <- -1; cc <- -1

    isigma <- 2                 # scale estimated by MAD
    dfrpar(x, "Mal-Std", b2, cc)

    # Weights
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y); nitw <- z$nit
    wgt <- Www(z$dist)          # See Www in Chpt. 14
    # Initial cov. matrix of coefficient estimates
    z <- kiedch(wgt)
    cov <- ktaskw(x, z$d, z$e, f = 1 / n)
    # Initial theta and sigma
    z <- rbmost(x, y, 1.5, userfd)
    theta0 <- z$theta; sigma0 <- z$sigma; nitw0 <- z$nitw
    # Final theta and sigma
    ribet0(wgt)
    z <- rywalg(x, y, theta0, wgt, cov$cov, sigmai = sigma0)
    theta1 <- z$theta[1:np]; sigma1 <- z$sigmaf; nit1 <- z$nit
    # covariance estimate
    z <- kfedcc(wgt, z$rs, sigma = sigma1)
    cov <- ktaskw(x, z$d, z$e, f = sigma1^2 / n)

    list(coef = theta1, scale = sigma1, wgt = wgt, cov = tri_complete(cov$cov))
}

# Schweppe type estimator
robeth_schweppe <- function(formula, data){
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y); np <- ncol(x)
    b2 <- -1; cc <- -1
    isigma <- 2 # scale estimated by MAD
    dfrpar(x, "Kra-Wel", b2, cc)
    # Weights
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y); nitw <- z$nit
    wgt <- Www(z$dist)   # See Www in Chpt. 14
    # Initial cov. matrix of coefficient estimates
    z <- kiedch(wgt)
    cov <- ktaskw(x, z$d, z$e, f = 1 / n)
    # Initial theta and sigma
    z <- rbmost(x, y, 1.5, userfd)
    theta0 <- z$theta; sigma0 <- z$sigma; nitw0 <- z$nitw
    # Final theta and sigma
    ribet0(wgt)
    z <- rywalg(x, y, theta0, wgt, cov$cov, sigmai = sigma0, tol = 0.00001)
    theta1 <- z$theta[1:np]; sigma1 <- z$sigmaf; nit1 <- z$nit

    #FIXME: this cov may be wrong...

    # covariance estimate
    z <- kfedcc(wgt, z$rs, sigma = sigma1)
    cov <- ktaskw(x, z$d, z$e, f = sigma1^2 / n)
    list(coef = theta1, scale = sigma1, wgt = wgt, cov = tri_complete(cov$cov))
}

# triangular matrix
tri <- function(x, lower = TRUE){
    p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    b <- matrix(0, ncol = p, nrow = p)
    if (lower)
        b[lower.tri(b) | row(b) == col(b)] <- x
    else
        b[upper.tri(b) | row(b) == col(b)] <- x
    b
}

# complete a triangular matrix to be a 'full' matrix
tri_complete <- function(x, lower = FALSE){
    b <- tri(x, lower)
    c <- b + t(b)
    diag(c) <- diag(b)
    c
}
