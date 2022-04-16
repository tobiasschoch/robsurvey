# NOTE: robustness tuning constant is fixed at k = 1.345 (default) for all
#       functions

# wrapper function for the coefficients of the Huber regression M-estimator
robeth_reg <- function(formula, data, tol = 1e-5, maxit = 50, W = TRUE)
{
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y)
    p <- ncol(x)

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
	        tol = tol,		    # relative precision of convergence crit.
	        gam = 1,			# relaxation factor
	        tau = 1e-5,		    # tolerance for determinaton of pseudo rank
	        itype = 1,		    # 1: Huber, 2: Mallows, 3: Schweppe
	        isigma = 2,		    # estimator of sigma (see RYSIGM)
	        icnv = 1,			# convergence criterion
	        maxit = maxit,      # max. iterations of IRWLS steps
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
    n      <- nrow(x)
    np <- ncol(x)
    dfcomn(xk=np)
    .dFvPut(1,"itw")
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y, exu = usext)
    nitw <- z$nit
    wgt <- 1 / z$dist
    wgt[wgt>1.e6] <- 1.e6
    z <- comval()
    bto <- z$bt0
    ipso <- z$ipsi
    co <- z$c
    z <- ribet0(wgt, itype = 2, isqw = 0)
    xt <- x * wgt
    yt <- y * wgt
    z <- rilars(xt, yt)
    theta0 <- z$theta
    sigma0 <- z$sigma
    rs <- z$rs / wgt
    r1 <- rs/sigma0
    dfcomn(ipsi = 1,c = cc)
    z <- liepsh(cc)
    den <- z$epsip
    g <- Psp(r1) / den # (see Psi in Chpt. 14)
    dfcomn(ipsi = ipso, c = co, bet0 = bto)
    list(theta = theta0, sigma = sigma0, rs = rs, g = g, nitw = nitw)
}

# ROBETH: Mallows-type standard estimator (with wyfalg and rywalg)
robeth_mallows <- function(formula, data, tol = 1e-5, maxit = 50){
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y)
    np <- ncol(x)
    b2 <- -1
    cc <- -1

    isigma <- 2                 # scale estimated by MAD
    dfrpar(x, "Mal-Std", b2, cc)

    # Weights
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y)
    nitw <- z$nit
    wgt <- Www(z$dist)          # See Www in Chpt. 14
    # Initial cov. matrix of coefficient estimates
    z <- kiedch(wgt)
    cov <- ktaskw(x, z$d, z$e, f = 1 / n)
    # Initial theta and sigma
    z <- rbmost(x, y, 1.5, userfd)
    theta0 <- z$theta
    sigma0 <- z$sigma
    nitw0 <- z$nitw
    # Final theta and sigma
    ribet0(wgt)
    z <- rywalg(x, y, theta0, wgt, cov$cov, sigmai = sigma0, tol = tol,
        maxit = maxit, itype = 2, isigma = 2, icnv = 1, maxis = 1)
    theta1 <- z$theta[1:np]
    sigma1 <- z$sigmaf
    nit1 <- z$nit
    # covariance estimate
    z <- kfedcc(wgt, z$rs, sigma = sigma1)
    cov <- ktaskw(x, z$d, z$e, f = sigma1^2 / n)

    list(coef = theta1, scale = sigma1, wgt = wgt, cov = tri_complete(cov$cov))
}

# Schweppe type estimator
robeth_schweppe <- function(formula, data, maxit = 50, tol = 1e-5){
    mf <- model.frame(formula, data)
    y <- as.numeric(model.response(mf))
    x <- model.matrix(terms(formula), mf)
    n <- length(y)
    np <- ncol(x)
    b2 <- -1
    cc <- -1
    isigma <- 2 # scale estimated by MAD
    dfrpar(x, "Kra-Wel", b2, cc)
    # Weights
    z <- wimedv(x)
    z <- wyfalg(x, z$a, y)
    nitw <- z$nit
    wgt <- Www(z$dist)   # See Www in Chpt. 14
    # Initial cov. matrix of coefficient estimates
    z <- kiedch(wgt)
    cov <- ktaskw(x, z$d, z$e, f = 1 / n)
    # Initial theta and sigma
    z <- rbmost(x, y, 1.5, userfd)
    theta0 <- z$theta
    sigma0 <- z$sigma
    nitw0 <- z$nitw
    # Final theta and sigma
    ribet0(wgt)
    z <- rywalg(x, y, theta0, wgt, cov$cov, sigmai = sigma0, tol = tol,
        maxit = maxit, itype = 3, isigma = 2, icnv = 1, maxis = 1)
    theta1 <- z$theta[1:np]
    sigma1 <- z$sigmaf
    nit1 <- z$nit

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
# M-estimator: compare our implementation with ROBETH and MASS
M_compare <- function(formula, data, digits = 3, tol = 1e-5,
    maxit = 50)
{
    design <- svydesign(id = ~1, weights = rep(1, nrow(data)),
        data = data)

    # compute the different estimators
    robsurvey1 <- svyreg_huberM(formula, design, k = 1.345,
        mad_center = TRUE, tol = tol, maxit = maxit)
    robsurvey2 <- svyreg_huberM(formula, design, k = 1.345,
        mad_center = FALSE, tol = tol, maxit = maxit)

    mass <- rlm(formula, data, k = 1.345, method = "M", scale.est = "MAD",
        acc = tol, maxit = maxit, test.vec = "coef")

    robeth_w <- robeth_reg(formula, data, tol = tol, maxit = maxit)

    coeff <- rbind(
        c(coef(robsurvey1), robsurvey1$scale),
        c(coef(robsurvey2), robsurvey2$scale),
        robeth_w,
        c(coef(mass), mass$s))

    colnames(coeff)[NCOL(coeff)] <- "scale"
    rownames(coeff) <- c("svyreg_huberM", "svyreg_huberM (mad0)",
        "rywalg (ROBETH)", "rlm (MASS)")

    print(round(coeff, digits))
}

# Mallows GM-estimator: compare our implementation with ROBETH
GM_mallows_compare <- function(formula, data, digits = 3, tol = 1e-5,
    maxit = 50)
{
    # ROBETH Mallows
    robeth <- robeth_mallows(formula, data)
    wgt_delivery_mallow <- robeth$wgt

    # our implementation
    design <- svydesign(id = ~1, weights = rep(1, nrow(data)),
        data = data)
    robsurvey1 <- svyreg_huberGM(formula, design, k = 1.345,
        xwgt = wgt_delivery_mallow, type = "Mallows", mad_center = TRUE)
    robsurvey2 <- svyreg_huberGM(formula, design, k = 1.345,
        xwgt = wgt_delivery_mallow, type = "Mallows", mad_center = FALSE)

    coeff <- rbind(
        c(coef(robsurvey1), robsurvey1$scale),
        c(coef(robsurvey2), robsurvey2$scale),
        c(robeth$coef, robeth$scale))

    colnames(coeff)[NCOL(coeff)] <- "scale"
    rownames(coeff) <- c(
        "svyreg_huberGM (Mallows)",
        "svyreg_huberGM (Mallows, mad0)",
        "rywalg (ROBETH, Mallows)")

    print(round(coeff, digits))
}

# Schweppe GM-estimator: compare our implementation with ROBETH
GM_schweppe_compare <- function(formula, data, digits = 3, tol = 1e-5,
    maxit = 50)
{
    # ROBETH Schweppe
    robeth <- robeth_schweppe(formula, data)
    wgt_delivery_schweppe <- robeth$wgt

    # our implementation
    design <- svydesign(id = ~1, weights = rep(1, nrow(data)),
        data = data)
    robsurvey1 <- svyreg_huberGM(formula, design, k = 1.345,
        xwgt = wgt_delivery_schweppe, type = "Schweppe", mad_center = TRUE)
    robsurvey2 <- svyreg_huberGM(formula, design, k = 1.345,
        xwgt = wgt_delivery_schweppe, type = "Schweppe", mad_center = FALSE)

    coeff <- rbind(
        c(coef(robsurvey1), robsurvey1$scale),
        c(coef(robsurvey2), robsurvey2$scale),
        c(robeth$coef, robeth$scale))

    colnames(coeff)[NCOL(coeff)] <- "scale"
    rownames(coeff) <- c(
        "svyreg_huberGM (Schweppe)",
        "svyreg_huberGM (Schweppe, mad0)",
        "rywalg (ROBETH, Schweppe)")

    print(round(coeff, digits))
}

M_compare_cov <-  function(formula, data, digits = 3, tol = 1e-5,
    maxit = 50)
{
    # robsurvey
    design <- svydesign(id = ~1, weights = rep(1, nrow(data)),
        data = data)
    robsurvey <- svyreg_huberM(formula, design, k = 1.345,
        mad_center = FALSE, tol = tol, maxit = maxit)
    cov_robsurvey <- diag(vcov(robsurvey))

    # MASS
    rlm_mod <- MASS:::rlm.default
    body(rlm_mod)[[22]][[4]][[3]][[3]][[2]][[3]][[3]][[3]] <-
        substitute(median(abs(resid)) * 1.482602)
    mass <- rlm_mod(robsurvey$model$x, robsurvey$model$y, k = 1.345,
        method = "M", scale.est = "MAD", acc = tol, maxit = maxit,
        test.vec = "coef")
    cov_mass <- diag(vcov(mass))

    print(cov_robsurvey)
    print(cov_mass)

    # absolute relative difference
    cat("\nabs_rel_DIFF: ",
        100 * max(abs(cov_mass / cov_robsurvey -  1)), "%\n")
}
