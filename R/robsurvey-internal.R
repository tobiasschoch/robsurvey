# some sanity checks (univariate)
.check <- function(x, w, na.rm)
{
    if (missing(w))
        stop("Argument 'w' is missing\n", call. = FALSE)

    if (is.factor(x) || is.factor(w) || is.data.frame(x))
        stop("Arguments data and weights must be numeric vectors\n",
            call. = FALSE)

    n <- length(x)
    if (n != length(w))
        stop("Data vector and weights are not of the same dimension\n",
	        call. = FALSE)
    if (n == 0)
        return(NULL)

    # check for missing values
    cc <- stats::complete.cases(x, w)
    if (sum(cc) != n) {
        if (na.rm) {
	        x <- x[cc]
            w <- w[cc]
        } else {
	        return(NULL)
        }
    }
    n <- length(x)

    # check if data vector and weights are finite
    if (sum(is.finite(c(x, w))) != 2 * n) {
        warning("Some observations are not finite\n", call. = FALSE,
	        immediate. = TRUE)
        return(NULL)
    }
    list(x = x, w = w, n = n)
}

# check and extract data from survey.design object (for regression)
.checkreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    if (!inherits(formula, "formula"))
        stop("Argument '", formula, "' must be a formula\n", call. = FALSE)

    # heteroscedasticity
    if (!is.null(var))
        var <- .checkformula(var, design)$y

    # extract the variables
    mf <- stats::model.frame(formula, design$variables,
        na.action = stats::na.pass)
    mt <- stats::terms(mf)
    response <- attr(mt, "response")
    if (response == 0)
        stop("The LHS of formula is not defined\n", call. = FALSE)
    yname <- names(mf)[response]
    y <- stats::model.response(mf)
    x <- stats::model.matrix(mt, mf)
    w <- as.numeric(1 / design$prob)

    # NA treatment
    cc <- stats::complete.cases(y, x, w, var)
    if (sum(cc) != length(y)) {
        if (na.rm) {
	        x <- x[cc, ]
	        y <- y[cc]
	        w <- w[cc]
	        if (!is.null(var))
	            var <- var[cc]
        } else {
            stop("Data must not contain missing values; see argument 'na.rm'\n",
	            call. = FALSE)
        }
    }
    n <- nrow(x); p <- ncol(x)

    # check if any element is not finite
    if (!is.null(var)) {
        chk <- sum(is.finite(c(x, y, w, var))) != (3 + p) * n
        if (any(var <= 0))
	        stop("Some of the variances are <= 0\n", call. = FALSE)
    } else {
        chk <- sum(is.finite(c(x, y, w))) != (2 + p) * n
    }

    if (chk)
        stop("Some observations are not finite\n", call. = FALSE)

    list(x = x, y = as.numeric(y), yname = yname, var = var, w = w,
        intercept = attr(mt, "intercept"))
}

# check and extract data from survey.design object
.checkformula <- function(f, design)
{
    if (inherits(f, "formula")) {
        mf <- stats::model.frame(f, design$variables,
            na.action = stats::na.pass)
        mt <- stats::terms(mf)
        if (attr(mt, "response") != 0)
	        stop("The LHS of the formula must not be defined\n", call. = FALSE)
        yname <- names(mf)[1]
        if (ncol(mf) > 1)
	        warning("No. of variables > 1; hence, only variable '", yname,
	            "' is considered\n", call. = FALSE)
        y <- mf[[1]]
        x <- NULL
    } else {
        if (is.character(f[1]) && length(f) == 1) {
	        if (f %in% names(design$variables)) {
	            y <- design$variables[, f]
	            x <- NULL
	            yname <- f
	        } else {
	            stop("Variable '", f,"' does not exist\n", call. = FALSE)
            }
        } else {
	        stop("Type of argument '", f, "' is not supported\n", call. = FALSE)
        }
    }
    list(yname = yname, y = as.numeric(y), w = as.numeric(1 / design$prob))
}

# influence function winsorized mean
.infl_winsorized <- function(x, w, LB, UB, wm, ngrid = 401)
{
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    # density estimates at 'LB' and '1 - UB'
    bwd <- KernSmooth::dpik(x, scalest = "minim", level = 2, kernel = "normal",
        canonical = FALSE, gridsize = ngrid, range.x = range(x),
        truncate = TRUE)
    at <- seq(min(x), max(x), length = ngrid)
    nx <- rowsum(c(rep(0, ngrid), w), c(1:ngrid, findInterval(x, at)))
    dens <- KernSmooth::locpoly(rep(1, ngrid), nx * ngrid / (diff(range(x)) *
        sum(w)), binned = TRUE, bandwidth = bwd, range.x = range(x))
    f_LB <- dens$y[floor(ngrid * LB) + 1]
    f_UB <- dens$y[floor(ngrid * UB)]
    # influence function
    infl <- pmin.int(qs[2] + (1 - UB) / f_UB, pmax.int(qs[1] - LB / f_LB, x))
    w <- (UB - LB) * wm + LB * qs[1] + (1 - UB) * qs[2]
    infl - w + LB^2 / f_LB + (1 - UB)^2 / f_UB
}

# influence function trimmed mean
.infl_trimmed <- function(x, w, LB, UB, tm)
{
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    # influence function
    infl <- pmin.int(qs[2], pmax.int(qs[1], x))
    w <- (UB - LB) * tm + LB * qs[1] + (1 - UB) * qs[2]
    (infl - w) / (UB - LB)
}

# onAttach function
.onAttach <- function(libname, pkgname)
{
    quietly <- TRUE

    # this code is taken from Henrik Bengtsson
    # https://github.com/HenrikBengtsson/R.methodsS3/ => pkgStartupMessage.R
    tryCatch({
        # The default, if not found
        quietly <- formals(base::library)$quietly

        # Identify the environment/frame of interest by making sure
        # it at least contains all the arguments of source().
        argsToFind <- names(formals(base::library))

        # Scan the call frames/environments backwards...
        srcfileList <- list()
        for (ff in sys.nframe():0) {
            env <- sys.frame(ff)

            # Does the environment look like a library() environment?
            exist <- sapply(argsToFind, FUN = exists, envir = env,
                inherits = FALSE)

            if (!all(exist))
                next

            # Was argument 'quietly' specified?
            missing <- eval(expression(missing(quietly)), envir = env)
            if (!missing) {
                quietly <- get("quietly", envir = env, inherits = FALSE)
                break
            }
            # ...otherwise keep searching due to nested library() calls.
        }
    }, error = function() {})

    if (!quietly){
    packageStartupMessage("\n                  88
		  88
		  88
     e8d88 .d8b.  8888b.   ___ _   _ _ ____   _____ _   _
     8P'  d8' '8b 88 '8b  / __| | | | '__\\ \\ / / _ \\ | | |
     88   Y8. .8P 88  dP  \\__ \\ |_| | |   \\ V /  __/ |_| |
     88    'Y8P'  88e8P'  |___/\\__,_|_|    \\_/ \\___|\\__, |
						     __/ |
					version 0.2 |___/\n")
   }
}
