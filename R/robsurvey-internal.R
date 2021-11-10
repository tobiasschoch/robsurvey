.check <- function(x, w, na.rm = FALSE, check_NA = TRUE)
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
    if (check_NA) {
        cc <- stats::complete.cases(x, w)
        if (sum(cc) != n) {
            if (na.rm) {
                x <- x[cc]
                w <- w[cc]
            } else {
                return(NULL)
            }
        }
    }

    # number of observations (having removed possible NA's)
    n <- length(x)

    # check if data vector and weights are finite
    if (sum(is.finite(c(x, w))) != 2 * n) {
        warning("Some observations are not finite\n", call. = FALSE,
	        immediate. = TRUE)
        return(NULL)
    }
    list(x = x, w = w, n = n)
}
# check and extract data from survey.design object
.checkformula <- function(f, design, na.rm = FALSE)
{
    if (inherits(f, "formula")) {
        if (length(all.vars(f)) > 1)
            stop("Formula must refer to one variable only\n", call. = FALSE)
        tf <- stats::terms.formula(f)
        if(attr(tf, "response") != 0)
            stop("The LHS of the formula must not be defined\n", call. = FALSE)
        yname <- attr(tf, "term.labels")
        y <- stats::model.frame(f, design$variables,
            na.action = stats::na.pass)[, 1]
    } else {
        if (is.character(f[1]) && length(f) == 1) {
	        if (f %in% names(design$variables)) {
	            y <- design$variables[, f]
	            yname <- f
	        } else {
	            stop("Variable '", f,"' does not exist\n", call. = FALSE)
            }
        } else {
	        stop("Type of argument '", f, "' is not supported\n",
                call. = FALSE)
        }
    }
    w <- as.numeric(1 / design$prob)
    # check for missing values
    failure <- FALSE
    cc <- stats::complete.cases(y, w)
    if (sum(cc) != length(y)) {
        if (na.rm) {
            y <- y[cc]
            w <- w[cc]
            # drop missing values from survey.design object
            design <- subset(design, cc)
        } else {
            failure <- TRUE
        }
    }
    list(failure = failure, yname = yname, y = y, w = w, design = design)
}
# check and extract data from survey.design object (for regression)
.checkreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    if (!inherits(formula, "formula"))
        stop("Argument '", formula, "' must be a formula\n", call. = FALSE)

    # heteroscedasticity (only one variable)
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
# check auxiliary totals and means
.checkauxiliary <- function(object, data, est = "mean", check.names = TRUE,
    na.action = stats::na.omit) #FIXME: default: na.fail
{
    names_beta <- names(object$estimate)
    names_data <- names(data)
    if (is.vector(data))
        if (NCOL(data) < NROW(data))
	        data <- t(data)

    # vector of population x-means or -totals
    if (NROW(data) == 1) {
        # drop intercept (if there is one)
        if (object$model$intercept) {
            names_data <- names_data[-1]
            names_beta <- names_beta[-1]
        }
        if (length(data) != object$model$p)
            stop("Length of auxiliary data does not match\n", call. = FALSE)

        if (check.names && !is.null(names_data)) {
            if (!all(names_data == names_beta))
                stop("Variable names do not match (check.names = TRUE)",
                    call. = FALSE)
        }
        x <- data
    # compute mean/total based on design matrix
    } else {
        mf <- stats::model.frame(object$call$formula, data,
            na.action = na.action)
        xmat <- stats::model.matrix(stats::terms(mf), mf)
        x <- switch(est, "mean" = colMeans(xmat), "total" = colSums(xmat))
    }
    unname(x)
}
# psi functions
.psi_function <- function(x, k, psi = c("Huber", "Huberasym", "Tukey"))
{
    stopifnot(is.numeric(x), k > 0)
    if (is.infinite(k))
        return(x)

    n <- length(x)
    type <- switch(match.arg(psi),
        "Huber" = 0,
        "Huberasym" = 1,
        "Tukey" = 2)
    if (any(is.na(x))) {
        rep(NA, n)
    } else {
        tmp <- .C("psi_function", x = as.double(x), k = as.double(k),
            n = as.integer(n), psi = as.integer(type),
            res = as.double(numeric(n)), PACKAGE = "robsurvey")
        tmp$res
    }
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
					version 0.2 |___/\n\ntype: package?robsurvey to learn more\n")
   }
}
