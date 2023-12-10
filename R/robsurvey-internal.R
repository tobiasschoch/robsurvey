# check data and weight
.check_data_weights <- function(x, w, na.rm = FALSE, check_NA = TRUE)
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
    # empty data
    if (n == 0)
        return(NULL)

    # check for missing values
    if (check_NA) {
        cc <- complete.cases(x, w)
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
.check_formula <- function(f, design, na.rm = FALSE, check_NA = TRUE)
{
    if (inherits(f, "formula")) {
        if (length(all.vars(f)) > 1)
            stop("Formula must refer to one r.h.s. variable only\n",
                call. = FALSE)
        tf <- terms.formula(f)
        if(attr(tf, "response") != 0)
            stop("The LHS of the formula must not be defined\n", call. = FALSE)
        yname <- attr(tf, "term.labels")
        y <- model.frame(f, design$variables, na.action = na.pass)[, 1]
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
    # special return for empty data
    if (length(y) == 0)
        return(list(failure = TRUE))

    # check for missing values
    failure <- FALSE
    if (check_NA) {
        is_not_NA <- complete.cases(y, w)
        if (sum(is_not_NA) != length(y)) {
            if (na.rm) {
                y <- y[is_not_NA]
                w <- w[is_not_NA]
                # drop missing values from survey.design object
                design <- subset(design, is_not_NA)
            } else {
                failure <- TRUE
            }
        }
    }
    list(failure = failure, yname = yname, y = y, w = w, design = design)
}
# check and extract data from survey.design object (for regression)
.check_regression <- function(formula, design, var = NULL, xwgt = NULL,
    na.rm = FALSE)
{
    if (!inherits(formula, "formula"))
        stop("Argument '", formula, "' must be a formula\n", call. = FALSE)

    # heteroscedasticity (only one variable); without NA handling; we will
    # deal with this together with x, y, etc.
    if (!is.null(var)) {
        var <- .check_formula(var, design, FALSE, FALSE)$y
        if (any(var <  .Machine$double.eps))
            stop("Some variances are zero or close to zero\n", call. = FALSE)
    }

    # extract the variables
    mf <- model.frame(formula, design$variables, na.action = na.pass)
    mt <- terms(mf)
    response <- attr(mt, "response")
    if (response == 0)
        stop("The LHS of formula is not defined\n", call. = FALSE)
    yname <- names(mf)[response]
    y <- as.numeric(model.response(mf))
    x <- model.matrix(mt, mf)
    w <- as.numeric(1 / design$prob)
    n <- length(y)

    if (is.null(xwgt))
        xwgt <- rep(1, n)

    # NA treatment
    failure <- FALSE
    is_not_NA <- complete.cases(y, x, w, var, xwgt)
    if (sum(is_not_NA) != n) {
        if (na.rm) {
	        x <- x[is_not_NA, ]
	        y <- y[is_not_NA]
	        w <- w[is_not_NA]
            xwgt <- xwgt[is_not_NA]
            n <- NROW(x)
            p <- NCOL(x)
	        if (!is.null(var))
	            var <- var[is_not_NA]
            # drop missing values from survey.design object
            design <- subset(design, is_not_NA)
            # check if any element is not finite
            if (!is.null(var)) {
                chk <- sum(is.finite(c(x, y, w, var, xwgt))) != (4 + p) * n
                if (any(var <= 0))
                    stop("Some of the variances are <= 0\n", call. = FALSE)
            } else {
                chk <- sum(is.finite(c(x, y, w, xwgt))) != (3 + p) * n
            }
            if (chk)
                stop("Some observations are not finite\n", call. = FALSE)
        } else {
            failure <- TRUE
        }
    }
    list(failure = failure, x = x, y = y, var = var, w = w, terms = mt,
        design = design, xwgt = xwgt, yname = yname)
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
        tmp <- .C(C_psi_function, x = as.double(x), k = as.double(k),
            n = as.integer(n), psi = as.integer(type),
            res = as.double(numeric(n)))
        tmp$res
    }
}
# weight function
.psi_wgt_function <- function(x, k, psi = c("Huber", "Huberasym", "Tukey"))
{
    .psi_function(x, k, psi) / x
}
# onAttach function
.onAttach <- function(libname, pkgname)
{
    quietly <- TRUE

    # this code is taken from Henrik Bengtsson
    # https://github.com/HenrikBengtsson/R.methodsS3/ => pkgStartupMessage.R
    tryCatch({
        # The default, if not found
        quietly <- formals(library)$quietly

        # Identify the environment/frame of interest by making sure
        # it at least contains all the arguments of source().
        argsToFind <- names(formals(library))

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
    packageStartupMessage("\n          88
          88
          88
     e8d88 .d8b.  8888b.   ___ _   _ _ ____   _____ _   _
     8P'  d8' '8b 88 '8b  / __| | | | '__\\ \\ / / _ \\ | | |
     88   Y8. .8P 88  dP  \\__ \\ |_| | |   \\ V /  __/ |_| |
     88    'Y8P'  88e8P'  |___/\\__,_|_|    \\_/ \\___|\\__, |
                                                     __/ |
                                         version 0.6 |___/\n
type: package?robsurvey to learn more
use:  library(robsurvey, quietly = TRUE) to suppress the
      start-up message\n")
   }
}
