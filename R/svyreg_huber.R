# robust Huber M-estimator of regression (depends on pkg survey)
svyreg_huber <- function(formula, design, k, var = NULL, na.rm = FALSE,
    asym = FALSE, verbose = TRUE, ...)
{
    dat <- .checkreg(formula, design, var, na.rm)
    res <- robsvyreg(dat$x, dat$y, dat$w, k, ifelse(asym, 1, 0), 0, NULL,
        dat$var, verbose, ...)
    res$design <- design
    res$call <- match.call()
    res$model$intercept <- dat$intercept
    res$model$yname <- dat$yname
    class(res) <- "svyreg_rob"
    res
}
# robust Huber GM-estimator of regression (depends on pkg survey)
svyreg_huberGM <- function(formula, design, k, type, xwgt, var = NULL,
    na.rm = FALSE, asym = FALSE, verbose = TRUE, ...)
{
    type_int <- switch(toupper(type), "MALLOWS" = 1, "SCHWEPPE" = 2)
    if (is.null(type_int))
        stop("Type '", type,"' is not defined\n", class. = FALSE)

    dat <- .checkreg(formula, design, var, na.rm)

    if (NCOL(xwgt) > 1) {
        xwgt <- as.numeric(xwgt[, 1])
        warning("Only first column of argument 'xwgt' is used\n",
            class. = FALSE)
    }
    if (length(xwgt) != length(dat$y))
        stop("Argument 'xwgt' is not of length n\n", class. = FALSE)

    res <- robsvyreg(dat$x, dat$y, dat$w, k, ifelse(asym, 1, 0), type_int,
        xwgt, dat$var, verbose, ...)
    res$design <- design
    res$call <- match.call()
    res$model$intercept <- dat$intercept
    res$model$yname <- dat$yname
    res$model$xwgt <- xwgt
    class(res) <- "svyreg_rob"
    res
}
