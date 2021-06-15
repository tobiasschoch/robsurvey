# robust Tukey biweight M-estimator of regression (depends on pkg survey)
svyreg_tukey <- function(formula, design, k, var = NULL, na.rm = FALSE,
    verbose = TRUE, ...)
{
    dat <- .checkreg(formula, design, var, na.rm)
    res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, 0, NULL, dat$var, verbose, ...)
    # add a 'reduced' survey.design2 object
    design$variables <- NULL
    res$design <- design
    res$call <- match.call()
    res$model$intercept <- dat$intercept
    res$model$yname <- dat$yname
    class(res) <- "svyreg_rob"
    res
}
# robust Tukey biweight GM-estimator of regression (depends on pkg survey)
svyreg_tukeyGM <- function(formula, design, k, type = c("Mallows", "Schweppe"),
    xwgt, var = NULL, na.rm = FALSE, verbose = TRUE, ...)
{
    if (missing(xwgt))
        stop("Argument 'xwgt' is missing\n", call. = FALSE)
    type <- match.arg(type)
    type_int <- switch(type, "Mallows" = 1L, "Schweppe" = 2L)
    dat <- .checkreg(formula, design, var, na.rm)

    if (NCOL(xwgt) > 1) {
        xwgt <- as.numeric(xwgt[, 1])
        warning("Only first column of argument 'xwgt' is used\n", call. = FALSE)
    }
    if (length(xwgt) != length(dat$y))
        stop("Argument 'xwgt' is not of length n\n", call. = FALSE)

    res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, type_int, xwgt, dat$var,
        verbose, ...)
    # add a 'reduced' survey.design2 object
    design$variables <- NULL
    res$design <- design
    res$call <- match.call()
    res$model$intercept <- dat$intercept
    res$model$yname <- dat$yname
    res$model$xwgt <- xwgt
    class(res) <- "svyreg_rob"
    res
}
