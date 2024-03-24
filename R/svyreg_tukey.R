# robust Tukey biweight M-estimator of regression (depends on pkg survey)
svyreg_tukeyM <- function(formula, design, k, var = NULL, na.rm = FALSE,
                          verbose = TRUE, ...)
{
    dat <- .check_regression(formula, design, var, NULL, na.rm)
    # add a 'reduced' survey.design2 object
    dat$design$variables <- NULL
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "regression",
            estimator = list(
                string = paste0("Survey regression M-estimator (Tukey psi k = ",
                k, ")"), psi = 2, psi_fun = "Tukey", k = k),
            estimate = rep(NA, NCOL(dat$x)),
            scale = NA, robust = NA, optim = NA, residuals = NA,
            model = list(x = dat$x, y = dat$y, w = dat$w, var = dat$var,
                xwgt = rep(1, length(dat$y)), n = length(dat$y),
                p = NCOL(dat$x), yname = dat$yname),
            design = dat$design, terms = dat$terms, call = match.call()),
            class = "svyreg_rob"))
    # otherwise
    res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, 0, dat$xwgt, dat$var, verbose,
                     ...)
    res$design <- dat$design
    res$terms <- dat$terms
    res$call <- match.call()
    res$model$yname <- dat$yname
    class(res) <- "svyreg_rob"
    res
}
# deprecated function kept for compatibility reasons
svyreg_tukey <- function(formula, design, k, var = NULL, na.rm = FALSE,
                         verbose = TRUE, ...)
{
    .Deprecated("svyreg_tukeyM")
    tmp <- svyreg_tukeyM(formula, design, k, var, na.rm, verbose, ...)
    tmp$call <- match.call()
    tmp
}
# robust Tukey biweight GM-estimator of regression (depends on pkg survey)
svyreg_tukeyGM <- function(formula, design, k, type = c("Mallows", "Schweppe"),
                           xwgt, var = NULL, na.rm = FALSE, verbose = TRUE,
                           ...)
{
    type <- match.arg(type)
    type_int <- switch(type, "Mallows" = 1L, "Schweppe" = 2L)
    if (missing(xwgt))
        stop("Argument 'xwgt' is missing\n", call. = FALSE)
    if (NCOL(xwgt) > 1) {
        xwgt <- as.numeric(xwgt[, 1])
        warning("Only first column of argument 'xwgt' is used\n",
                call. = FALSE)
    }
    dat <- .check_regression(formula, design, var, xwgt, na.rm)
    # add a 'reduced' survey.design2 object
    dat$design$variables <- NULL
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "regression",
            estimator = list(string = paste0("Survey regression ", type,
                " GM-estimator (Tukey psi, k = ", k, ")"),
                psi = 2, psi_fun = "Tukey", k = k),
            estimate = rep(NA, NCOL(dat$x)),
            scale = NA, robust = NA, optim = NA, residuals = NA,
            model = list(x = dat$x, y = dat$y, w = dat$w, var = dat$var,
                xwgt = rep(1, length(dat$y)), n = length(dat$y),
                p = NCOL(dat$x), yname = dat$yname),
            design = dat$design, terms = dat$terms, call = match.call()),
            class = "svyreg_rob"))
    # otherwise
    res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, type_int, dat$xwgt, dat$var,
                     verbose, ...)
    res$design <- dat$design
    res$terms <- dat$terms
    res$call <- match.call()
    res$model$xwgt <- dat$xwgt
    res$model$yname <- dat$yname
    class(res) <- "svyreg_rob"
    res
}
