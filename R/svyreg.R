# Regression estimator of the mean (depends on pkg survey)
svyreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    dat <- .check_regression(formula, design, var, NULL, na.rm)
    n <- length(dat$y); p <- NCOL(dat$x)
    # account for heteroscedasticity
    if (is.null(var)) {
        x <- dat$x
        y <- dat$y
    } else {
        x <- dat$x / sqrt(dat$var)
        y <- dat$y / sqrt(dat$var)
    }
    if (n < p)
        stop("Number of observations cannot be smaller than no. of variables\n")
    tmp <- .C(C_wlslm, x = as.double(x), y = as.double(y),
        w = as.double(dat$w), resid = as.double(numeric(n)),
        n = as.integer(n), p = as.integer(p), beta = as.double(numeric(p)),
        scale = as.double(numeric(1)))
    # Note: The residuals are (y - x*beta) / sqrt(var)
    names(tmp$beta) <- colnames(dat$x)
    # return
    structure(list(characteristic = "regression",
        estimator = list(string = "Weighted least squares", type = 0, psi = 0,
        k = Inf), estimate = tmp$beta, scale = tmp$scale,
        optim = list(converged = TRUE), residuals = tmp$resid,
        model = list(x = dat$x, y = dat$y, w = dat$w, var = dat$var,
        n = n, p = p, yname = dat$yname), design = dat$design,
        terms = dat$terms, call = match.call()), class = "svyreg_rob")
}
