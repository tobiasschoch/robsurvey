# Regression estimator of the mean (depends on pkg survey)
svyreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    dat <- .checkreg(formula, design, var, NULL, na.rm)
    # account for heteroscedasticity
    if (!is.null(var)) {
        dat$x <- dat$x / sqrt(dat$var)
        dat$y <- dat$y / sqrt(dat$var)
    }
    n <- length(dat$y); p <- NCOL(dat$x)
    if (n < p)
        stop("Number of observations cannot be smaller than no. of variables\n")
    tmp <- .C("wlslm", x = as.double(dat$x), y = as.double(dat$y),
        w = as.double(dat$w), resid = as.double(numeric(n)),
        n = as.integer(n), p = as.integer(p), beta = as.double(numeric(p)),
        scale = as.double(numeric(1)), PACKAGE = "robsurvey")
    names(tmp$beta) <- colnames(dat$x)
    # return
    structure(list(characteristic = "regression",
        estimator = list(string = "Weighted least squares", type = 0, psi = 0,
        k = Inf), estimate = tmp$beta, scale = tmp$scale,
        optim = list(converged = TRUE), residuals = tmp$resid,
        model = list(x = dat$x, y = dat$y, w = dat$w, var = dat$var,
        n = n, p = p), design = dat$design, call = match.call()),
        class = "svyreg_rob")
}
