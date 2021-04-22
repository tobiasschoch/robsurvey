huber2 <- function(x, w, k, na.rm = FALSE, maxit = 50, tol = 1e-4,
    info = FALSE)
{
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    kk <- if (is.finite(k))
        k
    else
        1e5

    tmp <- .C("huberm", x = as.double(x), w = as.double(w),
        robwgt = as.double(numeric(dat$n)), k = as.double(kk),
        loc = as.double(numeric(1)), scale = as.double(numeric(1)),
        n = as.integer(dat$n), maxit = as.integer(maxit), tol = as.double(tol),
        PACKAGE = "robsurvey")
    if (info) {
        resid <- dat$x - tmp$loc
        res <- list(characteristic = "mean",
	        estimator = paste0("Weighted Huber proposal 2 estimator (k=", k,
                ")"),
            estimate = tmp$loc, variance = NA,
            robust = list(k = k, robweights = tmp$robwgt),
                residuals = tmp$resid,
            model = list(y = dat$x, w = dat$w), design = NA,
            call = match.call())
        return(res)
    } else {
        return(tmp$loc)
    }
}
