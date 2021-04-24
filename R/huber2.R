huber2 <- function(x, w, k = 1.5, na.rm = FALSE, maxit = 50, tol = 1e-4,
    info = FALSE, k_Inf = 1e5, df_cor = TRUE)
{
    stopifnot(k_Inf > 0)
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    kk <- ifelse(is.finite(k), k, k_Inf)

    tmp <- .C("huberm", x = as.double(dat$x), w = as.double(dat$w),
        robwgt = as.double(numeric(dat$n)), k = as.double(kk),
        loc = as.double(numeric(1)), scale = as.double(numeric(1)),
        n = as.integer(dat$n), maxit = as.integer(maxit), tol = as.double(tol),
        df_cor = as.integer(df_cor), success = as.integer(0),
        PACKAGE = "robsurvey")

    if (tmp$success == 0) {
        warning("Initial estimate of scale (IQR) is zero\n", call. = FALSE)
        return(NA)
    }

    if (tmp$maxit == maxit)
        warning(paste0("Failure of convergence\n"), call. = FALSE)

    if (info) {
        res <- list(characteristic = "mean",
	        estimator = paste0("Weighted Huber proposal 2 estimator (k=", k,
                ")"),
            estimate = tmp$loc, variance = NA, scale = tmp$scale,
            robust = list(k = k, robweights = tmp$robwgt),
            optim = list(converged = (tmp$maxit < maxit), niter = tmp$maxit,
                tol = tmp$tol),
            residuals = dat$x - tmp$loc,
            model = list(y = dat$x, w = dat$w), design = NA,
            call = match.call())
        return(res)
    } else {
        return(tmp$loc)
    }
}
