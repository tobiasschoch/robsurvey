# Tukey biweight M-estimator of the weighted mean
weighted_mean_tukey <- function(x, w, k, type = "rwm", info = FALSE,
    na.rm = FALSE, verbose = TRUE, ...)
{
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    # select method
    if (type == "rwm") {
        res <- robsvyreg(rep(1, dat$n), dat$x, dat$w, k, 2, 0, NULL, NULL,
            verbose, ...)
    } else if (type == "rht") {
        xvar <- mean(dat$w) / dat$w
        res <- robsvyreg(xvar, dat$x, dat$w, k, 2, 0, NULL, xvar, verbose, ...)
    } else {
        stop(paste0("Method '", type, "' does not exist\n"), call. = FALSE)
    }

    if (length(res) == 1)
        return(NA)

    # check for failure of convergence
    if (!res$optim$converged) {
        res$estimate <- NA
        res$scale <- NA
    }

    if (info) {
        res$model[c("n", "p")] <- NULL
        if (type == "rwm")
            res$model$x <- NULL
        res$characteristic <- "mean"
        res$estimator = paste0("Tukey biweight M-estimator (type = ", type, ")")
        res$robust[c("Epsi2", "Epsiprime")] <- NULL
        res$call <- match.call()
        return(res)
    } else {
        return(res$estimate)
    }
}
# Tukey biweight M-estimator of the weighted total
weighted_total_tukey <- function(x, w, k, type = "rwm", info = FALSE,
    na.rm = FALSE, verbose = TRUE, ...)
{
    res <- weighted_mean_tukey(x, w, k, type, info, na.rm, verbose, ...)
    if (length(res) == 1) {
        return(res * sum(w))
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        return(res)
    }
}
