# Huber M-estimator of the weighted mean
weighted_mean_huber <- function(x, w, k, type = "rhj", asym = FALSE,
    info = FALSE, na.rm = FALSE, verbose = TRUE, ...)
{
    string <- paste0("Huber M-estimator (type = ", type, ifelse(asym,
        "; asym. psi", ""), ")")
    psi <- ifelse(asym, 1, 0)
    dat <- .check(x, w, na.rm)
    # empty data
    if (is.null(dat))
        return(NA)
    # only one observation
    if (dat$n == 1) {
        if (info)
            return(list(characteristic = "mean", estimator = list(string =
                string, type = type, psi = psi, psi_fun = ifelse(asym,
                "asymHuber", "Huber"), k = k), estimate = dat$x, scale = NA,
                residuals = 0, model = list(y = dat$x, w = dat$w),
                call = match.call()))
        else
            return(dat$x)
    }
    # otherwise
    if (type == "rwm") {
        warning("The 'rwm' argument is deprecated; please use 'rhj' instead",
            call. = FALSE)
        type <- "rhj"
    }
    # select method
    if (type == "rhj") {
        res <- robsvyreg(rep(1, dat$n), dat$x, dat$w, k, psi, 0, NULL, NULL,
            verbose, ...)
    } else if (type == "rht") {
        xvar <- mean(dat$w) / dat$w
        res <- robsvyreg(xvar, dat$x, dat$w, k, psi, 0, NULL, xvar, verbose,
            ...)
    } else {
        stop(paste0("Method '", type, "' does not exist\n"), call. = FALSE)
    }

    # check for failure of convergence
    if (!res$optim$converged) {
        res$estimate <- NA
        res$scale <- NA
    }
    # return
    if (info) {
        res$model[c("n", "p")] <- NULL
        if (type == "rhj")
            res$model$x <- NULL
        res$characteristic <- "mean"
        res$estimator$string <- string
        res$estimator$type <- type
        res$call <- match.call()
        return(res)
    } else {
        return(res$estimate)
    }
}
# Huber M-estimator of the weighted total
weighted_total_huber <- function(x, w, k, type = "rhj", asym = FALSE,
    info = FALSE, na.rm = FALSE, verbose = TRUE, ...)
{
    res <- weighted_mean_huber(x, w, k, type, asym, info, na.rm, verbose, ...)
    if (length(res) == 1) {
        res <- res * sum(w)
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
    }
    return(res)
}
