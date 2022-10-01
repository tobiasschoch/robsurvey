# Tukey biweight M-estimator of the weighted mean
weighted_mean_tukey <- function(x, w, k, type = "rwm", info = FALSE,
    na.rm = FALSE, verbose = TRUE, ...)
{
    type <- match.arg(type, c("rht", "rhj", "rwm"))
    if (type == "rhj")
        type <- "rwm"

    string <- paste0("Tukey M-estimator (type = ", type, ")")
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    # only one observation
    if (dat$n == 0)
        return(dat$x)
    # only one observation
    if (dat$n == 1) {
        if (info)
            return(list(characteristic = "mean", estimator = list(string =
                string, type = type, psi = 2, psi_fun = "Tukey", k = k),
                estimate = dat$x, scale = NA, residuals = 0, model =
                list(y = dat$x, w = dat$w), call = match.call()))
        else
            return(dat$x)
    }
    # otherwise
    if (type == "rwm") {
        res <- robsvyreg(rep(1, dat$n), dat$x, dat$w, k, 2, 0, NULL, NULL,
            verbose, ...)
    }
    if (type == "rht") {
        xvar <- mean(dat$w) / dat$w
        res <- robsvyreg(xvar, dat$x, dat$w, k, 2, 0, NULL, xvar, verbose, ...)
    }

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
        res$estimator$string <- string
        res$estimator$type <- type
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
