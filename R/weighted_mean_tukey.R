weighted_mean_tukey <- function(x, w, k = 1.5, type = "rwm", info = FALSE,
    na.rm = FALSE, ...)
{
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    res <- if (type == "rwm")
        robsvyreg(rep(1, dat$n), dat$x, dat$w, k, 2, 0, NULL, NULL, ...)
    else if (type == "rht")
        robsvyreg(mean(dat$w) / dat$w, dat$x, dat$w, k, 2, 0, NULL, x, ...)
    else
        stop(paste0("Method '", type, "' does not exist\n"), call. = FALSE)

    if (length(res) == 1)
        return(NA)

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

weighted_total_tukey <- function(x, w, k = 1.5, type = "rwm", info = FALSE,
    na.rm = FALSE, ...)
{
    res <- weighted_mean_tukey(x, w, k, type, info, na.rm, ...)
    if (length(res) == 1) {
        return(res * sum(w))
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        return(res)
    }
}
