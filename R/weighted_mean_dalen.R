# Dalén's weight reduction estimating method of the weighted mean
weighted_mean_dalen <- function(x, w, censoring, type = "Z2", info = FALSE,
    na.rm = FALSE, verbose = TRUE)
{
    res <- weighted_total_dalen(x, w, censoring, type, na.rm, verbose,
        info = TRUE)
    if (is.na(res$estimate))
        return(NA)

    res$characteristic <- "mean"
    res$estimate <- res$estimate / sum(res$model$w)
    res$call <- match.call()
    if (info)
        res
    else
        res$estimate
}
# Dalén's weight reduction estimating method of the weighted total
weighted_total_dalen <- function(x, w, censoring, type = "Z2", info = FALSE,
    na.rm = FALSE, verbose = TRUE)
{
    stopifnot(censoring > 0)
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    xw <- dat$x * dat$w
    at <- xw > censoring
    n_censored <- sum(at)
    if (n_censored > 0) {
        if (type == "Z2") {             # Z2 estimator
            xw[at] <- censoring
            if (verbose)
                cat(paste0(n_censored, " of ", length(dat$w),
                           " observations censored\n"))
        } else if (type == "Z3") {      # Z3 estimator
            xw[at] <- censoring + (dat$x[at] - censoring / dat$w[at])
        } else {
            stop("Argument 'type' must either 'Z2' or 'Z3'\n", call. = FALSE)
        }
    } else if (verbose) {
        cat("No observations have been censored\n")
    }
    estimate <- sum(xw)
    if (info) {
        res <- list(characteristic = "total",
            estimator = list(string = paste0("Dalen ", type,
                " estimator (censored at ", censoring, ")"),
                censoring = censoring),
	        estimate = estimate, variance = NA,
	        robust = list(xw = xw),
	        residuals = NA,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        res
    } else {
        estimate
    }
}
