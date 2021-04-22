weighted_mean_dalen <- function(x, w, censored, na.rm = FALSE, verbose = TRUE,
    info = FALSE)
{
    res <- robsurvey::weighted_total_dalen(x, w, censored, na.rm, verbose,
        info = TRUE)
    res$characteristic <- "mean"
    res$estimate <- res$estimate / sum(res$robust$weightsmod)
    res$call <- match.call()
    if (info)
        return(res)
    else
        return(res$estimate)
}

weighted_total_dalen <- function(x, w, censored, na.rm = FALSE, verbose = TRUE,
    info = FALSE)
{
    stopifnot(censored > 0)
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    xw <- dat$x * dat$w
    if (verbose)
        cat(paste0(sum(xw > censored), " of ", length(x),
            " observations censored\n"))
    at <- xw > censored
    if (sum(at) > 0) {
        xw[at] <- censored + (xw[at] - censored) / dat$w[at]
        weightsmod <- xw / dat$x
    }
    if (info) {
        res <- list(characteristic <- "total",
            estimator = paste0("Dalen estimator (censored at ", censored, ")"),
	        estimate = sum(xw), variance = NA,
	        robust = list(censored = censored, weightsmod = weightsmod),
	        residuals = NA,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        return(res)
    } else {
        return(sum(xw))
    }
}
