# Dalén's weight reduction estimating method of the weighted mean
weighted_mean_dalen <- function(x, w, censoring, verbose = TRUE, info = FALSE,
    na.rm = FALSE)
{
    res <- robsurvey::weighted_total_dalen(x, w, censoring, na.rm, verbose,
        info = TRUE)
    res$characteristic <- "mean"
    res$estimate <- res$estimate / sum(res$robust$weightsmod)
    res$call <- match.call()
    if (info)
        return(res)
    else
        return(res$estimate)
}
# Dalén's weight reduction estimating method of the weighted total
weighted_total_dalen <- function(x, w, censoring, verbose = TRUE, info = FALSE,
    na.rm = FALSE)
{
    stopifnot(censoring > 0)
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    xw <- dat$x * dat$w
    if (verbose)
        cat(paste0(sum(xw > censoring), " of ", length(x),
            " observations censored\n"))
    at <- xw > censoring
    if (sum(at) > 0) {
        xw[at] <- censoring + (xw[at] - censoring) / dat$w[at]
        weightsmod <- xw / dat$x
    } else {
        weightsmod <- w
    }
    if (info) {
        res <- list(characteristic <- "total",
            estimator = paste0("Dalen estimator (censored at ", censoring, ")"),
	        estimate = sum(xw), variance = NA,
	        robust = list(censoring = censoring, weightsmod = weightsmod),
	        residuals = NA,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        return(res)
    } else {
        return(sum(xw))
    }
}
