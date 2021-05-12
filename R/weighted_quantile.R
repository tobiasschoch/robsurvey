# weighted quantile
weighted_quantile <- function(x, w, probs, na.rm = FALSE)
{
    dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
    if (any(probs < 0) || any(probs > 1))
        stop("Argument 'probs' not in [0, 1]\n", call. = FALSE)
    res <- NULL
    for (i in 1:length(probs)) {
        tmp <- .C("wquantile", x = as.double(dat$x), w = as.double(dat$w),
	        n = as.integer(dat$n), probs = as.double(probs[i]),
	        q = as.double(numeric(1)), PACKAGE = "robsurvey")
        res <- c(res, tmp$q)
    }
    names(res) <- paste0(probs * 100, "%")
    return(res)
}
