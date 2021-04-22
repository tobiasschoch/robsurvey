weighted_mean_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
    na.rm = FALSE)
{
    dat <- .check(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    if (LB >= UB)
        stop("Argument 'LB' must be smaller than 'UB'!", call. = FALSE)
    if (LB < 0)
        stop("Argument 'LB' must not be < 0!", call. = FALSE)
    if (UB > 1)
        stop("Argument 'UB' must not be > 1!", call. = FALSE)

    tmp <- .C("wtrimmedmean", x = as.double(dat$x), w = as.double(dat$w),
        lb = as.double(LB), ub = as.double(UB), loc = as.double(numeric(1)),
        n = as.integer(dat$n), PACKAGE = "robsurvey")
     if (info) {
        resid <- dat$x - tmp$loc
        res <- list(characteristic = "mean",
	        estimator = paste0("Weighted trimmed estimator (", LB, ", ",
                UB, ")"),
	        estimate = tmp$loc, variance = NA,
	        robust = list(UB = UB, LB = LB),
	        residuals = resid,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        return(res)
    } else {
        return(tmp$loc)
    }
}

weighted_total_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
    na.rm = FALSE)
{
    res <- weighted_mean_trimmed(x, w, LB, UB, info, na.rm)
    if (length(res) == 1) {
        return(res * sum(w))
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        return(res)
    }
}
