# weighted trimmed mean
weighted_mean_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
    na.rm = FALSE)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)

    if (LB >= UB)
        stop("Argument 'LB' must be smaller than 'UB'!", call. = FALSE)
    if (LB < 0)
        stop("Argument 'LB' must not be < 0!", call. = FALSE)
    if (UB > 1)
        stop("Argument 'UB' must not be > 1!", call. = FALSE)

    tmp <- .C(C_wtrimmedmean, x = as.double(dat$x), w = as.double(dat$w),
        lb = as.double(LB), ub = as.double(UB), loc = as.double(numeric(1)),
        n = as.integer(dat$n), success = as.integer(1))

    if (tmp$success == 0) {
        warning("Division by zero\n", call. = FALSE)
        return(NA)
    }

    if (info) {
        resid <- dat$x - tmp$loc
        res <- list(characteristic = "mean",
	        estimator = list(string = paste0("Weighted trimmed estimator (",
                LB, ", ", UB, ")"), LB = LB, UB = UB),
	        estimate = tmp$loc, variance = NA,
	        residuals = resid,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        res
    } else {
        tmp$loc
    }
}
# weighted trimmed total
weighted_total_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
    na.rm = FALSE)
{
    res <- weighted_mean_trimmed(x, w, LB, UB, info, na.rm)
    if (length(res) == 1) {
        res * sum(w)
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        res
    }
}
