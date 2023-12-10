# weighted winsorized mean
weighted_mean_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB,
    info = FALSE, na.rm = FALSE)
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

     tmp <- .C(C_wwinsorizedmean, x = as.double(dat$x), w = as.double(dat$w),
        lb = as.double(LB), ub = as.double(UB), loc = as.double(numeric(1)),
        n = as.integer(dat$n))

    if (info) {
        resid <- dat$x - tmp$loc
        res <- list(characteristic = "mean",
	        estimator = list(string = paste0("Weighted winsorized estimator (",
                LB, ", ", UB, ")"), LB = LB, UB = UB),
	        estimate = tmp$loc, variance = NA,
	        residuals = resid,
	        model = list(y = dat$x, w = dat$w),
	        design = NA, call = match.call())
        return(res)
    } else {
        return(tmp$loc)
    }
}
# one-sided weighted k winsorized mean
weighted_mean_k_winsorized <- function(x, w, k, info = FALSE, na.rm = FALSE)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    if (k %% 1 > 0){
        k <- as.integer(k)
        cat(paste0("Argument 'k' is casted to integer: k = ", k,"\n"))
    }
    n <- dat$n
    if (k >= n)
        stop("k must be smaller than n\n", call. = FALSE)
    if (k < 1)
        stop("k must larger than 1\n", call. = FALSE)

    tmp <- .C(C_wkwinsorizedmean, x = as.double(dat$x), w = as.double(dat$w),
        k = as.integer(k), loc = as.double(numeric(1)),
        n = as.integer(n), prob = as.double(numeric(1)))

    if (info) {
        res <- list(characteristic = "mean",
	        estimator = list(string =
                paste0("weighted k winsorized estimator (k = ", k, ")"),
                k = k, UB = tmp$prob),
	        estimate = tmp$loc,
	        variance = NA,
	        residuals = dat$x - tmp$loc,
	        model = list(y = x, w = w),
	        design = NA, call = match.call())
        return(res)
    } else {
        return(tmp$loc)
    }
}
# weighted winsorized total
weighted_total_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB,
    info = FALSE, na.rm = FALSE)
{
    res <- weighted_mean_winsorized(x, w, LB, UB, info, na.rm)
    if (length(res) == 1){
        return(res * sum(w))
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        return(res)
    }
}
# one-sided weighted k winsorized total
weighted_total_k_winsorized <- function(x, w, k, info = FALSE, na.rm = FALSE)
{
    res <- weighted_mean_k_winsorized(x, w, k, info, na.rm)
    if (length(res) == 1){
        return(res * sum(w))
    } else {
        res$characteristic <- "total"
        res$estimate <- res$estimate * sum(w)
        res$call <- match.call()
        return(res)
   }
}
