# Minimum estimated risk estimator for location M-estimates
mer <- function(object, verbose = TRUE, max_k = 1000, optim_args = list())
{
    if (!inherits(object, "mer_capable"))
        stop("MER-estimator cannot be compute for this class of estimators\n",
            call. = FALSE)
    stopifnot(is.numeric(max_k) && max_k > 1)

    # search interval
    range_k <- c(0.5, max_k)

    # weighted mean or total (reference)
    reference_estimator <- object$call
    reference_estimator$k <- Inf
    # estimate of the design-unbiased estimator (mean or total)
    ref <- eval(reference_estimator)

    est <- object$call

    # function to be minimized (estimated mean square error)
    estimated_mse <- function(k, est, ref)
    {
        ref_location <- ref$estimate
        est$k <- k
        tmp <- eval(est)
        (tmp$variance + (tmp$estimate - ref_location)^2) / ref$variance
    }

    # minimize mse
    opt <- stats::optim(1, estimated_mse, est = est, ref = ref,
        method = "L-BFGS-B", lower = range_k[1], upper = range_k[2],
        control = optim_args)

    if (verbose)
        cat(paste0("Search interval: [", range_k[1], ", ", round(range_k[2], 1),
            "]\nMinimum found for k = ", round(opt$par, 4), "\n"))
    if (any(abs(opt$par - range_k) < .Machine$double.eps^0.25))
        cat("Minimum is attained at the boundary of the search interval!\n")
    cat("\n")
    est$k <- round(opt$par, 4)
    result <- eval(est)
    result$estimator <- paste("MER:", result$estimator)
    result$call <- match.call()
    result
}
