# Minimum estimated risk estimator for location M-estimates
mer <- function(object, verbose = TRUE, max_k = 10, init = 1, method = "Brent",
    optim_args = list())
{
    if (!inherits(object, "mer_capable"))
        stop("MER-estimator cannot be compute for this class of estimators\n",
            call. = FALSE)
    stopifnot(is.numeric(max_k), max_k > 1, is.numeric(init), init >= 0,
        init < max_k)

    # search interval
    range_k <- c(init, max_k)

    # weighted mean or total (reference estimator)
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
    opt <- optim(init, estimated_mse, est = est, ref = ref, method = method,
        lower = range_k[1], upper = range_k[2], control = optim_args)

    if (verbose)
        cat(paste0("Search interval: [", range_k[1], ", ",
            round(range_k[2], 1), "]\n"))

    if (opt$par < 1) {
        cat("No minimum found\n")
        result <- ref
        result$estimate <- NA
        result$scale <- NA
        result[c("robust", "residuals")] <- NULL
        result$optim <- list(converged = FALSE, niter = 0)
    } else {
        if (verbose) {
            cat("Minimum found for k = ", round(opt$par, 4), "\n")
            cat(paste0("Rel. efficiency gain: ", -100 * round(opt$val - 1, 2),
                "%\n"))
        }
        # compute mer
        est$k <- round(opt$par, 4)
        result <- eval(est)
        result$robust$rel_mse <- opt$value - 1
    }
    result$estimator$string <- paste("MER:", result$estimator$string)
    result$call <- match.call()
    result
}
