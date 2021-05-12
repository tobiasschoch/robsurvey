mer <- function(object, verbose = TRUE, k = c(1, 50))
{
    if (!inherits(object, "mer_capable"))
        stop("MER-estimator cannot compute for this class of estimators\n")

    # weighted mean or total (reference)
    reference_estimator <- object$call
    reference_estimator$k <- Inf
    ref <- eval(reference_estimator)$estimate

    est <- object$call

    cc <- function(k, est, ref)
    {
        est$k <- k
        tmp <- eval(est)
        tmp$variance + (tmp$estimate - ref)^2
    }
    iter <- 1
    while (1) {
        opt <- stats::optimize(cc, interval = k, est = est, ref = ref)
        if (abs(opt$minimum - k[2]) < 0.01) {
            iter <- iter + 1
            if (0.8 * k[2] > k[1]){
                k <- c(k[1], k[2] * 0.8)
                cat(paste0("range of k is modified to: [", k[1], ", ",
                round(k[2], 1), "]\n"))
            } else {
                cat("\nfailed...\n")
                return(NA)
            }
        } else {
            break;
        }
    }
    if (verbose)
        cat(paste0("\nminimum found for k = ", round(opt$minimum, 4), "\n\n"))

    est$k <- opt$minimum
    result <- eval(est)
    result$estimator <- paste0("MER: ", result$estimator)
    result$call <- match.call()
    result$robust$niter <- iter
    result$robust$k <- opt$minimum
    result
}
