# Regression estimator of the mean, Huber M-estimator (depends on pkg survey)
svymean_reg_huber <- function(object, mean_auxiliary, k)
{
    stopifnot(k > 0)
    call <- match.call()
    # check dimensions
    if (length(mean_auxiliary) != object$model$p)
        stop("Dimension of argument 'mean_auxiliary' is not correct\n",
            call. = FALSE)
    # robust greg estimate
    w <- object$model$w; sum_w <- sum(w)
    resid_winsorized <- pmin.int(k, pmax.int(-k, object$residuals))
    est <- sum(mean_auxiliary * object$estimate) + sum(w * resid_winsorized) /
        sum_w
    names(est) <- object$model$yname
    # compute variance
    design <- object$design
    v <- survey::svyrecvar(resid_winsorized * w / sum_w, design$cluster,
        design$strata, design$fpc, postStrata = design$postStrata)
    res <- list(characteristic = "mean", estimator = "GREG M-estiamtor",
        estimate = est, variance = v, design = design, call = call)
    class(res) <- "svystat.rob"
    res
}
