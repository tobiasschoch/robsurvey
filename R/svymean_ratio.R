# ratio estimator of the total
svytotal_ratio <- function(object, total, variance = "wu", keep_object = TRUE)
{
    .check_class(object)
    stopifnot(is.numeric(total), total >= 0)

    # standard error
    variance <- match.arg(variance, c("base", "wu", "hajek"))
    stderr <- .stderr_ratio(object, total, variance)

    res <- structure(list(characteristic = "total",
        estimator = object$estimator, estimate = object$estimate * total,
        robust = object$robust, residuals = object$residuals,
        model = list(object$model, coef = object$estimate,
            total = total, variance = variance),
        design = object$design, call = match.call(),
        variance = stderr^2), class = c("svystat_rob", "ratio_est"))
    if (keep_object)
        res$object <- object
    res
}
# ratio estimator of the mean
svymean_ratio <- function(object, total, N = NULL, variance = "wu",
    keep_object = TRUE, N_unknown = FALSE)
{
    .check_class(object)
    stopifnot(is.numeric(total), total >= 0)
    if (is.null(N)) {
        if (N_unknown)
            N <- sum(object$model$w)
        else
            stop("Argument 'N' is missing\n")
    }
    # standard error
    variance <- match.arg(variance, c("base", "wu", "hajek"))
    stderr <- .stderr_ratio(object, total, variance) / N

    res <- structure(list(characteristic = "mean",
        estimator = object$estimator, estimate = object$estimate * total / N,
        robust = object$robust, residuals = object$residuals,
        model = list(object$model, coef = object$estimate,
            total = total, variance = variance),
        design = object$design, call = match.call(),
        variance = stderr^2), class = c("svystat_rob", "ratio_est"))
    if (keep_object)
        res$object <- object
    res
}
# standard error; estimators of Wu (1982, Biometrika): v0, v1, and v2
.stderr_ratio <- function(object, total, variance)
{
    tx_hat <- sum(object$model$w * object$model$x)
    stderr <- SE.svyreg_rob(object, "design") * tx_hat
    exponent <- switch(variance, "base" = 0, "wu" = 0.5, "hajek" = 1)
    stderr * (total / tx_hat)^exponent
}
# check class
.check_class <- function(object)
{
    if (!inherits(object, "ratio")) {
        if (inherits(object, "svyreg_rob"))
            stop("Object is a regression object, use function with suffix '_reg'\n",
                 call. = FALSE)
        else
            stop(paste0("The function cannot be used for an object of class '",
                class(object), "'\n"), call. = FALSE)
    }
}
