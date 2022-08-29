# ratio estimator of the total
svytotal_ratio <- function(object, total, variance = "wu", keep_object = TRUE)
{
    if (!inherits(object, "ratio"))
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))

    .check_class(object)
    stopifnot(is.numeric(total), total >= 0)

    # estimate
    estimate <- object$estimate * total
    names(estimate) <- object$model$yname

    # standard error
    variance <- match.arg(variance, c("base", "wu", "hajek"))
    stderr <- .stderr_ratio(object, total, variance)

    res <- structure(list(characteristic = "total",
        estimator = object$estimator, estimate = estimate,
        robust = object$robust, residuals = object$residuals,
        model = list(object$model, coef = object$estimate,
            total = total, variance = variance, call = object$call),
        design = object$design, call = match.call(),
        variance = stderr^2), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    res
}
# ratio estimator of the mean
svymean_ratio <- function(object, total, N = NULL, variance = "wu",
    keep_object = TRUE, N_unknown = FALSE)
{
    if (!inherits(object, "ratio"))
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))

    .check_class(object)
    stopifnot(is.numeric(total), total >= 0)
    if (is.null(N)) {
        if (N_unknown)
            N <- sum(object$model$w)
        else
            stop("Argument 'N' is missing\n")
    }
    # estimate
    estimate <- object$estimate * total / N
    names(estimate) <- object$model$yname

    # standard error
    variance <- match.arg(variance, c("base", "wu", "hajek"))
    stderr <- .stderr_ratio(object, total, variance) / N

    res <- structure(list(characteristic = "mean",
        estimator = object$estimator, estimate = estimate,
        robust = object$robust, residuals = object$residuals,
        model = list(object$model, coef = object$estimate,
            total = total, variance = variance, call = object$call),
        design = object$design, call = match.call(),
        variance = stderr^2), class = "svystat_rob")
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
