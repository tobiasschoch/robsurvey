# summary method for robust survey statistic object
summary.svystat_rob <- function(object, digits = max(3L, getOption("digits") -
    3L), ...)
{
    cat(object$estimator, " of the sample ", object$characteristic,
        "\n")
    cat("\n")
    est <- cbind(object$estimate, sqrt(object$variance))
    colnames(est) <- c(object$characteristic, "SE")
    print(est, digits)
    cat("\n")
    if (!is.null(object$optim)) {
        cat("Robustness:\n")
        cat("  Psi-function:", object$robust$psifunction, "with k =",
	        object$robust$k, "\n")
        cat("  mean of robustness weights:",
	        round(mean(object$robust$robweights), digits), "\n")
        cat("\n")
        cat("Algorithm performance:\n")
        if (object$optim$converged) {
            cat("  converged in", object$optim$niter, "iterations\n")
	        cat("  with residual scale (weighted MAD):",
	            format(object$scale, digits = digits), "\n")
        } else {
	        cat("  FAILURE of convergence in", object$optim$niter,
	            " iterations\n")
        }
    }
    cat("\nSampling design:\n")
    print(object$design)
}
# extract estimate from robust survey statistic object
coef.svystat_rob <- function(object, ...)
{
    object$estimate
}
# extract standard error from robust survey statistic object
SE.svystat_rob <- function(object, ...)
{
    sqrt(object$variance)
}
# extract variance from robust survey statistic object
vcov.svystat_rob <- function(object, ...)
{
    v <- as.matrix(object$variance)
    rownames(v) <- names(object$estimate)
    colnames(v) <- "Variance"
    v
}
# extract residuals from robust survey statistic object
residuals.svystat_rob <- function(object, ...)
{
    object$residuals
}
# extract fitted values from robust survey statistic object
fitted.svystat_rob <- function(object, ...)
{
    object$model$y - object$residuals
}
# extract robustness weights from robust survey statistic object, generic
robweights <- function(object)
{
    UseMethod("robweights", object)
}
# extract robustness weights from robust survey statistic object
robweights.svystat_rob <- function(object)
{
    tmp <- object$robust$robweights
    if (is.null(tmp))
        stop("Robustness weights are not available\n")
    else
        tmp
}
# print method for robust survey statistic object
print.svystat_rob <- function(x, digits = max(3L, getOption("digits") - 3L),
    ...)
{
    conv <- TRUE
    if (!is.null(x$optim))
        conv <- x$optim$converged

    if (conv) {
        m <- cbind(x$estimate, sqrt(x$variance))
        colnames(m) <- c(x$characteristic, "SE")
        print(m, digits)
    } else {
        cat(x$call[[1]], ": failure of convergence in", x$optim$niter,
            "steps\n")
        cat("(you may use the 'summary' method to learn more)\n")
    }
}

# extract estimate of scale
scale.svystat_rob <- function(x, ...)
{
    x$scale
}
