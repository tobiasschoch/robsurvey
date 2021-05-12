# print method for robust regression object
print.svyreg_rob <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    converged <- x$optim$converged
    if (is.null(converged) || converged) {
        cat("\n", x$estimator$string, "\n")
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        if (!is.null(x$optim))
            cat("\nIRWLS converged in", x$optim$niter, "iterations\n\n")
        else
            cat("\n")

        if (length(x$estimate)) {
            cat("Coefficients:\n")
            print.default(format(x$estimate, digits = digits), print.gap = 2L,
                quote = FALSE)
        }
        cat("\nScale estimate:", format(signif(x$scale, digits)), "\n")
    } else {
        cat("\nIRWLS not converged in", x$optim$niter, "iterations (with tol =",
            x$optim$tol, ")\n")
    }
    invisible(x)
}
# summary method for robust regression object
summary.svyreg_rob <- function(object, digits = max(3L, getOption("digits")
    - 3L), ...)
{
    converged <- object$optim$converged
    if (is.null(converged) || converged) {
        n <- object$model$n; p <- object$model$p
        r <- object$residuals; N <- sum(object$model$w)

        cat("\nCall:\n")
        print(object$call)
        cat("\nResiduals:\n")
        print(summary(r))
        cat("\nCoefficients:\n")

        # covariance matrix
        tmp <- .C("cov_rwlslm", resid = as.double(r),
            x = as.double(object$model$x), xwgt = as.double(object$model$xwgt),
            robwgt = as.double(object$robust$robweights),
            w = as.double(object$model$w),
            k = as.double(object$estimator$k), scale = as.double(object$scale),
            scale2 = as.double(numeric(1)),
            n = as.integer(n), p = as.integer(p),
            psi = as.integer(object$estimator$psi),
            type = as.integer(object$estimator$type))
        res <- list(stddev = sqrt(tmp$scale2), covmat = matrix(tmp$x[1:(p * p)],
            ncol = p), n = n, p = p, N = N)
        colnames(res$covmat) <- colnames(object$model$x)
        rownames(res$covmat) <- colnames(object$model$x)

        # print out
        est <- object$estimate
        se <- sqrt(diag(res$covmat))
        tval <- est / se
        stats::printCoefmat(cbind(Estimate = est, 'Std. Error' = se,
            't value' = tval, 'Pr(>|t|)' = 2 * stats::pt(abs(tval), N - p,
            lower.tail = FALSE)), digits = digits)

        cat("\nResidual standard error:", format(signif(res$stddev, digits)),
            "on", N - p, "degrees of freedom\n")

        if (!is.null(object$robust))
            cat("\nRobustness weights:\n")
        print(summary(object$robust$robweights))

        invisible(res)
    } else {
        cat("\nIRWLS not converged in", object$optim$niter,
            "iterations (with tol =", object$optim$tol, ")\n")
    }
}
# extract coefficients from robust regression object
coef.svyreg_rob <- function(object, ...)
{
    object$estimate
}
# extract variance from robust regression object
vcov.svyreg_rob <- function(object, ...)
{
    converged <- object$optim$converged
    if (is.null(converged) || converged) {
        n <- object$model$n; p <- object$model$p
        tmp <- .C("cov_rwlslm", resid = as.double(object$residuals),
            x = as.double(object$model$x), xwgt = as.double(object$model$xwgt),
            robwgt = as.double(object$robust$robweights),
            w = as.double(object$model$w), k = as.double(object$estimator$k),
            scale = as.double(object$scale), scale2 = as.double(numeric(1)),
            n = as.integer(n), p = as.integer(p),
            psi = as.integer(object$estimator$psi),
            type = as.integer(object$estimator$type))
        mat <- matrix(tmp$x[1:(p * p)], ncol = p)
        colnames(mat) <- colnames(object$model$x)
        rownames(mat) <- colnames(object$model$x)
        mat
    } else {
        warning("\nnot avaliable (algorithm did not converge)\n")
        NA
    }
}
# extract residuals from robust regression object
residuals.svyreg_rob <- function(object, ...)
{
    object$residuals
}
# extract fitted values from robust regression object
fitted.svyreg_rob <- function(object, ...)
{
    object$model$y - object$residuals
}
# extract robustness weights from robust regression object
robweights.svyreg_rob <- function(object)
{
    tmp <- object$robust$robweights
    if (is.null(tmp)) {
        warning("not available (for this estimator)\n")
        NA
    } else
        tmp
}
# plot method for robust regression object
plot.svyreg_rob <- function (x, which = 1:5, ...)
{
}


# b <- lm(Petal.Width ~ Sepal.Width + Sepal.Length, iris, x = TRUE, y = TRUE)
# m <- list(coefficients = b$coefficients,
#    model = list(x = b$x, y = b$y, w = rep(1, 150), var = rep(1, 150)),
#    robust = list(robweights = rep(1, 150), scale = 1),
#    fitted.values = fitted.values(b),
#    residuals = residuals(b),
#    qr = b$qr,
#    call = "a"
# )
#
# plot.svyreg(m)
#
#
