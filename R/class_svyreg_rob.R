# print method for robust regression object
print.svyreg_rob <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    if (any(is.na(x$estimate))) {
        cat(paste0("\n", x$estimator$string, "\n"))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        cat("\nCoefficients:\n")
        print.default(format(x$estimate, digits = digits), print.gap = 2L,
            quote = FALSE)
        cat("\nScale estimate:", NA, "\n")
        return(invisible(x))
    }
    converged <- x$optim$converged
    if (is.null(converged) || converged) {
        cat(paste0("\n", x$estimator$string, "\n"))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        if (is.finite(x$estimator$k))
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
summary.svyreg_rob <- function(object, mode = c("design", "model", "compound"),
    digits = max(3L, getOption("digits") - 3L), ...)
{
    if (any(is.na(object$estimate))) {
        cat("\nCall:\n")
        print(object$call)
        cat("\nResiduals:\n")
        print(NA)
        cat("\nCoefficients:\n")
        NAs <- object$estimate
        stats::printCoefmat(cbind(Estimate = NAs, 'Std. Error' = NAs,
            't value' = NAs, 'Pr(>|t|)' = NAs), digits = digits)
        cat("\nResidual standard error:", NA, "on", NA, "degrees of freedom\n")
        return(invisible(object))
    }
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
        tmp <- switch(match.arg(mode),
            "model" = .cov_reg_model(object),
            "design" = .cov_reg_design(object),
            "compound" = .cov_reg_compound(object))
        if (tmp$ok == 0)
            stop("\nCovariance estimation failed\n", call. = FALSE)

        res <- list(stddev = tmp$scale, covmat = tmp$cov, n = n, p = p, N = N)
        ns <- colnames(object$model$x)
        dimnames(res$covmat) <- list(ns, ns)

        # print out
        est <- object$estimate
        se <- sqrt(diag(res$covmat))
        tval <- est / se
        stats::printCoefmat(cbind(Estimate = est, 'Std. Error' = se,
            't value' = tval, 'Pr(>|t|)' = 2 * stats::pt(abs(tval), N - p,
            lower.tail = FALSE)), digits = digits)

        cat("\nResidual standard error:", format(signif(res$stddev, digits)),
            "on", N - p, "degrees of freedom\n")

        if (is.finite(object$estimator$k)) {
            cat("\nRobustness weights:\n")
            print(summary(object$robust$robweights))
        }
        invisible(res)
    } else {
        warning(" covariance is not avaliable because",
            "\n regression algorithm did not converge\n", call. = FALSE)
    }
}
# extract variance from robust regression object
vcov.svyreg_rob <- function(object, mode = c("design", "model", "compound"),
    ...)
{
    converged <- object$optim$converged
    if (is.null(converged) || converged) {
        mat <- switch(match.arg(mode),
            "model" = .cov_reg_model(object)$cov,
            "design" = .cov_reg_design(object)$cov,
            "compound" = .cov_reg_compound(object)$cov)
    } else {
        warning(" covariance is not avaliable because",
            "\n regression algorithm did not converge\n", call. = FALSE)
        mat <- matrix(NA, object$model$p, object$model$p)
    }
    ns <- colnames(object$model$x)
    dimnames(mat) <- list(ns, ns)
    mat
}

# model-based covariance matrix of M- and GM-regression estimators
.cov_reg_model <- function(object)
{
    # robustness tuning constant
    k <- object$estimator$k
    if (is.infinite(k))
        k <- svyreg_control()$k

    # account for heteroscedasticity
    r <- object$residuals
    v <- object$model$var
    if (!is.null(v))
        r <- r / sqrt(v)

    tmp <- .C("cov_reg_model", resid = as.double(r),
        x = as.double(object$model$x), xwgt = as.double(object$model$xwgt),
        robwgt = as.double(object$robust$robweights),
        w = as.double(object$model$w), k = as.double(k),
        scale = as.double(object$scale), scale = as.double(numeric(1)),
        n = as.integer(object$model$n), p = as.integer(object$model$p),
        psi = as.integer(object$estimator$psi),
        type = as.integer(object$estimator$type), ok = as.integer(0),
        PACKAGE = "robsurvey")
    p <- object$model$p
    if (tmp$ok == 0) {
        warning("Covariance estimation failed", call. = FALSE)
        cov_mat <- matrix(NA, p, p)
    } else {
        cov_mat <- matrix(tmp$x[1:(p * p)], ncol = p)
    }
    list(ok = tmp$ok, cov = cov_mat, scale = tmp$scale)
}

# design-based covariance matrix of M- and GM-regression estimators
.cov_reg_design <- function(object)
{
    # robustness tuning constant
    k <- object$estimator$k
    if (is.infinite(k))
        k <- svyreg_control()$k

    x <- object$model$x
    r <- object$residuals
    n <- NROW(x); p <- NCOL(x)

    # account for heteroscedasticity
    v <- object$model$var
    if (!is.null(v))
        r <- r / sqrt(v)

    # weights in the model's design space
    xwgt <- object$model$xwgt
    if (is.null(xwgt))
        xwgt <- rep(1, n)

    # account for heteroscedasticity
    if (!is.null(object$model$var))
        x <- x / sqrt(object$model$var)

    # scale the weights (prevent overflow)
    w <- object$model$w / sum(object$model$w)

    # Q matrix
    Q <- survey::svyrecvar(w * xwgt * x * object$robust$robweights * r,
        object$design$cluster, object$design$strata, object$design$fpc)

    # covariance matrix (takes Q matrix as one of its argument)
    tmp <- .C("cov_reg_design", x = as.double(x), w = as.double(w),
        xwgt = as.double(xwgt), resid = as.double(r),
        scale = as.double(object$scale), k = as.double(k),
        psi = as.integer(object$estimator$psi),
        type = as.integer(object$estimator$type), n = as.integer(n),
        p = as.integer(p), ok = as.integer(0), mat = as.double(Q),
        PACKAGE = "robsurvey")

    if (tmp$ok == 0) {
        warning("Covariance estimation failed", call. = FALSE)
        cov_mat <- matrix(NA, p, p)
    } else {
        cov_mat <- matrix(tmp$mat, ncol = p, nrow = p)
    }
    list(ok = 1, cov = cov_mat, scale = object$scale)
}

# compound design-model-based covariance matrix of M- and GM-regression est.
.cov_reg_compound <- function(object)
{
    gamma_m <- .cov_reg_model(object)
    gamma_d <- .cov_reg_design(object)
    # sampling fraction
    f <- object$model$n / sum(object$model$w)
    if (gamma_d$ok * gamma_m$ok == 0)
        cov_mat <- matrix(NA, object$model$p, object$model$p)
    else
        cov_mat <- gamma_d$cov + f * gamma_m$cov
    list(ok = 1, cov = cov_mat, scale = object$scale)
}

# extract coefficients from robust regression object
coef.svyreg_rob <- function(object, ...)
{
    object$estimate
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
    } else {
        tmp
    }
}
# plot method for robust regression object
plot.svyreg_rob <- function (x, which = 1:5, ...)
{
    .NotYetImplemented()
}
