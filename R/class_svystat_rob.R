#FIXME: remove the paste0 in cat


summary.svystat_rob <- function(object, digits = max(3L, getOption("digits") -
    3L), ...)
{
    cat(paste0(object$estimator, " of the sample ", object$characteristic,
        "\n"))
    cat("\n")
    est <- cbind(object$estimate, sqrt(object$variance))
    colnames(est) <- c(object$characteristic, "SE")
    print(est, digits)
    cat("\n")
    if (!is.null(object$optim)) {
        cat("Robustness:\n")
        cat(paste0("  Psi-function: ", object$robust$psifunction, " with k = ",
	        object$robust$k, "\n"))
        cat(paste0("  mean of robustness weights: ",
	        round(mean(object$robust$robweights), digits), "\n"))
        cat("\n")
        cat("Algorithm performance:\n")
        if (object$optim$converged) {
            cat(paste0("  converged in ", object$optim$niter, " iterations \n"))
	        cat(paste0("  with residual scale (weighted MAD): ",
	            format(object$robust$scale, digits = digits), "\n"))
        } else {
	        cat(paste0("  FAILURE of convergence in ", object$optim$niter,
	            " iterations \n"))
	        cat(paste0("  with residual scale (weighted MAD): ",
	            round(object$robust$scale, digits), "\n"))
        }
    }
#FIXME:
    #cat("Sampling design:\n")
    #print(object$design)
}

coef.svystat_rob <- function(object, ...)
{
    object$estimate
}

SE.svystat_rob <- function(object, ...)
{
    sqrt(object$variance)
}

vcov.svystat_rob <- function(object, ...)
{
    v <- as.matrix(object$variance)
    rownames(v) <- names(object$estimate)
    colnames(v) <- "Variance"
    v
}

residuals.svystat_rob <- function(object, ...)
{
    object$residuals
}

fitted.svystat_rob <- function(object, ...)
{
    object$model$y - object$residuals
}

robweights <- function(object)
{
    UseMethod("robweights", object)
}

robweights.svystat_rob <- function(object)
{
    tmp <- object$robust$robweights
    if (is.null(tmp))
        stop("Robustness weights are not available\n")
    else
        tmp
}

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
        cat(paste0(x$call[[1]], ": failure of convergence in ", x$optim$niter,
            " steps\n"))
        cat("(you may use the 'summary' method to see more details)\n")
    }
}