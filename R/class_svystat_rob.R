# Slots of objects of class 'svystat_rob'
#  + characteristic: "regression"
#  + estimator [list]
#     + string: [char]
#     + type: [int]
#     + psi: [int]
#     + psi_fun: [char]
#     + k: [numeric]
#  + estimate: [numeric] vector of estimated regression coefficients
#  + scale: [numeric] scale estimate
#  + robust [list]
#     + robweights: [numeric] robustness weights
#     + outliers: [numeric] indicator variable
#  + optim [list]
#     + converged: [logical]
#     + niter: [int] number of IRWLS iterations
#     + tol: [numeric] numerical tolerance criterion (IRLWS)
#     + used_iqr: [int] 1 = scale estimated by IQR not MAD
#  + residuals: [numeric]
#  + model [list]
#     + x: [matrix] design matrix (only for GREG type estimators)
#     + y: [numeric] response variable
#     + w: [numeric] sampling weights
#     + var: [numeric] heteroscedasticity variances
#     + xwgt: [numeric] weights in the model's design space (GM-estimator)
#     + n [int] number of observations
#     + p [int] number of independent variables
#     + [others]
#  + design: [survey.design object without 'variables']
#  + variance: [numeric]
#  + call: [call object]
#
# constructor for empty object of class svystat_rob
.new_svystat_rob <- function(characteristic, yname, string, call,
    design, class = NULL, ...)
{
    structure(list(characteristic = characteristic,
        estimator = list(string = string, ...),
        estimate = stats::setNames(NA, yname), variance = NA,
        residuals = NA, model = NA, design = design, call = call),
        class = c("svystat_rob", class))
}
# summary method for robust survey statistic object
summary.svystat_rob <- function(object, digits = max(3L, getOption("digits") -
    3L), ...)
{
    cat(object$estimator$string, "of the population", object$characteristic,
        "\n")
    cat("\n")
    est <- cbind(object$estimate, sqrt(object$variance))
    colnames(est) <- c(object$characteristic, "SE")
    print(est, digits)
    cat("\n")
    if (!is.null(object$optim)) {
        cat("Robustness:\n")
        cat("  Psi-function:", object$robust$psifunction, "with k =",
	        object$estimator$k, "\n")
        cat("  mean of robustness weights:",
	        round(mean(object$robust$robweights), digits), "\n")
        cat("\n")
        cat("Algorithm performance:\n")
        if (object$optim$converged) {
            cat("  converged in", object$optim$niter, "iterations\n")
	        cat(paste0("  with residual scale ", ifelse(object$optim$used_iqr,
                "(weighted IQR): ", "(weighted MAD): "),
	            format(object$scale, digits = digits), "\n\n"))
        } else {
	        cat("  FAILURE of convergence in", object$optim$niter,
	            " iterations\n\n")
        }
    }
    cat("Sampling design:\n")
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
        cat(x$call[[1]], ": failure of convergence\n")
        cat("(you may use the 'summary' method to learn more)\n")
    }
}

# extract estimate of scale
scale.svystat_rob <- function(x, ...)
{
    x$scale
}
# compute estimated mse, more precisely, estimated risk
mse <- function(object)
{
    call <- object$call
    # reference estimator
    if (grepl("_reg$", call[[1]])) {                # robust GREG
        reg_call <- object$model$call
        tmp <- eval(as.call(list(substitute(svyreg), reg_call$formula,
            reg_call$design, reg_call$var, !is.null(reg_call$na.rm))))
        call$object <- substitute(tmp)
        call$k <- NULL
        call$type <- "ADU"
        ref <- coef.svystat_rob(eval(call))
    } else if (grepl("_ratio$", call[[1]])) {       # robust ratio
        rat_call <- object$model$call
        rat_call$k <- Inf
        rat_call$asym <- FALSE
        rat_call$verbose <- FALSE
        tmp <- eval(rat_call)
        call$object <- substitute(tmp)
        ref <- coef.svystat_rob(eval(call))
    } else {                                        # otherwise
        ref <- weighted_total(object$model$y, object$model$w,
            na.rm = !is.null(call$na.rm))
        if (object$characteristic == "mean")
            ref <- ref / sum(object$model$w)
    }
    # mse
    as.numeric(object$variance + (object$estimate - ref)^2)
}
