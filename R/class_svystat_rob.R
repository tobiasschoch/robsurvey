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
#     + x: [matrix] design matrix
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
    cat(object$estimator$string, "of the sample", object$characteristic,
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
# compute estimated mse, more precisely, estimated risk; see mer()
mse <- function(object)
{
    if (!inherits(object, "svystat_rob"))
        stop("MSE is defined only for object of class 'svystat_rob'\n",
            call. = FALSE)

    # consistent reference estimator (mean or total)
    reference_estimator <- object$call
    if (inherits(object, "mest"))
        reference_estimator$k <- Inf
    if (inherits(object, "dalen"))
        reference_estimator$censoring <- Inf
    if (inherits(object, c("wins", "trim"))) {
        reference_estimator$LB<- 0
        reference_estimator$UB<- 1
    } else
        reference_estimator$verbose <- FALSE

    reference_location <- coef.svystat_rob(eval(reference_estimator))
    # estimated mse
    as.numeric(object$variance + (object$estimate - reference_location)^2)
}
