# weight associated with the Huber psi-function
huberWgt <- function(x, k)
{
    stopifnot(k > 0)
    pmin.int(1, k / abs(x))
}

# weight associated with the Tukey biweight psi-function
tukeyWgt <- function(x, k)
{
    stopifnot(k > 0)
    (1 - (x / k)^2)^2 * (abs(x) <= k)
}

# weight for regression GM-estiamtors (Simpson, D.G., Ruppert, D., and Carroll,
# R.J. (1992). On One-Step GM Estimates and Stability of Inferences in Linear
# Regression, JASA 87, pp. 439-450).
simpsonWgt <- function(object, ...)
{
    UseMethod("simpsonWgt", object)
}

simpsonWgt.default <- function(object, a = 1, b)
{
    stopifnot(a >= 0, b >= 0)
    pmin.int(1, (b / object)^(a / 2))
}

simpsonWgt.robmv <- function(object, a = 1, b = NULL)
{
    b <- if (is.null(b))
        object$cutoff
    simpsonWgt.default(distance(object), a, b)
}
