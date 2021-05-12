# weight associated with the Huber psi-function
huberWgt <- function(x, k = 1.345)
{
    stopifnot(k > 0)
    pmin.int(1, k / abs(x))
}
# weight associated with the Tukey biweight psi-function
tukeyWgt <- function(x, k = 4.685)
{
    stopifnot(k > 0)
    (1 - (x / k)^2)^2 * (abs(x) <= k)
}
# weight for regression GM-estiamtors (Simpson, D.G., Ruppert, D., and Carroll,
# R.J. (1992). On One-Step GM Estimates and Stability of Inferences in Linear
# Regression, JASA 87, pp. 439-450).
simpsonWgt <- function(x, a, b)
{
    stopifnot(a >= 0, b >= 0)
    pmin.int(1, (b / x)^(a / 2))
}
