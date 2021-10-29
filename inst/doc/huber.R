library(MASS)
library(robustbase)
data(iris)



covar <- function(est, cor = TRUE)
{
    k <- est$call$k
    std_res <- as.numeric(residuals(est)) / est$s
    x <- est$x
    p <- NCOL(x)
    n <- NROW(x)
    R <- qr(x)$qr[1L:p, 1L:p, drop = FALSE]
    R[lower.tri(R)] <- 0
    Rinv <- solve(R, diag(p))
    XtX_inv <- tcrossprod(Rinv, Rinv)
    Epsi2 <- sum(huberPsi@psi(std_res, k)^2) / (n - p)
    Epsiprime <- sum(huberPsi@Dpsi(std_res, k)) / n
    varEpsi <- sum((huberPsi@Dpsi(std_res, k) - Epsiprime)^2) / n
    if (cor)
        kappa <- 1 + p / n + varEpsi / Epsiprime^2
    else
        kappa <- 1
    kappa^2 * Epsi2 / Epsiprime^2 * XtX_inv
}

est <- rlm(Sepal.Length ~ Petal.Length + Petal.Width, data = iris,
    method = "M", k = 1)
ref <- lm(Sepal.Length ~ Petal.Length + Petal.Width, data = iris)

covar(est)
vcov(est)

# test covariance
library(testthat)
library(MASS)
library(robsurvey, quietly = TRUE)
library(survey)
data(workplace)
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
# regression estimate
k <- 2
ref <- rlm(payroll ~ employment, weights = weights(dn),
    data = dn$variables, method = "M", wt.method = "case", k = k)
est <- svyreg_huber(payroll ~ employment, dn, k = k,  mad_center = FALSE, tol = 1e-5)
expect_equal(vcov(est, "model"), vcov(ref))
# max. difference in %
100 * max(abs(vcov(est, "model") / vcov(ref) - 1))

est <- svyreg_huberGM(payroll ~ employment, dn, k = k, xwgt = rep(1, 142),
    type = "Mallows", mad_center = FALSE, tol = 1e-5)

