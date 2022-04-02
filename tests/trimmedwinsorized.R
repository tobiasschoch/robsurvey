library(testthat)
library(robsurvey, quietly = TRUE)
attach(workplace)

#-------------------------------------------------------------------------------
# bare-bone methods
#-------------------------------------------------------------------------------
R_weighted_mean_trimmed <- function(x, w, LB = 0.05, UB = 1 - LB)
{
    if (length(x) == 1)
        return(x)
    q <- weighted_quantile(x, w, c(LB, UB))
    at <- x >= q[1] & x <= q[2]
    sum(x[at] * w[at]) / sum(w[at])
}

R_weighted_mean_winsorized <- function(x, w, LB = 0.05, UB = 1 - LB)
{
    if (length(x) == 1)
        return(x)
    q <- weighted_quantile(x, w, c(LB, UB))

    at <- x < q[1]
    if (sum(at) > 0)
        x[at] <- rep(q[1], sum(at))

    at <- x >= q[2]
    if (sum(at) > 0)
        x[at] <- rep(q[2], sum(at))

    sum(w * x) / sum(w)
}

expect_identical(weighted_mean_trimmed(payroll, weight, LB = 0),
    sum(weight * payroll) / sum(weight))

expect_identical(weighted_mean_trimmed(payroll, weight),
    R_weighted_mean_trimmed(payroll, weight))

expect_identical(weighted_mean_trimmed(payroll, weight, LB = 0, UB = 0.95),
    R_weighted_mean_trimmed(payroll, weight, LB  = 0, UB = 0.95))

expect_identical(weighted_mean_winsorized(payroll, weight, LB = 0),
    sum(weight * payroll) / sum(weight))

expect_identical(weighted_mean_winsorized(payroll, weight),
    R_weighted_mean_winsorized(payroll, weight))

expect_identical(weighted_mean_winsorized(payroll, weight, LB = 0, UB = 0.95),
    R_weighted_mean_winsorized(payroll, weight, LB = 0, UB = 0.95))

expect_identical(weighted_mean_k_winsorized(payroll, weight, k = 2), {
        k <- 2
        ord <- order(payroll)
        x <- payroll[ord]; w <- weight[ord]
        n <- length(x)
        x[(n - k):n] <- rep(x[n - k], k + 1)
        sum(w * x) / sum(w)})

# unweighted
expect_identical(weighted_mean_trimmed(payroll, rep(1, length(payroll)),
        LB = 0.10), mean(payroll, trim = 0.1))

#-------------------------------------------------------------------------------
# survey methods
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(survey))
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# location
expect_identical(unname(coef(svymean_trimmed(~payroll, dn))),
    weighted_mean_trimmed(payroll, weight))

expect_identical(unname(coef(svymean_winsorized(~payroll, dn))),
    weighted_mean_winsorized(payroll, weight))

expect_identical(unname(coef(svymean_k_winsorized(~payroll, dn, k = 2))),
    weighted_mean_k_winsorized(payroll, weight, k = 2))

# variance
expect_identical(unname(vcov(svymean_trimmed(~payroll, dn, LB = 0))),
    unname(vcov(svymean(~payroll, dn))))

expect_identical(unname(vcov(svymean_winsorized(~payroll, dn, LB = 0))),
    unname(vcov(svymean(~payroll, dn))))

# standard error of totals
expect_identical(as.numeric(SE(svytotal_trimmed(~payroll, dn, LB = 0))),
    as.numeric(SE(svytotal(~payroll, dn))))

expect_identical(as.numeric(SE(svytotal_winsorized(~payroll, dn, LB = 0))),
    as.numeric(SE(svytotal(~payroll, dn))))
