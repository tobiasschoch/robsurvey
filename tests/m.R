library(testthat)
suppressPackageStartupMessages(library(survey))
library(robsurvey, quietly = TRUE)
attach(workplace)

expect_equal(weighted_mean_huber(payroll, weight, k = 1.345), 94838.4395,
    tolerance = 1e-5, label = "Huber M-est")

expect_equal(weighted_mean_tukey(payroll, weight, k = 4.685), 78499.4842,
    tolerance = 1e-5, label = "Tukey M-est")

dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
#-------------------------------------------------------------------------------
# Huber regression M-estimator
expect_equal(
    coef(svyreg_huber(payroll ~ employment, dn, k = 1.345)),
    c("(Intercept)" = 8915.9505, employment = 14214.2246),
    tolerance = 1e-4, label = "Huber regression M-est")

# Tukey regression M-estimator
expect_equal(
    coef(svyreg_tukey(payroll ~ employment, dn, k = 4.685)),
    c("(Intercept)" = 5471.2246, employment =  14467.8600),
    tol = 1e-4, label = "Tukey regression M-est")

#-------------------------------------------------------------------------------
xwgt <- simpsonWgt(employment / weighted_median(employment, weight), a = 1,
    b = 20)

# Mallows regression GM-estimator: Huber psi
expect_equal(
    coef(svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Mallows",
        xwgt = xwgt)),
    c("(Intercept)" = 7843.4909, employment = 14425.3658),
    tolerance = 1e-4, label = "Mallows regression GM-est: Huber")

# Schweppe regression GM-estimator: Huber psi
expect_equal(
    coef(svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Schweppe",
        xwgt = xwgt)),
    c("(Intercept)" = 7769.1899, employment = 14424.6617),
    tolerance = 1e-4, label = "Schweppe regression GM-est: Huber")

expect_equal(1, 2, label = "not a real test")

