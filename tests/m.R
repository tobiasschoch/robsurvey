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
est <- svyreg_huber(payroll ~ employment, dn, k = 1.345)
expect_equal(coef(est),
    c("(Intercept)" = 8915.9505, employment = 14214.2246),
    tolerance = 1e-4, label = "Huber regression M-est")
# covariance matrix
ref <- matrix(c(36075.4946, -540.54394, -540.5439, 49.04239), ncol = 2)
colnames(ref) <- c("(Intercept)", "employment")
rownames(ref) <- c("(Intercept)", "employment")
expect_equal(vcov(est), ref, tolerance = 1e-4,
    label = "Huber regression M-est: cov")

#-------------------------------------------------------------------------------
# Tukey regression M-estimator
est <- svyreg_tukey(payroll ~ employment, dn, k = 4.685)
expect_equal(coef(est),
    c("(Intercept)" = 5471.2246, employment =  14467.8600),
    tol = 1e-4, label = "Tukey regression M-est")
# covariance matrix
ref <- matrix(c(27846.8367, -417.24830, -417.24830, 37.85604), ncol = 2)
colnames(ref) <- c("(Intercept)", "employment")
rownames(ref) <- c("(Intercept)", "employment")
expect_equal(vcov(est), ref, tolerance = 1e-4,
    label = "Tukey regression M-est: cov")

#-------------------------------------------------------------------------------
# Mallows regression GM-estimator: Huber psi
xwgt <- simpsonWgt(employment / weighted_median(employment, weight), a = 1,
    b = 20)
est <- svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Mallows",
    xwgt = xwgt)
expect_equal(coef(est),
    c("(Intercept)" = 7843.4909, employment = 14425.3658),
    tolerance = 1e-4, label = "Huber Mallows regression GM-est")
# covariance matrix
#FIXME:
#vcov(est)

#-------------------------------------------------------------------------------
# Schweppe regression GM-estimator: Huber psi
xwgt <- simpsonWgt(employment / weighted_median(employment, weight), a = 1,
    b = 20)
est <- svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Schweppe",
    xwgt = xwgt)
expect_equal(coef(est),
    c("(Intercept)" = 7769.1899, employment = 14424.6617),
    tolerance = 1e-4, label = "Huber Schweppe regression GM-est")
# covariance matrix
ref <- matrix(c(35697.6180, -620.1554, -620.1554, 57.4449), ncol = 2)
colnames(ref) <- c("(Intercept)", "employment")
rownames(ref) <- c("(Intercept)", "employment")
expect_equal(vcov(est), ref, tolerance = 1e-4,
    label = "Huber Schweppe regression GM-est: cov")

# planned to fail...
expect_equal(1, 2, label = "not a real test")
