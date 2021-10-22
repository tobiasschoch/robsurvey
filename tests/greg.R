library(testthat)
suppressPackageStartupMessages(library(survey))
library(robsurvey, quietly = TRUE)
data(workplace)
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

#-------------------------------------------------------------------------------
# GREG
est_beta <- svyreg(payroll ~ employment, dn)
est_mean <- svymean_reg(est_beta, c(1, 11.022))
# estimate of mean
expect_equal(unname(coef(est_mean)), 179591.9542025, label = "GREG")
# variance of variance
expect_equal(unname(SE(est_mean)), 11038.1134323, label = "GREG")

#FIXME:
pop.totals <- c(`(Intercept)` = sum(weights(dn)), coef(svytotal(~employment, dn)))
dn_cal <- calibrate(dn, payroll ~ employment, pop.totals)
svymean(~payroll, dn_cal)

# planned to fail...
# expect_equal(1, 2, label = "not a real test")
