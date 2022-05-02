library(testthat)
suppressPackageStartupMessages(library(survey))
library(robsurvey, quietly = TRUE)

#===============================================================================
# MU284strat data
#===============================================================================
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
         weights = ~weights, data = MU284strat)
f <- RMT85 ~ P85 + S82 + CS82
# total
aux_totals <- c("(Intercept)" = 284, P85 = 8339, S82 = 13500, CS82 = 2583)
ref <- svytotal(~RMT85, calibrate(dn, f, aux_totals))
reg <- svyreg(f, dn)
est <- svytotal_reg(reg, aux_totals, type = "ADU")
expect_equal(coef(ref), coef(est),
    label = "MU284strat: GREG: total: coef")
expect_equal(as.numeric(SE(ref)), SE(est),
    label = "MU284strat: GREG: total: SE")
# mean
aux_means <- aux_totals / 284
ref <- svymean(~RMT85, calibrate(dn, f, aux_means))
reg <- svyreg(f, dn)
est <- svymean_reg(reg, aux_means, type = "ADU")
expect_equal(coef(ref), coef(est),
    label = "MU284strat: GREG: mean: coef")
expect_equal(as.numeric(SE(ref)), SE(est),
    label = "MU284strat: GREG: mean: SE")

#===============================================================================
# workplace data
#===============================================================================
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
f <- payroll ~ employment
# total
aux_totals <- c("(Intercept)" = 86514, employment = 953555)
ref <- svytotal(~payroll, calibrate(dn, f, aux_totals))
reg <- svyreg(f, dn)
est <- svytotal_reg(reg, aux_totals, type = "ADU")
expect_equal(coef(ref), coef(est),
    label = "workplace: GREG: total: coef")
expect_equal(as.numeric(SE(ref)), SE(est),
    label = "workplace: GREG: total: SE")
# mean
aux_means <- aux_totals / 86514
ref <- svytotal(~payroll, calibrate(dn, f, aux_means))
reg <- svyreg(f, dn)
est <- svymean_reg(reg, aux_means, type = "ADU")
expect_equal(coef(ref), coef(est),
    label = "workplace: GREG: total: coef")
expect_equal(as.numeric(SE(ref)), SE(est),
    label = "workplace: GREG: total: SE")
