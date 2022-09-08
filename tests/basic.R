library("testthat")
suppressPackageStartupMessages(library(survey))
library("robsurvey", quietly = TRUE)

#===============================================================================
# 1 MU284 data
#===============================================================================
data("MU284strat")
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
    weights = ~weights, data = MU284strat)
sm <- svymean(~RMT85, dn)
st <- svytotal(~RMT85, dn)
#-------------------------------------------------------------------------------
# Huber M-estimator
est_m <- svymean_huber(~RMT85, dn, k = Inf)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "MU284: mean: Huber: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "MU284: mean: Huber: coef")
est_t <- svytotal_huber(~RMT85, dn, k = Inf)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "MU284: total: Huber: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "MU284: total: Huber: coef")
#-------------------------------------------------------------------------------
# Tukey M-estimator
est_m <- svymean_tukey(~RMT85, dn, k = Inf)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "MU284: mean: Tukey: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "MU284: mean: Tukey: coef")
est_t <- svytotal_tukey(~RMT85, dn, k = Inf)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "MU284: total: Tukey: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "MU284: total: Tukey: coef")
#-------------------------------------------------------------------------------
# Trimming
est_m <- svymean_trimmed(~RMT85, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "MU284: mean: Trimmed: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "MU284: mean: Trimmed: coef")
est_t <- svytotal_trimmed(~RMT85, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "MU284: total: Trimmed: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "MU284: total: Trimmed: coef")
#-------------------------------------------------------------------------------
# Winsorized
est_m <- svymean_winsorized(~RMT85, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "MU284: mean: Winsorized: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "MU284: mean: Winsorized: coef")
est_t <- svytotal_winsorized(~RMT85, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "MU284: total: Winsorized: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "MU284: total: Winsorized: coef")
#-------------------------------------------------------------------------------
# Dalen
est_m <- svymean_dalen(~RMT85, dn, censoring = 1e10, verbose = FALSE)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "MU284: mean: Dalen: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "MU284: mean: Dalen: coef")
est_t <- svytotal_dalen(~RMT85, dn, censoring = 1e10, verbose = FALSE)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "MU284: total: Dalen: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "MU284: total: Dalen: coef")

#===============================================================================
# 2 workplace data
#===============================================================================
data("workplace")
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
sm <- svymean(~payroll, dn)
st <- svytotal(~payroll, dn)
#-------------------------------------------------------------------------------
# Huber M-estimator
est_m <- svymean_huber(~payroll, dn, k = Inf)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "workplace: mean: Huber: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "workplace: mean: Huber: coef")
est_t <- svytotal_huber(~payroll, dn, k = Inf)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "workplace: total: Huber: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "workplace: total: Huber: coef")
#-------------------------------------------------------------------------------
# Tukey M-estimator
est_m <- svymean_tukey(~payroll, dn, k = Inf)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "workplace: mean: Tukey: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "workplace: mean: Tukey: coef")
est_t <- svytotal_tukey(~payroll, dn, k = Inf)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "workplace: total: Tukey: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "workplace: total: Tukey: coef")
#-------------------------------------------------------------------------------
# Trimming
est_m <- svymean_trimmed(~payroll, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "workplace: mean: Trimmed: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "workplace: mean: Trimmed: coef")
est_t <- svytotal_trimmed(~payroll, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "workplace: total: Trimmed: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "workplace: total: Trimmed: coef")
#-------------------------------------------------------------------------------
# Winsorized
est_m <- svymean_winsorized(~payroll, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "workplace: mean: Winsorized: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "workplace: mean: Winsorized: coef")
est_t <- svytotal_winsorized(~payroll, dn, LB = 0, UB = 1)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "workplace: total: Winsorized: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "workplace: total: Winsorized: coef")
#-------------------------------------------------------------------------------
# Dalen
est_m <- svymean_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE)
expect_equal(as.numeric(SE(est_m)), as.numeric(SE(sm)),
    label = "workplace: mean: Dalen: SE")
expect_equal(as.numeric(coef(est_m)), as.numeric(coef(sm)),
    label = "workplace: mean: Dalen: coef")
est_t <- svytotal_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE)
expect_equal(as.numeric(SE(est_t)), as.numeric(SE(st)),
    label = "workplace: total: Dalen: SE")
expect_equal(as.numeric(coef(est_t)), as.numeric(coef(st)),
    label = "workplace: total: Dalen: coef")
