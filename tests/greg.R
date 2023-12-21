# NOTE: The estimators of the population mean for a model without a
#       regression intercept are excluded from the test because they
#       differ from the implementation based on survey::calibrate
suppressPackageStartupMessages(library("survey"))
library("robsurvey", quietly = TRUE)
source("check_functions.R")

#===============================================================================
# 1 MU284strat data (weights are calibrated to N)
#===============================================================================
dataset <- "MU284strat"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
         weights = ~weights, data = MU284strat)

f0 <- RMT85 ~ -1 + P85 + S82 + CS82
f <- RMT85 ~ P85 + S82 + CS82

N <- 284
aux_totals <- c(P85 = 8339, S82 = 13500, CS82 = 2583)

#--------------------------------------
# 1.1 total
#--------------------------------------

# GREG total with intercept
ref <- svytotal(~RMT85, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, dataset, "total with intercept")
# type "projective" and "ADU" coincide
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "projective")
check(ref, est, dataset, "total with intercept: projective")

# GREG total with intercept + heteroscedasticity
ref <- svytotal(~RMT85, calibrate(dn, f, c("(Intercept)" = N, aux_totals),
    variance = MU284strat$P85))
est <- svytotal_reg(svyreg(f, dn, var = ~ P85), aux_totals, N, type = "ADU")
check(ref, est, dataset, "total with intercept: hetero")

# GREG total without intercept
ref <- svytotal(~RMT85, calibrate(dn, f0, aux_totals))
est <- svytotal_reg(svyreg(f0, dn), aux_totals, type = "ADU")
check(ref, est, dataset, "total w/o intercept")

#--------------------------------------
# 1.2 mean
#--------------------------------------
# GREG mean with intercept
ref <- svymean(~RMT85, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svymean_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, dataset, "mean with intercept")

# GREG mean with intercept + heteroscedasticity
ref <- svymean(~RMT85, calibrate(dn, f, c("(Intercept)" = N, aux_totals),
    variance = MU284strat$P85))
est <- svymean_reg(svyreg(f, dn, var = ~ P85), aux_totals, N, type = "ADU")
check(ref, est, dataset, "mean with intercept: hetero")

# GREG mean without intercept, see NOTE:
if (FALSE) {
    ref <- svymean(~RMT85, calibrate(dn, f0, aux_totals))
    est <- svymean_reg(svyreg(f0, dn), aux_totals, N, type = "ADU")
    check(ref, est, dataset, "mean w/o intercept")
}
#===============================================================================
# 2 workplace data (weights are calibrated to N)
#===============================================================================
dataset <- "workplace"
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

f0 <- payroll ~ -1 + employment
f <- payroll ~ employment

N <- 86514
aux_totals <- c(employment = 953555 * 1.05)

#--------------------------------------
# 2.1 total
#--------------------------------------
# GREG total with intercept
ref <- svytotal(~payroll, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, dataset, "total with intercept")
# type "projective" and "ADU" coincide
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "projective")
check(ref, est, dataset, "total with intercept: projective")

# GREG total without intercept
ref <- svytotal(~payroll, calibrate(dn, f0, aux_totals))
est <- svytotal_reg(svyreg(f0, dn), aux_totals, type = "ADU")
check(ref, est, dataset, "total w/o intercept")

# GREG total without intercept + heteroscedasticity (ratio estimator)
ref <- svytotal(~payroll, calibrate(dn, f0, aux_totals,
    variance = workplace$employment))
est <- svytotal_reg(svyreg(f0, dn, var = ~employment), aux_totals,
    type = "ADU")
check(ref, est, dataset, "total w/o intercept: hetero")

#--------------------------------------
# 2.2 mean
#--------------------------------------
# GREG mean with intercept
ref <- svymean(~payroll, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svymean_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, dataset, "mean with intercept")

# GREG mean without intercept; see NOTE:
if (FALSE) {
    ref <- svymean(~payroll, calibrate(dn, f0, aux_totals))
    est <- svymean_reg(svyreg(f0, dn), aux_totals, N, type = "ADU")
    check(ref, est, dataset, "mean w/o intercept")
}

# GREG mean without intercept + heteroscedasticity (ratio estimator); see NOTE:
if (FALSE) {
    ref <- svymean(~payroll, calibrate(dn, f0, aux_totals,
        variance = workplace$employment))
    est <- svymean_reg(svyreg(f0, dn, var = ~employment), aux_totals, N,
        type = "ADU")
    check(ref, est, dataset, "mean w/o intercept: hetero")
}

#===============================================================================
# 3 workplace data (weights are NOT calibrated to N)
#===============================================================================
dataset <- "workplace"
# modify workplace data (mulitply weights of the first unit in each stratum
# by 0.25; thus, the weights are not calibrated)
wp <- workplace
wp[wp$ID == 2, "weight"] <- 0.25 * wp[wp$ID == 2, "weight"]
wp[wp$ID == 1, "weight"] <- 0.25 * wp[wp$ID == 1, "weight"]
wp[wp$ID == 7, "weight"] <- 0.25 * wp[wp$ID == 7, "weight"]
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = wp)

#--------------------------------------
# 3.1 total
#--------------------------------------
# GREG total with intercept
ref <- svytotal(~payroll, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, paste0(dataset, " (not calibrated)"), "total with intercept")
# type "projective" and "ADU" coincide
est <- svytotal_reg(svyreg(f, dn), aux_totals, N, type = "projective")
check(ref, est, paste0(dataset, " (not calibrated)"),
    "total with intercept: projective")

# GREG total without intercept
ref <- svytotal(~payroll, calibrate(dn, f0, aux_totals))
est <- svytotal_reg(svyreg(f0, dn), aux_totals, type = "ADU")
check(ref, est, paste0(dataset, " (not calibrated)"), "total w/o intercept")

# GREG total without intercept + heteroscedasticity (ratio estimator)
ref <- svytotal(~payroll, calibrate(dn, f0, aux_totals,
    variance = wp$employment))
est <- svytotal_reg(svyreg(f0, dn, var = ~employment), aux_totals,
    type = "ADU")
check(ref, est, paste0(dataset, " (not calibrated)"),
    "total w/o intercept: hetero")

#--------------------------------------
# 3.2 mean
#--------------------------------------
# GREG mean with intercept
ref <- svymean(~payroll, calibrate(dn, f, c("(Intercept)" = N, aux_totals)))
est <- svymean_reg(svyreg(f, dn), aux_totals, N, type = "ADU")
check(ref, est, paste0(dataset, " (not calibrated)"), "mean with intercept")

# GREG mean without intercept; see NOTE:
if (FALSE) {
    ref <- svymean(~payroll, calibrate(dn, f0, aux_totals))
    est <- svymean_reg(svyreg(f0, dn), aux_totals, N, type = "ADU")
    check(ref, est, paste0(dataset, " (not calibrated)"), "mean w/o intercept")
}

# GREG mean without intercept + heteroscedasticity (ratio estimator); see NOTE:
if (FALSE) {
    ref <- svymean(~payroll, calibrate(dn, f0, aux_totals,
        variance = wp$employment))
    est <- svymean_reg(svyreg(f0, dn, var = ~employment), aux_totals, N,
        type = "ADU")
    check(ref, est, paste0(dataset, " (not calibrated)"),
        "mean w/o intercept: hetero")
}
