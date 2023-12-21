suppressPackageStartupMessages(library("survey"))
library("robsurvey", quietly = TRUE)
source("check_functions.R")

#===============================================================================
# workplace data
#===============================================================================
dataset <- "workplace"
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
TOT_employment <- 1001233

#--------------------------------------
# ratio estimator (Huber)
ref <- svyratio(~payroll, ~employment, dn)
est <- svyratio_huber(~payroll, ~employment, dn, k = Inf)
check(ref, est, dataset, "svyratio_huber")

# ratio estimator of the total (with "hajek" variance estimator)
ref_tot <- predict(ref, TOT_employment)
est_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: SE"))

# ratio estimator using svyreg_huberM (with "hajek" variance estimator)
est_reg <- svyreg_huberM(payroll ~ -1 + employment, dn, var = ~employment,
    k = Inf)
check(ref, est_reg, dataset, "svyreg_huberM")

# ratio estimator of the total using svyreg_huberM
est_reg_tot <- svytotal_reg(est_reg, TOT_employment, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
check(est_reg_tot, est_rat_tot, dataset, "svyreg_huberM: svytotal_reg")

#--------------------------------------
# ratio estimator (Tukey)
ref <- svyratio(~payroll, ~employment, dn)
est <- svyratio_tukey(~payroll, ~employment, dn, k = Inf)
check(est, ref, dataset, "svyratio_tukey")

# ratio estimator of the total (with "hajek* variance estimator)
ref_tot <- predict(ref, TOT_employment)
est_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: SE"))

# ratio estimator using svyreg_tukeyM (with "hajek" variance estimator)
est_reg <- svyreg_tukeyM(payroll ~ -1 + employment, dn, var = ~employment,
    k = Inf)
check(est, est_reg, dataset, "svyreg_tukeyM")

# ratio estimator of the total using svyreg_tukeyM
est_reg_tot <- svytotal_reg(est_reg, TOT_employment, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
check(est_reg_tot, est_rat_tot, dataset, "svyreg_tukeyM: svytotal_reg")

#===============================================================================
# MU284strat data
#===============================================================================
dataset <- "MU284strat"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
         weights = ~weights, data = MU284strat)
TOT_P85 <- 8339

#--------------------------------------
# ratio estimator (Huber)
ref <- svyratio(~RMT85, ~P85, dn)
est <- svyratio_huber(~RMT85, ~P85, dn, k = Inf)
check(est, ref, dataset, "svyratio_huber")

# ratio estimator of the total (with "hajek" variance estimator)
ref_tot <- predict(ref, TOT_P85)
est_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: SE"))

# ratio estimator using svyreg_huberM (with "hajek" variance estimator)
est_reg <- svyreg_huberM(RMT85~ -1 + P85, dn, var = ~P85, k = Inf)
check(est, est_reg, dataset, "svyreg_huberM")

# ratio estimator of the total using svyreg_huberM
est_reg_tot <- svytotal_reg(est_reg, TOT_P85, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
check(est_reg_tot, est_rat_tot, dataset, "svyreg_huberM: svytotal_reg")

#--------------------------------------
# ratio estimator (Tukey)
ref <- svyratio(~RMT85, ~P85, dn)
est <- svyratio_tukey(~RMT85, ~P85, dn, k = Inf)
check(est, ref, dataset, "svyratio_tukey")

# ratio estimator of the total (with "hajek" variance estimator)
ref_tot <- predict(ref, TOT_P85)
est_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: SE"))

# ratio estimator using svyreg_tukeyM (with "hajek" variance estimator)
est_reg <- svyreg_tukeyM(RMT85~ -1 + P85, dn, var = ~P85, k = Inf)
check(est, est_reg, dataset, "svyreg_tukeyM")

# ratio estimator of the total using svyreg_tukeyM
est_reg_tot <- svytotal_reg(est_reg, TOT_P85, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
check(est_reg_tot, est_rat_tot, dataset, "svyreg_tukeyM: svytotal_reg")
