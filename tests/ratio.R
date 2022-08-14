all_equal <- function(target, current, label,
    tolerance = sqrt(.Machine$double.eps), scale = NULL,
    check.attributes = FALSE)
{
    if (missing(label))
        stop("Argument 'label' is missing\n")
    res <- all.equal(target, current, tolerance, scale,
        check.attributes = check.attributes)
    if (is.character(res))
        cat(paste0(label, ": ", res, "\n"))
}
suppressPackageStartupMessages(library(survey))
library(robsurvey, quietly = TRUE)

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
all_equal(coef(ref), coef(est), paste0(dataset, ": svyratio_huber: coef"))
all_equal(SE(ref), SE(est), paste0(dataset, ": svyratio_huber: SE"))

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
all_equal(coef(est), coef(est_reg), paste0(dataset, ": svyreg_huberM: coef"))
all_equal(SE(est), SE(est_reg), paste0(dataset, ": svyreg_huberM: SE"))

# ratio estimator of the total using svyreg_huberM
est_reg_tot <- svytotal_reg(est_reg, TOT_employment, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
all_equal(coef(est_reg_tot), coef(est_rat_tot),
    paste0(dataset, ": svyreg_huberM: svytotal_reg: coef"))
all_equal(SE(est_reg_tot), SE(est_rat_tot),
    paste0(dataset, ": svyreg_huberM: svytotal_reg: SE"))

#--------------------------------------
# ratio estimator (Tukey)
ref <- svyratio(~payroll, ~employment, dn)
est <- svyratio_tukey(~payroll, ~employment, dn, k = Inf)
all_equal(coef(ref), coef(est), paste0(dataset, ": svyratio_tukey: coef"))
all_equal(SE(ref), SE(est), paste0(dataset, ": svyratio_tukey: SE"))

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
all_equal(coef(est), coef(est_reg), paste0(dataset, ": svyreg_tukeyM: coef"))
all_equal(SE(est), SE(est_reg), paste0(dataset, ": svyreg_tukeyM: SE"))

# ratio estimator of the total using svyreg_tukeyM
est_reg_tot <- svytotal_reg(est_reg, TOT_employment, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_employment, variance = "hajek")
all_equal(coef(est_reg_tot), coef(est_rat_tot),
    paste0(dataset, ": svyreg_tukeyM: svytotal_reg: coef"))
all_equal(SE(est_reg_tot), SE(est_rat_tot),
    paste0(dataset, ": svyreg_tukeyM: svytotal_reg: SE"))

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
all_equal(coef(ref), coef(est), paste0(dataset, ": svyratio_huber: coef"))
all_equal(SE(ref), SE(est), paste0(dataset, ": svyratio_huber: SE"))

# ratio estimator of the total (with "hajek" variance estimator)
ref_tot <- predict(ref, TOT_P85)
est_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_huber: svytotal_ratio: SE"))

# ratio estimator using svyreg_huberM (with "hajek" variance estimator)
est_reg <- svyreg_huberM(RMT85~ -1 + P85, dn, var = ~P85, k = Inf)
all_equal(coef(est), coef(est_reg), paste0(dataset, ": svyreg_huberM: coef"))
all_equal(SE(est), SE(est_reg), paste0(dataset, ": svyreg_huberM: SE"))

# ratio estimator of the total using svyreg_huberM
est_reg_tot <- svytotal_reg(est_reg, TOT_P85, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(coef(est_reg_tot), coef(est_rat_tot),
    paste0(dataset, ": svyreg_huberM: svytotal_reg: coef"))
all_equal(SE(est_reg_tot), SE(est_rat_tot),
    paste0(dataset, ": svyreg_huberM: svytotal_reg: SE"))

#--------------------------------------
# ratio estimator (Tukey)
ref <- svyratio(~RMT85, ~P85, dn)
est <- svyratio_tukey(~RMT85, ~P85, dn, k = Inf)
all_equal(coef(ref), coef(est), paste0(dataset, ": svyratio_tukey: coef"))
all_equal(SE(ref), SE(est), paste0(dataset, ": svyratio_tukey: SE"))

# ratio estimator of the total (with "hajek" variance estimator)
ref_tot <- predict(ref, TOT_P85)
est_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(as.numeric(ref_tot$total), coef(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: coef"))
all_equal(as.numeric(ref_tot$se), SE(est_tot),
    paste0(dataset, ": svyratio_tukey: svytotal_ratio: SE"))

# ratio estimator using svyreg_tukeyM (with "hajek" variance estimator)
est_reg <- svyreg_tukeyM(RMT85~ -1 + P85, dn, var = ~P85, k = Inf)
all_equal(coef(est), coef(est_reg), paste0(dataset, ": svyreg_tukeyM: coef"))
all_equal(SE(est), SE(est_reg), paste0(dataset, ": svyreg_tukeyM: SE"))

# ratio estimator of the total using svyreg_tukeyM
est_reg_tot <- svytotal_reg(est_reg, TOT_P85, type = "projective")
est_rat_tot <- svytotal_ratio(est, TOT_P85, variance = "hajek")
all_equal(coef(est_reg_tot), coef(est_rat_tot),
    paste0(dataset, ": svyreg_tukeyM: svytotal_reg: coef"))
all_equal(SE(est_reg_tot), SE(est_rat_tot),
    paste0(dataset, ": svyreg_tukeyM: svytotal_reg: SE"))
