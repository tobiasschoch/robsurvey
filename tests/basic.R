suppressPackageStartupMessages(library(survey))
library("robsurvey", quietly = TRUE)
source("check_functions.R")

#===============================================================================
# 1 MU284 data
#===============================================================================
data("MU284strat"); data_name <- "MU284"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
    weights = ~weights, data = MU284strat)

# Reference estimates (against which we check)
sm <- svymean(~RMT85, dn)
st <- svytotal(~RMT85, dn)

#-------------------------------------------------------------------------------
# Huber M-estimator
check(sm, svymean_huber(~RMT85, dn, k = Inf), data_name, "svymean_huber")
check(st, svytotal_huber(~RMT85, dn, k = Inf), data_name, "svytotal_huber")

#-------------------------------------------------------------------------------
# Tukey M-estimator
check(sm, svymean_tukey(~RMT85, dn, k = Inf), data_name, "svymean_tukey")
check(st, svytotal_tukey(~RMT85, dn, k = Inf), data_name, "svytotal_tukey")

#-------------------------------------------------------------------------------
# Trimming
check(sm, svymean_trimmed(~RMT85, dn, LB = 0, UB = 1), data_name,
    "svymean_trimmed")
check(st, svytotal_trimmed(~RMT85, dn, LB = 0, UB = 1), data_name,
    "svytotal_trimmed")

#-------------------------------------------------------------------------------
# Winsorized
check(sm, svymean_winsorized(~RMT85, dn, LB = 0, UB = 1), data_name,
    "svymean_winsorized")
check(st, svytotal_winsorized(~RMT85, dn, LB = 0, UB = 1), data_name,
    "svytotal_winsorized")

#-------------------------------------------------------------------------------
# Dalen
check(sm, svymean_dalen(~RMT85, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svymean_dalen")
check(st, svytotal_dalen(~RMT85, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svytotal_dalen")

#===============================================================================
# 2 workplace data
#===============================================================================
data("workplace"); data_name <- "workplace"
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Reference estimates (against which we check)
sm <- svymean(~payroll, dn)
st <- svytotal(~payroll, dn)

#-------------------------------------------------------------------------------
# Huber M-estimator
check(sm, svymean_huber(~payroll, dn, k = Inf), data_name, "svymean_huber")
check(st, svytotal_huber(~payroll, dn, k = Inf), data_name, "svytotal_huber")

#-------------------------------------------------------------------------------
# Tukey M-estimator
check(sm, svymean_tukey(~payroll, dn, k = Inf), data_name, "svymean_tukey")
check(st, svytotal_tukey(~payroll, dn, k = Inf), data_name, "svytotal_tukey")

#-------------------------------------------------------------------------------
# Trimming
check(sm, svymean_trimmed(~payroll, dn, LB = 0, UB = 1), data_name,
    "svymean_trimmed")
check(st, svytotal_trimmed(~payroll, dn, LB = 0, UB = 1), data_name,
    "svytotal_trimmed")

#-------------------------------------------------------------------------------
# Winsorized
check(sm, svymean_winsorized(~payroll, dn, LB = 0, UB = 1), data_name,
    "svymean_winsorized")
check(st, svytotal_winsorized(~payroll, dn, LB = 0, UB = 1), data_name,
    "svytotal_winsorized")

#-------------------------------------------------------------------------------
# Dalen
check(sm, svymean_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svymean_dalen")
check(st, svytotal_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svytotal_dalen")

