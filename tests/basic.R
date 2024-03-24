suppressPackageStartupMessages(library(survey))
library("robsurvey", quietly = TRUE)
source("check_functions.R")

#===============================================================================
# 1 MU284 data
#===============================================================================
data("MU284strat"); data_name <- "MU284"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
    weights = ~weights, data = MU284strat, calibrate.formula = ~-1 + Stratum)
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
    data = workplace, calibrate.formula = ~-1 + strat)

# Reference estimates (against which we check)
sm <- svymean(~payroll, dn)
st <- svytotal(~payroll, dn)

#-------------------------------------------------------------------------------
# Huber M-estimator
check(sm, svymean_huber(~payroll, dn, k = Inf), data_name, "svymean_huber")
check(st, svytotal_huber(~payroll, dn, k = Inf), data_name, "svytotal_huber")

ref <- structure(list(estimate = 187670.8012254323985,
                      variance = 21858.046815482484817^2),
                 class = "svystat_rob")
check(ref, svymean_huber(~payroll, dn, k = 5, type = "rht"),
      data_name, "svymean_huber(k: 5, type: rht")


ref <- structure(list(estimate = 131136.23506721309968,
                      variance = 11475.091817854499823^2),
                 class = "svystat_rob")
check(ref, svymean_huber(~payroll, dn, k = 5, type = "rwm"),
      data_name, "svymean_huber(k: 5, type: rwm")

#-------------------------------------------------------------------------------
# Tukey M-estimator
check(sm, svymean_tukey(~payroll, dn, k = Inf), data_name, "svymean_tukey")
check(st, svytotal_tukey(~payroll, dn, k = Inf), data_name, "svytotal_tukey")

ref <- structure(list(estimate = 265306.64748404570855,
                      variance = 17947.151740625951788^2),
                 class = "svystat_rob")
check(ref, svymean_tukey(~payroll, dn, k = 5, type = "rht"),
      data_name, "svymean_tukey(k: 5, type: rht")

ref <- structure(list(estimate = 80131.795461709378287,
                      variance = 7146.9351801305610934^2),
                 class = "svystat_rob")
check(ref, svymean_tukey(~payroll, dn, k = 5, type = "rwm"),
      data_name, "svymean_tukey(k: 5, type: rwm")

#-------------------------------------------------------------------------------
# Trimming
check(sm, svymean_trimmed(~payroll, dn, LB = 0, UB = 1), data_name,
    "svymean_trimmed")
check(st, svytotal_trimmed(~payroll, dn, LB = 0, UB = 1), data_name,
    "svytotal_trimmed")

ref <- structure(list(estimate = 79140.83699598700332,
                      variance = 11370.017507516560727^2),
                 class = "svystat_rob")
check(ref, svymean_trimmed(~payroll, dn, LB = 0.1, UB = 0.8),
      data_name, "svymean_trimmed(LB: 0.1, UB: 0.8")

#-------------------------------------------------------------------------------
# Winsorized
check(sm, svymean_winsorized(~payroll, dn, LB = 0, UB = 1), data_name,
    "svymean_winsorized")
check(st, svytotal_winsorized(~payroll, dn, LB = 0, UB = 1), data_name,
    "svytotal_winsorized")

ref <- structure(list(estimate = 178730.03213352750754,
                      variance = 13810.51902068527852^2),
                 class = "svystat_rob")
check(ref, svymean_k_winsorized(~payroll, dn, k = 3),
      data_name, "svymean_k_winsorized(k: 3")

ref <- structure(list(estimate = 93888.43424185681215,
                      variance = 8834.4427633549585153^2),
                 class = "svystat_rob")
check(ref, svymean_winsorized(~payroll, dn, LB = 0.1, UB = 0.8),
      data_name, "svymean_winsorized(LB: 0.1, UB: 0.8")

#-------------------------------------------------------------------------------
# Dalen
check(sm, svymean_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svymean_dalen")
check(st, svytotal_dalen(~payroll, dn, censoring = 1e10, verbose = FALSE),
    data_name, "svytotal_dalen")

ref <- structure(list(estimate = 159624.55787502601743,
                      variance = 9054.0606693132558576^2),
                 class = "svystat_rob")
check(ref, svymean_dalen(~payroll, dn, censoring = 3e8, verbose = FALSE),
      data_name, "svymean_dalen(censoring: 3e8")
