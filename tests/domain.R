suppressPackageStartupMessages(library(survey))
library("robsurvey", quietly = TRUE)
source("check_functions.R")

# population and sample size (by strata)
N <- c(8, 20, 100, 200)
n <- c(6, 10, 12, 190)

# variables
strat <- unlist(Map(rep, seq_along(n), n))
fpc <- unlist(Map(rep, N, n))

set.seed(1024)
x <- unlist(Map(rnorm, n, seq_along(n) * 10))

# domain indicators
select <- function(n, prop)
{
    sel <- rep(FALSE, n)
    sel[seq_len(floor(prop * n))] <- TRUE
    sel
}

d60 <- unlist(Map(select, n, 0.6))
d80 <- unlist(Map(select, n, 0.8))
d90 <- unlist(Map(select, n, 0.9))

# design weights
w <- unlist(Map(rep, N / n, n))

# data set and design
dat <- data.frame(strat, x, d60, d80, d90, w, fpc)
dn <- svydesign(~1, strat = ~strat, fpc = ~fpc, data = dat)
dn_cal <- svydesign(~1, strat = ~strat, fpc = ~fpc, data = dat,
                    calibrate.formula = ~strat)

# set missing values in variable x (of a given design)
set_missing <- function(design, variable)
{
    design$variables$x[design$variables[, variable] == FALSE] <- NA
    design
}

#-------------------------------------------------------------------------------
# domain
dn_d60 <- subset(dn, d60 == TRUE)                       # original design
check(svymean(~x, dn_d60),
      svymean_trimmed(~x, dn_d60, LB = 0, UB = 1),
      "domain d60", "mean")

check(svyby(~x, ~strat, dn_d60, svymean),
      svyby(~x, ~strat, dn_d60, svymean_trimmed, LB = 0, UB = 1),
      "domain d60", "mean")

dn_cal_d60 <- subset(dn_cal, d60 == TRUE)               # calibrated
check(svymean(~x, dn_cal_d60),
      svymean_trimmed(~x, dn_cal_d60, LB = 0, UB = 1),
      "domain d60", "mean, calibrated")

check(svyby(~x, ~strat, dn_cal_d60, svymean),
      svyby(~x, ~strat, dn_cal_d60, svymean_trimmed, LB = 0, UB = 1),
      "domain d60", "mean, calibrated")

#-------------------------------------------------------------------------------
# missing values
dn_miss_d60 <- set_missing(dn, "d60")                   # original design
check(svymean(~x, dn_miss_d60, na.rm = TRUE),
      svymean_trimmed(~x, dn_miss_d60, LB = 0, UB = 1, na.rm = TRUE),
      "missing values: d60", "mean")

check(svyby(~x, ~strat, dn_miss_d60, svymean, na.rm = TRUE),
      svyby(~x, ~strat, dn_miss_d60, svymean_trimmed, LB = 0, UB = 1,
            na.rm = TRUE),
      "missing values: d60", "mean")

dn_cal_miss_d60 <- set_missing(dn_cal, "d60")           # calibrated
check(svymean(~x, dn_cal_miss_d60, na.rm = TRUE),
      svymean_trimmed(~x, dn_cal_miss_d60, LB = 0, UB = 1, na.rm = TRUE),
      "missing values: d60", "mean, calibrated")

check(svyby(~x, ~strat, dn_cal_miss_d60, svymean, na.rm = TRUE),
      svyby(~x, ~strat, dn_cal_miss_d60, svymean_trimmed, LB = 0, UB = 1,
            na.rm = TRUE),
      "missing values: d60", "mean, calibrated")

#-------------------------------------------------------------------------------
# missing values and domain estimation
# original design
dn_d90_miss_d60 <- subset(dn_miss_d60, d90 == TRUE)
check(svymean(~x, dn_d90_miss_d60, na.rm = TRUE),
      svymean_trimmed(~x, dn_d90_miss_d60, LB = 0, UB = 1, na.rm = TRUE),
      "domain d90, missing values: d60", "mean")

check(svyby(~x, ~strat, dn_d90_miss_d60, svymean, na.rm = TRUE),
      svyby(~x, ~strat, dn_d90_miss_d60, svymean_trimmed, LB = 0, UB = 1,
            na.rm = TRUE),
      "domain d90, missing values: d60", "mean")

# calibrated
dn_cal_d90_miss_d60 <- subset(dn_cal_miss_d60, d90 == TRUE)
check(svymean(~x, dn_cal_d90_miss_d60, na.rm = TRUE),
      svymean_trimmed(~x, dn_cal_d90_miss_d60, LB = 0, UB = 1, na.rm = TRUE),
      "domain d90, missing values: d60", "mean, calibrated")

check(svyby(~x, ~strat, dn_cal_d90_miss_d60, svymean, na.rm = TRUE),
      svyby(~x, ~strat, dn_cal_d90_miss_d60, svymean_trimmed, LB = 0, UB = 1,
            na.rm = TRUE),
      "domain d90, missing values: d60", "mean, calibrated")

