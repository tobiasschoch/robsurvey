suppressPackageStartupMessages(library("survey"))
library("robsurvey", quietly = TRUE)
attach(workplace)

# design object
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

source("check_functions.R")

#===============================================================================
# 1 Test against MASS::rlm
#===============================================================================
if (requireNamespace("MASS", quietly = TRUE)) {
    library(MASS)
    # make a copy of MASS::rlm
    rlm_mod <- MASS:::rlm.default
    # replace wmad-function with our weighted mad
    body(rlm_mod)[[4]] <- substitute(wmad <- function(x, w) weighted_mad(x, w))

    #---------------------------------------------------------------------------
    # Huber regression M-estimator
    est <- svyreg_huberM(payroll ~ employment, dn, k = 1.345, tol = 1e-9)
    ref <- rlm_mod(x = est$model$x, y = est$model$y, weights = est$model$w,
        k = 1.345, scale.est = "MAD", method = "M", wt.method = "case",
        acc = 1e-9, maxit = 50, test.vec = "coef")
    all_equal(coef(est), coef(ref), tolerance = 1e-7,
        label = "Huber regression regression: coefficients")
    # model-based covariance matrix
    all_equal(vcov(est, "model"), vcov(ref), tolerance = 1e-7,
        label = "Huber regression M-est: model-based cov")
    # design-based covariance matrix (test against fixed values)
    ref <- matrix(c(61704317.084028, -3367311.1407500, -3367311.140750,
        870905.6027803), ncol = 2)
    colnames(ref) <- rownames(ref) <- c("(Intercept)", "employment")
    all_equal(vcov(est), ref,
        label = "Huber regression M-est: design-based cov")

    #---------------------------------------------------------------------------
    # Tukey regression M-estimator
    est <- svyreg_tukeyM(payroll ~ employment, dn, k = 4.6, tol = 1e-9)
    ref <- rlm_mod(x = est$model$x, y = est$model$y, weights = est$model$w,
        c = 4.6, psi = psi.bisquare,  scale.est = "MAD", method = "M",
        wt.method = "case", acc = 1e-9, maxit = 50, test.vec = "coef")
    all_equal(coef(est), coef(ref), tolerance = 1e-7,
        label = "Tukey regression regression: coefficients")
    # model-based covariance matrix
    all_equal(vcov(est, "model"), vcov(ref), tolerance = 1e-7,
        label = "Tukey regression M-est: model-based cov")
    # design-based covariance matrix (test against fixed values)
    ref <- matrix(c(43513457.824988, -1421272.808406, -1421272.808406,
        589121.583227), ncol = 2)
    colnames(ref) <- rownames(ref) <- c("(Intercept)", "employment")
    all_equal(vcov(est), ref,
        label = "Tukey regression M-est: design-based cov")
}

#===============================================================================
# 2 Test against survey::svyglm
#===============================================================================
# survey weighted regression
ref <- svyglm(payroll ~ employment, dn)
est <- svyreg(payroll ~ employment, dn)

all_equal(coef(est), coef(ref),
    label = "Survey weighted regression: coefficients")
# design-based covariance matrix (test against survey::svyglm)
all_equal(vcov(est, "design"), vcov(ref),
    label = "Survey weighted regression: design-based cov")
# model-based covariance matrix (test against MASS::rlm)
ref <- rlm(payroll ~ employment, weights = weights(dn), data = dn$variables,
    method = "M", wt.method = "case", k = 1000)
all_equal(vcov(est, "model"), vcov(ref),
    label = "Survey weighted regression: model-based cov")

#-------------------------------------------------------------------------------
# Mallows regression GM-estimator: Huber psi
xwgt <- simpsonWgt(employment / weighted_median(employment, weight), a = 1,
    b = 20)
est <- svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Mallows",
    xwgt = xwgt)
all_equal(coef(est),
    c("(Intercept)" = 7843.4909, employment = 14425.3658),
    tolerance = 1e-4, label = "Huber Mallows regression GM-est")
# model-based covariance matrix
ref <- matrix(c(35768.4163, -572.71408, -572.71408, 52.40226), ncol = 2)
colnames(ref) <- c("(Intercept)", "employment")
rownames(ref) <- c("(Intercept)", "employment")
all_equal(vcov(est, "model"), ref, tolerance = 1e-4,
    label = "Huber Mallows regression GM-est: model-based cov")

#-------------------------------------------------------------------------------
# Schweppe regression GM-estimator: Huber psi
xwgt <- simpsonWgt(employment / weighted_median(employment, weight), a = 1,
    b = 20)
est <- svyreg_huberGM(payroll ~ employment, dn, k = 1.345, type = "Schweppe",
    xwgt = xwgt)
all_equal(coef(est),
    c("(Intercept)" = 7769.1899, employment = 14424.6617),
    tolerance = 1e-4, label = "Huber Schweppe regression GM-est")
# model-based covariance matrix
ref <- matrix(c(35697.6180, -620.1554, -620.1554, 57.4449), ncol = 2)
colnames(ref) <- c("(Intercept)", "employment")
rownames(ref) <- c("(Intercept)", "employment")
all_equal(vcov(est, "model"), ref, tolerance = 1e-4,
    label = "Huber Schweppe regression GM-est: model-based cov")
