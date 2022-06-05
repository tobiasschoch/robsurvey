library(testthat)
suppressPackageStartupMessages(library(survey))
library(robsurvey, quietly = TRUE)

check <- function(ref, est, dataset, characteristic)
{
    expect_equal(as.numeric(coef(ref)), as.numeric(coef(est)),
        label = paste0(dataset, ": GREG: ", characteristic, ": coef"))
    expect_equal(as.numeric(SE(ref)), as.numeric(SE(est)),
        label = paste0(dataset, ": GREG: ", characteristic, ": SE"))
}

greg <- function(formula, design, aux, type = c("total", "mean"),
    gweights = FALSE)
{
    type <- match.arg(type, c("total", "mean"))
    mf <- stats::model.frame(formula, design$variables)
    w <- weights(design)
    x <- stats::model.matrix(stats::terms(mf), mf)
    m <- lm.wfit(x, as.numeric(stats::model.response(mf)), w)

    res <- if (gweights) {
        t_hat <- colSums(w * x)
        if (type == "mean")
            t_hat <- t_hat / sum(w)
        g <- 1 + as.numeric(t(aux - t_hat) %*% tcrossprod(backsolve(qr.R(m$qr),
            diag(m$rank))) %*% t(x))
        g * residuals(m)
    } else {
        residuals(m)
    }

    dn_update <- update(dn, res = res)

    est <- sum(aux * coef(m))
    switch(type,
        "total" = {
            est <- est + sum(w * residuals(m))
            se <- SE(svytotal(~res, dn_update))
            attr(est, "statistic") <- "total"

        },
        "mean" = {
            est <- est + sum(w * residuals(m)) / sum(w)
            se <- SE(svymean(~res, dn_update))
            attr(est, "statistic") <- "mean"
        })
    attr(est, "var") <- se^2
    class(est) <- "svystat"
    est
}

#===============================================================================
# MU284strat data
#===============================================================================
dataset <- "MU284strat"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
         weights = ~weights, data = MU284strat)

# GREG total with intercept
f <- RMT85 ~ P85 + S82 + CS82
aux_totals <- c("(Intercept)" = 284, P85 = 8339, S82 = 13500, CS82 = 2583)
#ref <- svytotal(~RMT85, calibrate(dn, f, aux_totals))
ref <- greg(f, dn, aux_totals, "total")
est <- svytotal_reg(svyreg(f, dn), aux_totals, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "total with intercept")

# GREG total without intercept
f <- RMT85 ~ -1 + P85 + S82 + CS82
aux_totals <- c(P85 = 8339, S82 = 13500, CS82 = 2583)
#ref <- svytotal(~RMT85, calibrate(dn, f, aux_totals))
ref <- greg(f, dn, aux_totals, "total")
est <- svytotal_reg(svyreg(f, dn), aux_totals, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "total w/o intercept")

# GREG mean with intercept
f <- RMT85 ~ P85 + S82 + CS82
aux_means <- c("(Intercept)" = 284, P85 = 8339, S82 = 13500, CS82 = 2583) / 284
#ref <- svymean(~RMT85, calibrate(dn, f, aux_means))
ref <- greg(f, dn, aux_means, "mean")
est <- svymean_reg(svyreg(f, dn), aux_means, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "mean with intercept")

# GREG mean without intercept
f <- RMT85 ~ -1 + P85 + S82 + CS82
aux_means <- c(P85 = 8339, S82 = 13500, CS82 = 2583) / 284
#ref <- svymean(~RMT85, calibrate(dn, f, aux_means))
ref <- greg(f, dn, aux_means, "mean")
est <- svymean_reg(svyreg(f, dn), aux_means, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "mean w/o intercept")

#===============================================================================
# workplace data
#===============================================================================
dataset <- "workplace"
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# GREG total with intercept
f <- payroll ~ employment
aux_totals <- c("(Intercept)" = 86514, employment = 953555 * 1.05)
#ref <- svytotal(~payroll, calibrate(dn, f, aux_totals))
ref <- greg(f, dn, aux_totals, "total")
est <- svytotal_reg(svyreg(f, dn), aux_totals, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "total with intercept")

# GREG total without intercept
f <- payroll ~ -1 + employment
aux_totals <- c(employment = 953555 * 1.05)
#ref <- svytotal(~payroll, calibrate(dn, f, aux_totals))
ref <- greg(f, dn, aux_totals, "total")
est <- svytotal_reg(svyreg(f, dn), aux_totals, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "total w/o intercept")

# GREG mean with intercept
f <- payroll ~ employment
aux_means <- c("(Intercept)" = 86514, employment = 953555) / 86514
#ref <- svymean(~payroll, calibrate(dn, f, aux_means))
ref <- greg(f, dn, aux_means, "mean")
est <- svymean_reg(svyreg(f, dn), aux_means, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "mean with intercept")

# GREG mean without intercept
f <- payroll ~ -1 + employment
aux_means <- c(employment = 953555) / 86514
#ref <- svymean(~payroll, calibrate(dn, f, aux_means))
ref <- greg(f, dn, aux_means, "mean")
est <- svymean_reg(svyreg(f, dn), aux_means, type = "ADU", verbose = FALSE)
check(ref, est, dataset, "mean w/o intercept")

