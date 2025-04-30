# weighted winsorized mean (depends on pkg survey)
svymean_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB,
                               na.rm = FALSE, trim_var = FALSE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Weighted winsorized estimator (", LB, ", ", UB, ")"),
            dat$design, match.call(), "wins", LB = LB, UB = UB))

    # population- vs. domain-level estimate (matters for calibrated designs)
    res <- weighted_mean_winsorized(dat$y, dat$w, LB, UB, TRUE, FALSE)

    # influence function
    infl <- if (trim_var)
        .infl_trimmed(res$model$y, res$model$w, LB, UB, res$estimate)
    else
        .infl_winsorized(res$model$y, res$model$w, LB, UB, res$estimate)

    infl <- infl * res$model$w / sum(res$model$w)
    if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# weighted one-sided k winsorized mean (depends on pkg survey)
svymean_k_winsorized <- function(x, design, k, na.rm = FALSE,
                                 trim_var = FALSE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Weighted winsorized estimator (", "k = ", k),
            dat$design, match.call(), "wins", k = k))

    # population- vs. domain-level estimate (matters for calibrated designs)
    res <- weighted_mean_k_winsorized(dat$y, dat$w, k, TRUE, FALSE)

    # influence function
    infl <- if (trim_var)
        .infl_trimmed(res$model$y, res$model$w, 0, res$estimator$UB,
                      res$estimate)
    else
        .infl_winsorized(res$model$y, res$model$w, 0, res$estimator$UB,
                         res$estimate)

    infl <- infl * res$model$w / sum(res$model$w)
    if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# weighted winsorized total (depends on pkg survey)
svytotal_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB,
                                na.rm = FALSE, trim_var = FALSE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Weighted winsorized estimator (", LB, ", ", UB, ")"),
            dat$design, match.call(), "wins", LB = LB, UB = UB))

    # population- vs. domain-level estimate (matters for calibrated designs)
    res <- weighted_total_winsorized(dat$y, dat$w, LB, UB, TRUE, FALSE)

    # influence function
    infl <- if (trim_var)
        .infl_trimmed(res$model$y, res$model$w, LB, UB, 0)
    else
        .infl_winsorized(res$model$y, res$model$w, LB, UB, 0)

    infl <- infl * res$model$w
    if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# weighted one-sided k winsorized total (depends on pkg survey)
svytotal_k_winsorized <- function(x, design, k, na.rm = FALSE,
                                  trim_var = FALSE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Weighted winsorized estimator (", "k = ", k),
            dat$design, match.call(), "wins", k = k))

    # population- vs. domain-level estimate(matters for calibrated designs )
    res <- weighted_total_k_winsorized(dat$y, dat$w, k, TRUE, FALSE)

    # influence function
    infl <- if (trim_var)
        .infl_trimmed(res$model$y, res$model$w, 0, res$estimator$UB, 0)
    else
        .infl_winsorized(res$model$y, res$model$w, 0, res$estimator$UB, 0)

    infl <- infl * res$model$w
    if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# influence function, Huber (1981, p. 58-59)
.infl_winsorized <- function(x, w, LB, UB, wm, ngrid = 401)
{
    x <- as.double(x); w <- as.double(w)
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    # density estimates at 'LB' and '1 - UB'
    bwd <- dpik(x, scalest = "minim", level = 2, kernel = "normal",
                canonical = FALSE, gridsize = ngrid, range.x = range(x),
                truncate = TRUE)
    at <- seq(min(x), max(x), length = ngrid)
    nx <- rowsum(c(rep(0, ngrid), w), c(1:ngrid, findInterval(x, at)))
    dens <- locpoly(rep(1, ngrid), nx * ngrid / (diff(range(x)) * sum(w)),
                    binned = TRUE, bandwidth = bwd, range.x = range(x))
    f_LB <- dens$y[min(which(dens$x >= qs[1]))]
    f_UB <- dens$y[max(which(dens$x <= qs[2]))]
    # influence function
    infl <- pmin.int(qs[2] + (1 - UB) / f_UB, pmax.int(qs[1] - LB / f_LB, x))
    # functional W corresponding to winsorized mean
    W <- unname((UB - LB) * wm + LB * qs[1] + (1 - UB) * qs[2])
    infl <- infl - W
    if (f_LB > 0)
        infl <- infl + LB^2 / f_LB
    if (f_UB > 0)
        infl <- infl + (1 - UB)^2 / f_UB
    infl
}
