# weighted trimmed mean (depends on pkg survey)
svymean_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE,
                            ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Weighted trimmed estimator (", LB, ", ", UB, ")"),
            dat$domain, dat$design, match.call(), "trim", LB = LB, UB = UB))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_mean_trimmed(dat$y[dat$in_domain], dat$w[dat$in_domain], LB,
                              UB, TRUE, FALSE)
    else
        weighted_mean_trimmed(dat$y, dat$w, LB, UB, TRUE, FALSE)

    # influence function
    infl <- .infl_trimmed(res$model$y, res$model$w, LB, UB, res$estimate) *
                res$model$w / sum(res$model$w)
    if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$estimator$domain <- dat$domain
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# weighted trimmed total (depends on pkg survey)
svytotal_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE,
                             ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Weighted trimmed estimator (", LB, ", ", UB, ")"),
            dat$domain, dat$design, match.call(), "trim", LB = LB, UB = UB))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_total_trimmed(dat$y[dat$in_domain], dat$w[dat$in_domain], LB,
                               UB, TRUE, FALSE)
    else
        weighted_total_trimmed(dat$y, dat$w, LB, UB, TRUE, FALSE)

    # influence function
    infl <- .infl_trimmed(res$model$y, res$model$w, LB, UB, 0) * res$model$w
    if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }

    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$estimator$domain <- dat$domain
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
# influence function, Huber (1981, p. 58)
.infl_trimmed <- function(x, w, LB, UB, tm)
{
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    x_wins <- pmin.int(qs[2], pmax.int(qs[1], x))
    # functional W corresponding to winsorized mean
    W <- (UB - LB) * tm + LB * qs[1] + (1 - UB) * qs[2]
    # return
    (x_wins - W) / (UB - LB)
}
