# Dalen the weighted mean (depends on pkg survey)
svymean_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
                          verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname, paste0("Dalen ", type,
            " estimator (censored at ", censoring, ")"), dat$domain,
            dat$design, match.call(), "dalen", censoring = censoring))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_mean_dalen(dat$y[dat$in_domain], dat$w[dat$in_domain],
                            censoring, type, TRUE, FALSE, verbose)
    else
        weighted_mean_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
                            verbose)
    # influence function
    infl <- if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- (res$robust$xw - res$model$w * res$estimate) /
            sum(res$model$w)
        tmp
    } else {
        (res$robust$xw - res$model$w * res$estimate) / sum(res$model$w)
    }

    # compute variance
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
# Dalen the weighted total (depends on pkg survey)
svytotal_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
                           verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname, paste0("Dalen ", type,
            " estimator (censored at ", censoring, ")"), dat$domain,
            dat$design, match.call(), "dalen", censoring = censoring))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_total_dalen(dat$y[dat$in_domain], dat$w[dat$in_domain],
                             censoring, type, TRUE, FALSE, verbose)
    else
        weighted_total_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
                             verbose)
    # influence function
    infl <- if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- res$robust$xw
        tmp
    } else {
        res$robust$xw
    }

    # compute variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata,
                              design$fpc, postStrata = design$postStrata)
    # return
    names(res$estimate) <- dat$yname
    res$estimator$domain <- dat$domain
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- "svystat_rob"
    res
}
