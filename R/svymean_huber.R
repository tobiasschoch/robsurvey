# Huber M-estimator of the weighted mean (depends on pkg survey)
svymean_huber <- function(x, design, k, type = "rwm", asym = FALSE,
                          na.rm = FALSE, verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)

    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Huber M-estimator (type = ", type,
            if (asym) "; asym. psi" else "", ")"), dat$domain,
            dat$design, match.call(), "mest", type = type,
            psi = if (asym) 1 else 0, psi_fun = "Huber", k = k))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_mean_huber(dat$y[dat$in_domain], dat$w[dat$in_domain], k,
                            type, asym, TRUE, FALSE, verbose, ...)
    else
        weighted_mean_huber(dat$y, dat$w, k, type, asym, TRUE, FALSE, verbose,
                            ...)

    # modify residuals for type 'rht' (only for variance estimation)
    r <- if (type == "rht")
        sqrt(res$model$var) * res$model$y - res$estimate
    else
        res$residuals

    # influence function
    infl <- if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- res$robust$robweights * r * res$model$w /
            sum(res$model$w)
        tmp
    } else {
        res$robust$robweights * r * res$model$w / sum(res$model$w)
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
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
# Huber M-estimator of the weighted total (depends on pkg survey)
svytotal_huber <- function(x, design, k, type = "rwm", asym = FALSE,
                           na.rm = FALSE, verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Huber M-estimator (type = ", type,
            if (asym) "; asym. psi" else "", ")"), dat$domain,
            dat$design, match.call(), "mest", type = type,
            psi = if (asym) 1 else 0, psi_fun = "Huber", k = k))

    # population- vs. domain-level estimate
    res <- if (dat$domain)
        weighted_total_huber(dat$y[dat$in_domain], dat$w[dat$in_domain], k,
                            type, asym, TRUE, FALSE, verbose, ...)
    else
        weighted_total_huber(dat$y, dat$w, k, type, asym, TRUE, FALSE, verbose,
                            ...)
    # influence function
    infl <- if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- res$robust$robweights * res$model$w *
            res$model$y
        tmp
    } else {
        res$robust$robweights * res$model$y * res$model$w
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
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
