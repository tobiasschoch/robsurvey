# Tukey biweight M-estimator of the weighted mean (depends on pkg survey)
svymean_tukey <- function(x, design, k, type = "rwm", na.rm = FALSE,
                          verbose = TRUE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Tukey M-estimator (type = ", type, ")"), dat$design,
            match.call(), "mest", type = type, psi = 2, psi_fun = "Tukey",
            k = k))

    # population- vs. domain-level estimate (matters for calibrated designs)
    res <- weighted_mean_tukey(dat$y, dat$w, k, type, TRUE, FALSE, verbose,
                               ...)

    # modify residuals for type 'rht' (only for variance estimation)
    r <- if (type == "rht")
        sqrt(res$model$var) * res$model$y - res$estimate
    else
        res$residuals

    # influence function
    infl <- if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
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
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
# Tukey biweight M-estimator of the weighted total (depends on pkg survey)
svytotal_tukey <- function(x, design, k, type = "rwm", na.rm = FALSE,
                           verbose = TRUE, ...)
{
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Tukey M-estimator (type = ", type, ")"), dat$design,
            match.call(), "mest", type = type, psi = 2, psi_fun = "Tukey",
            k = k))

    # population- vs. domain-level estimate (matters for calibrated designs)
    res <- weighted_total_tukey(dat$y, dat$w, k, type, TRUE, FALSE, verbose,
                                ...)

    # influence function
    infl <- if (dat$calibrated) {
        tmp <- numeric(length(dat$in_domain))
        tmp[dat$in_domain] <- res$robust$robweights * res$model$w * res$model$y
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
    res$design <- dat$design
    res$call <- match.call()
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
