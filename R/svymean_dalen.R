# Dalen the weighted mean (depends on pkg survey)
svymean_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname, paste0("Dalen ", type,
            " estimator (censored at ", censoring, ")"), match.call(),
            design, "dalen", censoring = censoring))
    # otherwise
    design <- dat$design
    res <- weighted_mean_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
        verbose)
    # compute variance
    infl <- (res$robust$xw - dat$w * res$estimate) / sum(dat$w)
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob", "dalen")
    res
}
# Dalen the weighted total (depends on pkg survey)
svytotal_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname, paste0("Dalen ", type,
            " estimator (censored at ", censoring, ")"), match.call(),
            design, "dalen", censoring = censoring))
    # otherwise
    design <- dat$design
    res <- weighted_total_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
        verbose)
    # compute variance
    res$variance <- survey::svyrecvar(res$robust$xw, design$cluster,
        design$strata, design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob", "dalen")
    res
}
