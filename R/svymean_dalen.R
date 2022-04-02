# Dalen the weighted mean (depends on pkg survey)
svymean_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .checkformula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.empty_svystat_rob("mean", dat$yname, paste0("Dalen ", type,
            " estimator (censored at ", censoring, ")"), match.call(),
            design, censoring = censoring))
    # otherwise
    design <- dat$design
    res <- weighted_mean_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
        verbose)
    # compute variance
    infl <- res$robust$xw / sum(dat$w)
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- "svystat_rob"
    res
}
# Dalen the weighted total (depends on pkg survey)
svytotal_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    res <- svymean_dalen(x, design, censoring, type, na.rm, verbose)
    res$call <- match.call()
    res$characteristic <- "total"
    if (is.na(res$estimate))
        return(res)
    sum_w <- sum(res$model$w)
    res$estimate <- res$estimate * sum_w
    res$variance <- res$variance * sum_w^2
    res
}
