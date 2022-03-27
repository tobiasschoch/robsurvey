# Dalen the weighted mean (depends on pkg survey)
svymean_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    dat <- .checkformula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "mean",
            estimator = list(string = paste0("Dalen ", type,
                " estimator (censored at ", censoring, ")"),
                censoring = censoring),
            estimate = stats::setNames(NA, dat$yname), variance = NA,
            residuals = NA, model = NA, design = design, call = match.call()),
            class = "svystat_rob"))
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
    class(res) <- c("svystat_rob")
    res
}
# Dalen the weighted total (depends on pkg survey)
svytotal_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = TRUE)
{
    dat <- .checkformula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "total",
            estimator = list(string = paste0("Dalen ", type,
                " estimator (censored at ", censoring, ")"),
                censoring = censoring),
            estimate = stats::setNames(NA, dat$yname), variance = NA,
            residuals = NA, model = NA, design = design, call = match.call()),
            class = "svystat_rob"))
    # otherwise
    design <- dat$design
    res <- weighted_total_dalen(dat$y, dat$w, censoring, type, TRUE, FALSE,
        verbose)
    # compute variance
    infl <- res$robust$xw
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob")
    res
}
