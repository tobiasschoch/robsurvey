# Dalen the weighted mean (depends on pkg survey)
svymean_dalen <- function(x, design, censoring, type = "Z2", na.rm = FALSE,
    verbose = FALSE)
{
    dat <- .checkformula(x, design)
    res <- weighted_mean_dalen(dat$y, dat$w, censoring, type, info = TRUE,
        na.rm, verbose)
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
    verbose = FALSE)
{
    dat <- .checkformula(x, design)
    res <- weighted_total_dalen(dat$y, dat$w, censoring, type, info = TRUE,
        na.rm, verbose)
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
