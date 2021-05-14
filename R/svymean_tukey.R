# Tukey biweight M-estimator of the weighted mean (depends on pkg survey)
svymean_tukey <- function(x, design, k, type = "rwm", na.rm = FALSE,
    verbose = TRUE, ...)
{
    dat <- .checkformula(x, design)
    res <- weighted_mean_tukey(dat$y, dat$w, k, type, asym, info = TRUE,
        na.rm, verbose, ...)
    # modify residuals for type 'rht' (only for variance estimation)
    r <- if (type == "rht")
        sqrt(res$model$var) * res$model$y - res$estimate
    else
        res$residuals
   # compute variance
    infl <- res$robust$robweights * r * dat$w / sum(dat$w)
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
# Tukey biweight M-estimator of the weighted total (depends on pkg survey)
svytotal_tukey <- function(x, design, k, type = "rwm", na.rm = FALSE,
        verbose = TRUE, ...)
{
    dat <- robsurvey:::.checkformula(x, design)
    res <- weighted_total_tukey(dat$y, dat$w, k, type, asym, info = TRUE,
        na.rm, verbose, ...)
    # compute variance
    infl <- res$robust$robweights * dat$y * dat$w
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob", "mer_capable")
    res
}
