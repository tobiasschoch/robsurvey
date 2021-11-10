# Tukey biweight M-estimator of the weighted mean (depends on pkg survey)
svymean_tukey <- function(x, design, k, type = "rhj", na.rm = FALSE,
    verbose = TRUE, ...)
{
    dat <- .checkformula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "mean",
            estimator = list(string = paste0("Tukey M-estimator (type = ",
                type, ")"), type = type, psi = 2, psi_fun = "Tukey", k = k),
            estimate = stats::setNames(NA, dat$yname), variance = NA,
            residuals = NA, model = NA, design = design, call = match.call()),
            class = "svystat_rob"))
    # otherwise
    design <- dat$design
    res <- weighted_mean_tukey(dat$y, dat$w, k, type, TRUE, FALSE, verbose,
        ...)
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
svytotal_tukey <- function(x, design, k, type = "rhj", na.rm = FALSE,
        verbose = TRUE, ...)
{
    dat <- .checkformula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(structure(list(characteristic = "mean",
            estimator = list(string = paste0("Tukey M-estimator (type = ",
                type, ")"), type = type, psi = 2, psi_fun = 2, k = k),
            estimate = stats::setNames(NA, dat$yname), variance = NA,
            residuals = NA, model = NA, design = design, call = match.call()),
            class = "svystat_rob"))
    design <- dat$design
    res <- weighted_total_tukey(dat$y, dat$w, k, type, TRUE, FALSE, verbose,
        ...)
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
