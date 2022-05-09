# Huber M-estimator of the weighted mean (depends on pkg survey)
svymean_huber <- function(x, design, k, type = "rhj", asym = FALSE,
    na.rm = FALSE, verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Huber M-estimator (type = ", type,
            ifelse(asym, "; asym. psi", ""), ")"), match.call(),
            design, "mest", type = type, psi = ifelse(asym, 1, 0),
            psi_fun = "Huber", k = k))
    # otherwise
    design <- dat$design
    res <- weighted_mean_huber(dat$y, dat$w, k, type, asym, TRUE, FALSE,
        verbose, ...)
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
    class(res) <- c("svystat_rob", "mer_capable", "mest")
    res
}
# Huber M-estimator of the weighted total (depends on pkg survey)
svytotal_huber <- function(x, design, k, type = "rhj", asym = FALSE,
    na.rm = FALSE, verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Huber M-estimator (type = ", type,
            ifelse(asym, "; asym. psi", ""), ")"), match.call(),
            design, "mest", type = type, psi = ifelse(asym, 1, 0),
            psi_fun = "Huber", k = k))
    # otherwise
    design <- dat$design
    res <- weighted_total_huber(dat$y, dat$w, k, type, asym, TRUE, FALSE,
        verbose, ...)
    # modify residuals for type 'rht' (only for variance estimation)
    r <- if (type == "rht")
        sqrt(res$model$var) * res$model$y - res$estimate
    else
        res$residuals
   # compute variance
    infl <- res$robust$robweights * dat$y * dat$w
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- c("svystat_rob", "mer_capable", "mest")
    res
}
