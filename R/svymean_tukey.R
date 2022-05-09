# Tukey biweight M-estimator of the weighted mean (depends on pkg survey)
svymean_tukey <- function(x, design, k, type = "rhj", na.rm = FALSE,
    verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("mean", dat$yname,
            paste0("Tukey M-estimator (type = ", type, ")"), match.call(),
            design, "mest", type = type, psi = 2, psi_fun = "Tukey", k = k))
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
    class(res) <- c("svystat_rob", "mer_capable", "mest")
    res
}
# Tukey biweight M-estimator of the weighted total (depends on pkg survey)
svytotal_tukey <- function(x, design, k, type = "rhj", na.rm = FALSE,
        verbose = TRUE, ...)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)
    dat <- .check_formula(x, design, na.rm)
    # in the presence of NA's
    if (dat$failure)
        return(.new_svystat_rob("total", dat$yname,
            paste0("Tukey M-estimator (type = ", type, ")"), match.call(),
            design, "mest", type = type, psi = 2, psi_fun = "Tukey", k = k))
    # otherwise
    design <- dat$design
    res <- weighted_total_tukey(dat$y, dat$w, k, type, TRUE, FALSE, verbose,
        ...)
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
