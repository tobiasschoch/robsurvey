svymean_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
{
    dat <- .checkformula(x, design)
    res <- weighted_mean_trimmed(dat$y, dat$w, LB, UB, info = TRUE, na.rm)
    # influence function
    infl <- .infl_trimmed(dat$y, dat$w, LB, UB, res$estimate)
    # variance
    infl <- infl * dat$w / sum(dat$w)
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- "svystat_rob"
    res
}

svytotal_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
{
    res <- svymean_trimmed(x, design, LB, UB, na.rm)
    sum_w <- sum(res$model$w)
    res$estimate <- res$estimate * sum_w
    res$variance <- res$variance * sum_w^2
    res$characteristic <- "total"
    res$call <- match.call()
    res
}
