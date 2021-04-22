svymean_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB,
    na.rm = FALSE, simple_var = FALSE)
{
    dat <- .checkformula(x, design)
    res <- weighted_mean_winsorized(dat$y, dat$w, LB, UB, info = TRUE, na.rm)
    # influence function
    infl if (simple_var)
        .infl_trimmed(dat$y, dat$w, LB, UB, res$estimate)
    else
        .infl_winsorized(dat$y, dat$w, LB, UB, res$estimate)
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

svymean_k_winsorized <- function(x, design, k, na.rm = FALSE,
    simple_var = FALSE)
{
    dat <- .checkformula(x, design)
    res <- weighted_mean_k_winsorized(dat$y, dat$w, k, info = TRUE, na.rm)
    w <- res$model$w; x <- res$model$y
    # influence function
    infl <- if (simple_var)
        .infl_trimmed(dat$y, dat$w, 0, res$robust$UB, res$estimate)
    else
        .infl_winsorized(dat$y, dat$w, 0, res$robust$UB, res$estimate)
    # variance
    infl <- infl * w / sum(w)
    res$variance <- survey::svyrecvar(infl, design$cluster, design$strata,
        design$fpc, postStrata = design$postStrata)
    names(res$estimate) <- dat$yname
    res$call <- match.call()
    res$design <- design
    class(res) <- "svystat_rob"
    res
}

svytotal_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB,
    na.rm = FALSE, simple_var = FALSE)
{
    res <- svymean_winsorized(x, design, LB, UB, na.rm, simple_var)
    sum_w <- sum(res$model$w)
    res$estimate <- res$estimate * sum_w
    res$variance <- res$variance * sum_w^2
    res$characteristic <- "total"
    res$call <- match.call()
    res
}

svytotal_k_winsorized <- function(x, design, k, na.rm = FALSE,
    simple_var = FALSE)
{
    res <- svymean_k_winsorized(x, design, k, na.rm, simple_var)
    sum_w <- sum(res$model$w)
    res$estimate <- res$estimate * sum_w
    res$variance <- res$variance * sum_w^2
    res$characteristic <- "total"
    res$call <- match.call()
    res
}
