# weighted trimmed mean
svymean_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)

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

# weighted trimmed total
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

# influence function, Huber (1981, p. 58)
.infl_trimmed <- function(x, w, LB, UB, tm)
{
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    x_wins <- pmin.int(qs[2], pmax.int(qs[1], x))
    # functional W corresponding to winsorized mean
    W <- (UB - LB) * tm + LB * qs[1] + (1 - UB) * qs[2]
    # return
    (x_wins - W) / (UB - LB)
}
