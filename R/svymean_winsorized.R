# weighted winsorized mean (depends on pkg survey)
svymean_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB,
    na.rm = FALSE, simple_var = FALSE)
{
    if (!is.language(x))
        stop("Argument 'x' must be a formula object\n", call. = FALSE)

    dat <- .checkformula(x, design)
    res <- weighted_mean_winsorized(dat$y, dat$w, LB, UB, info = TRUE, na.rm)
    # influence function
    infl <- if (simple_var)
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
# weighted one-sided k winsorized mean (depends on pkg survey)
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
# weighted winsorized total (depends on pkg survey)
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
# weighted one-sided k winsorized total (depends on pkg survey)
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
# influence function, Huber (1981, p. 58-59)
.infl_winsorized <- function(x, w, LB, UB, wm, ngrid = 401)
{
    x <- as.double(x); w <- as.double(w)
    qs <- weighted_quantile(x, w, probs = c(LB, UB))
    # density estimates at 'LB' and '1 - UB'
    bwd <- KernSmooth::dpik(x, scalest = "minim", level = 2, kernel = "normal",
        canonical = FALSE, gridsize = ngrid, range.x = range(x),
        truncate = TRUE)
    at <- seq(min(x), max(x), length = ngrid)
    nx <- rowsum(c(rep(0, ngrid), w), c(1:ngrid, findInterval(x, at)))
    dens <- KernSmooth::locpoly(rep(1, ngrid), nx * ngrid / (diff(range(x)) *
        sum(w)), binned = TRUE, bandwidth = bwd, range.x = range(x))
    f_LB <- dens$y[min(which(dens$x >= qs[1]))]
    f_UB <- dens$y[max(which(dens$x <= qs[2]))]
    # influence function
    infl <- pmin.int(qs[2] + (1 - UB) / f_UB, pmax.int(qs[1] - LB / f_LB, x))
    # functional W corresponding to winsorized mean
    W <- unname((UB - LB) * wm + LB * qs[1] + (1 - UB) * qs[2])
    infl <- infl - W
    if (f_LB > 0)
        infl <- infl + LB^2 / f_LB
    if (f_UB > 0)
        infl <- infl + (1 - UB)^2 / f_UB
    infl
}
