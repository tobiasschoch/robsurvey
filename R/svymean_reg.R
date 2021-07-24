# Regression estimator of the mean (depends on pkg survey)
svymean_reg <- function(object, auxiliary, check.names = TRUE)
{
    if (!inherits(object, "svyreg_rob"))
        stop("svymean_reg is not defined for this object\n", call. = FALSE)

    w <- object$model$w
    sum_w <- sum(w)

    # FIXME: argument N
    x <- .checkauxiliary(object, auxiliary, "mean", N = NULL,
        check.names, na.action = stats::na.omit)

    # GREG estimate
    est <- sum(x * object$estimate)
    names(est) <- object$model$yname
    # FIXME: more general GREG condition
    if (object$model$intercept == 0)
        est <- est +  sum(w * object$residuals) / sum_w

    # variance estimate  # NOTE:: check
    infl <- object$robust$robweights * object$residuals * w / sum_w
    design <- object$design
    v <- survey::svyrecvar(infl, design$cluster, design$strata, design$fpc,
        postStrata = design$postStrata)

    object$characteristic <- "mean"
    object$estimate <- est
    object$variance <- v
    object$call <- match.call()
    class(object) <- "svystat_rob"
    object
}

svytotal_reg <- function()
{

}


#
# library(robsurvey)
# data(workplace)
# library(survey)
# dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
#     data = workplace)
# m <- svyreg_huber(payroll ~ employment, dn, k = 8)
# svymean_reg(m, workplace, check.names = TRUE)
#
#
