# Regression estimator of the mean (depends on pkg survey)
svymean_reg <- function(object, auxiliary, k = Inf,
    psi = c("Huber", "Huberasym", "Tukey"), check.names = TRUE)
{
    stopifnot(k > 0)
    w <- object$model$w; sum_w <- sum(w)
    # check whether auxiliary matches with the regression estimate
    x <- .checkauxiliary(object, auxiliary, "mean", N = NULL,
        check.names, na.action = stats::na.omit)
    # greg estimate
    resid_winsorized <- .psi_function(object$residuals, k, psi)
    est <- sum(auxiliary * object$estimate) +
        sum(w * resid_winsorized) / sum_w
    names(est) <- object$model$yname
    # compute variance
    design <- object$design
    v <- survey::svyrecvar(resid_winsorized * w / sum_w, design$cluster,
        design$strata, design$fpc, postStrata = design$postStrata)
    res <- list(characteristic = "mean", estimator = "GREG M-estimator",
        estimate = est, variance = v, design = design, call = match.call())
    class(res) <- "svystat_rob"
    res
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
# svymean_reg(m, workplace)
#


    # if (!inherits(object, "svyreg_rob"))
    #     stop("svymean_reg is not defined for this object\n", call. = FALSE)
    #
    # w <- object$model$w
    # sum_w <- sum(w)
    #
    # # FIXME: argument N
    # x <- .checkauxiliary(object, auxiliary, "mean", N = NULL,
    #     check.names, na.action = stats::na.omit)
    #
    # # GREG estimate
    # est <- sum(x * object$estimate)
    # names(est) <- object$model$yname
    # # FIXME: more general GREG condition
    # if (object$model$intercept == 0)
    #     est <- est +  sum(w * object$residuals) / sum_w
    #
    # # variance estimate  # NOTE:: check
    # infl <- object$robust$robweights * object$residuals * w / sum_w
    # design <- object$design
    # v <- survey::svyrecvar(infl, design$cluster, design$strata, design$fpc,
    #     postStrata = design$postStrata)
    #
    # object$characteristic <- "mean"
    # object$estimate <- est
    # object$variance <- v
    # object$call <- match.call()
    # class(object) <- "svystat_rob"
    # object
