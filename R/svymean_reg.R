# Regression estimator of the mean (depends on pkg survey)
svymean_reg <- function(object, auxiliary, k = Inf,
    psi = c("Huber", "Huberasym", "Tukey"), check.names = TRUE)
{
    model <- object$model
    stopifnot(k > 0)
    w <- model$w
    # check whether auxiliary matches with the regression estimate
    x <- .checkauxiliary(object, auxiliary, "mean", check.names,
        na.action = stats::na.omit)

    # greg estimate
    mod_resid <- w * .psi_function(object$residuals, k, psi) / sum(w)
    est <- sum(x * object$estimate)
    # bias-correction term
    if (!is.null(model$var)) {

        qr(cbind(model$x, model$var))$rank
        est <- est + sum(mod_resid)
    } else if (model$intercept) {

    }
    names(est) <- model$yname

    # compute variance
    design <- object$design
    v <- survey::svyrecvar(mod_resid, design$cluster,
        design$strata, design$fpc, postStrata = design$postStrata)
    res <- list(characteristic = "mean", estimator = "GREG M-estimator",
        estimate = est, variance = v, design = design, call = match.call())
    class(res) <- "svystat_rob"
    res
}

svytotal_reg <- function()
{
    # mean * N_hat?
}

if(0){
#
library(robsurvey)
data(workplace)
library(survey)
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
m <- svyreg_huber(payroll ~ employment, dn, k = Inf)
m0 <- svyreg(payroll ~ -1 + employment, dn)

svymean_reg(m, workplace)

Xmean <- coef(svymean(~employment, dn))
Xmean <- mean(workplace$employment)
svymean_reg(m, c(1, Xmean))

}

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
