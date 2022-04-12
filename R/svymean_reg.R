svymean_reg <- function(object, auxiliary, k = Inf,
    psi = c("Huber", "Huberasym", "Tukey"), check.names = TRUE, na.rm = FALSE)
{
    if (class(object) != "svyreg_rob")
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))
    stopifnot(k > 0)
    model <- object$model
    model$coefficients <- object$estimate
    # check whether auxiliary matches with the regression estimate
    x <- .checkauxiliary(object, auxiliary, "mean", check.names,
            na.rm = na.rm)
    # greg estimate
    est <- sum(x * object$estimate)
    # linearized values
    ui <- .psi_wgt_function(object$residuals / object$scale, k, psi)
    zi <- model$w * ui * object$residuals / sum(model$w)
    # bias-correction term
    est <- est + sum(zi)
    names(est) <- as.character(attributes(object$terms)$variables[[2]])
    # compute variance
    design <- object$design
    v <- survey::svyrecvar(zi, design$cluster, design$strata, design$fpc,
        postStrata = design$postStrata)
    # prepare return value
    string <- ifelse(is.infinite(k), "GREG estimator",
        paste0("Robust GREG (", psi, " psi-function, k = ", k, ")"))
    return(structure(list(characteristic = "mean",
        estimator = list(string = string, k = k), estimate = est,
        robust = list(robweights = ui), residuals = model$y - est,
        variance = v, model = model, design = design, call = match.call()),
        class = "svystat_rob"))
}
