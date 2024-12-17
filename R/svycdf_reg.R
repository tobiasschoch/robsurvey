svycdf_reg <- function(object, at, non_sampled, type = "cd", ...)
{
    stopifnot(is.numeric(at))
    if (!inherits(object, "svyreg_rob"))
        stop(paste0("The function cannot be used for an object of class '",
                    class(object), "'\n"))

    # std. residuals of the sample-based fit
    res_sample <- if (is.null(object$model$var))
            object$residuals
        else
            object$residuals / sqrt(object$model$var)

    # linear predictor for the non-sampled part
    f <- object$call$formula; f[[2]] <- NULL
    linpred_nonsample <- as.numeric(tcrossprod(object$estimate,
        model.matrix(terms.formula(f), non_sampled)))

    # standard deviations (heteroscedasticity) for the non-sampled part
    sd_nonsample <- if (is.null(object$model$var))
            rep(1, NROW(non_sampled))
        else
            sqrt(as.numeric(model.matrix(object$model$var_terms, non_sampled)))

    # estimate of ECDF
    tmp <- .C(C_ecdf_cd, rs = as.double(res_sample),
        linpred_nonsample = as.double(linpred_nonsample),
        sd_nonsample = as.double(sd_nonsample), at = as.double(at),
        at_length = as.integer(length(at)), n = as.integer(length(res_sample)),
        N = as.integer(length(sd_nonsample)),
        result = as.double(rep(0, length(at))))
    tmp$result
}
