# GREG predictor of the total
svytotal_reg <- function(object, totals, N = NULL, type, k = NULL,
                         check.names = TRUE, keep_object = TRUE, ...)
{
    if (!inherits(object, "svyreg_rob"))
        stop(paste0("The function cannot be used for an object of class '",
                    class(object), "'\n"))
    model <- object$model
    # check whether 'type' and 'k' are correctly specified
    type <- .check_k(type, k)
    # check whether auxiliary matches with the regression estimate
    auxiliary <- .check_auxiliary(object, totals, N, check.names, FALSE)
    # g-weights
    gi <- .gi_weights(object, auxiliary, type, k)
    # estimate
    est <- sum(gi * model$y)
    names(est) <- model$yname
    # residuals (not accounting for heteroscedasticity)
    ri <- model$y - as.numeric(model$x %*% object$estimate)
    # variance estimate
    design <- object$design
    v <- svyrecvar(ri * gi, design$cluster, design$strata, design$fpc,
                   postStrata = design$postStrata)
    # return
    model$coef <- object$estimate
    model$call <- object$call
    res <- structure(list(characteristic = "total",
        estimator = .method_name(type, k), estimate = est,
        robust = list(robweights = .bi_weights(object, type, k)),
        residuals = ri, model = model, design = design,
        call = match.call(), variance = v, auxiliary = auxiliary,
        gweights = gi), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    # return
    res
}
# GREG predictor of the mean
svymean_reg <- function(object, totals, N = NULL, type, k = NULL,
                        check.names = TRUE, keep_object = TRUE,
                        N_unknown = FALSE, ...)
{
    if (!inherits(object, "svyreg_rob"))
        stop(paste0("The function cannot be used for an object of class '",
                    class(object), "'\n"))
    model <- object$model
    # check whether 'type' and 'k' are correctly specified
    type <- .check_k(type, k)
    # if N is unknown, it will be estimated
    if (N_unknown)
        N <- sum(model$w)

    # check whether auxiliary matches with the regression estimate
    auxiliary <- .check_auxiliary(object, totals, N, check.names, TRUE)
    # g-weights
    gi <- .gi_weights(object, auxiliary, type, k) / N
    # estimate
    est <- sum(gi * model$y)
    names(est) <- model$yname
    # residuals (not accounting for heteroscedasticity)
    ri <- model$y - as.numeric(model$x %*% object$estimate)
    # variance estimate
    design <- object$design
    v <- svyrecvar(ri * gi, design$cluster, design$strata, design$fpc,
                   postStrata = design$postStrata)
    # return
    model$coef <- object$estimate
    model$call <- object$call
    res <- structure(list(characteristic = "mean",
        estimator = .method_name(type, k), estimate = est,
        robust = list(robweights = .bi_weights(object, type, k)),
        residuals = ri, model = model, design = design,
        call = match.call(), variance = v, auxiliary = auxiliary,
        gweights = gi), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    # return
    res
}
# g-weights
.gi_weights <- function(object, auxiliary, type, k)
{
    # bi's (of QR-predictor; originally, ri's in Wright, 1983)
    bi <- .bi_weights(object, type, k) * object$model$w
    delta <- auxiliary - colSums(bi * object$model$x)
    # qi's (of QR-predictor)
    qi_sqrt <- sqrt(.qi_weights(object))
    QR_x <- qr(qi_sqrt * object$model$x)
    H <- backsolve(qr.R(QR_x), t(qr.Q(QR_x)))
    # g-weights
    bi + as.numeric(crossprod(H, delta)) * qi_sqrt
}
# qi's in Wright's QR-estimator
.qi_weights <- function(object)
{
    # robustness and design weights
    qi <- robweights(object) * object$model$w
    # account for heteroscedasticity
    if (!is.null(object$model$var))
        qi <- qi / object$model$var
    # return
    qi
}
# bi's without sampling weight (actually, ri's in Wright's QR-estimator)
.bi_weights <- function(object, type, k)
{
    n <- object$model$n
    # scaled residuals (already accounts for heteroscedasticity)
    ri <- object$residuals / object$scale
    # Schweppe type GM-estimator
    if (object$estimator$type == 2)
        ri <- ri / object$model$xwgt
    # robustness weights (with k used in prediction)
    bi <- switch(type,
        "projective" = rep(0, n),
        "ADU" = rep(1, n),
        "huber" = .psi_wgt_function(ri, k, "Huber"),
        "tukey" = .psi_wgt_function(ri, k, "Tukey"),
        "lee" = rep(k, n),
        "BR" = .psi_wgt_beaumont_rivest(ri, object$model$w, k),
        "duchesne" = .psi_duchesne(ri, k[1], k[2]) / ri)
    # multiplicative factor for Mallows type GM-estimator
    if (object$estimator$type == 1)
        bi <- bi * object$model$xwgt
    # return
    bi
}
# Modified Huber psi-function of Duchsene (1999)
.psi_duchesne <- function(x, a, b)
{
    pmin.int(a, pmax.int(-a, x)) +
        (x - pmin.int(a / b, pmax.int(-a / b, x))) * b
}
# Modified Huber psi-function of Beaumont and Rivest (2009)
.psi_wgt_beaumont_rivest <- function(x, w, k)
{
    (x / w + (1 - 1 / w) * .psi_function(x, k, "Huber")) / x
}
# Check whether 'type' and 'k' are correctly specified
.check_k <- function(type, k)
{
    if (missing(type))
        stop("Argument 'type' is missing\n", call. = FALSE)
    type <- match.arg(type, c("projective", "ADU", "huber", "tukey", "lee",
                              "BR", "duchesne"))
    switch(type,
        "ADU" = {
            if (!is.null(k))
                warning("Argument 'k' is ignored\n", call. = FALSE,
                        immediate. = TRUE)
        },
        "projective" = {
            if (!is.null(k))
                warning("Argument 'k' is ignored\n", call. = FALSE,
                        immediate. = TRUE)
        },
        "huber" = stopifnot(is.numeric(k), length(k) == 1, k > 0),
        "tukey" = stopifnot(is.numeric(k), length(k) == 1, k > 0),
        "lee" = stopifnot(is.numeric(k), length(k) == 1, k >= 0, k <= 1),
        "BR" = stopifnot(is.numeric(k), length(k) == 1, k > 0),
        "duchesne" = stopifnot(is.numeric(k),
            "k must be a 2-vector of positive real values" = length(k) == 2,
            all(k > 0)))
    type
}
# Name of the method
.method_name <- function(type, k)
{
    string <- switch(type,
        "projective" = "(projective)",
        "ADU" = "(ADU)",
        "huber" = paste0("(robust, Huber psi, k = ", k, ")"),
        "tukey" = paste0("(robust, Tukey psi, k = ", k, ")"),
        "lee" = paste0("robust (Lee, k = ", k, ")"),
        "BR" = paste0("(BR, Huber psi, k = ", k, ")"),
        "duchesne" = paste0("(Duchesne, k = (", k[1], ", ", k[2], "))"))
    if (type == "tukey") {
        psi <- 2
        psi_fun <- "Tukey"
    } else {
        psi <- 0
        psi_fun <- "Huber"
    }
    list(string = paste0("GREG predictor ", string), type = type, psi = psi,
         psi_fun = psi_fun, k = k)
}
# check auxiliary totals and means
.check_auxiliary <- function(object, totals, N = NULL, check.names, mean)
{
    has_intercept <- attr(object$terms, "intercept")
    name_totals <- names(totals)
    name_coef <- names(object$estimate)
    if (has_intercept)
        name_coef <- name_coef[-1]

    if (is.matrix(totals) || is.data.frame(totals)) {
        if (all(dim(totals) > 1))
            stop("Argument 'totals' must be a vector\n", call. = FALSE)
        totals <- as.numeric(totals)
    }

    if (is.null(N)) {
        if (mean)
            stop("Argument 'N' is missing\n", call. = FALSE)
        else if (has_intercept)
            stop("Argument 'N' is missing (model has an intercept)\n",
                call. = FALSE)
    } else {
        stopifnot(is.numeric(N), N > 0)
    }

    if (length(totals) != length(name_coef))
        stop("Length of vector of totals is not appropriate\n", call. = FALSE)

    if (check.names && !is.null(name_totals)) {
        if (!all(name_totals == name_coef))
            stop("Variable names do not match (check.names = TRUE)",
                 call. = FALSE)
    }

    if (has_intercept)
        c(N, totals)
    else
        totals
}
