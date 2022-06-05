# robust GREG for mean
svymean_reg <- function(object, auxiliary, type, k = NULL, check.names = TRUE,
    na.rm = FALSE, keep_object = TRUE, verbose = TRUE)
{
    if (verbose)
        message("Note: This functions is experimental!")

    if (!inherits(object, "svyreg_rob"))
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))
    # check whether 'type' and 'k' are correctly specified
    type <- .check_k(type, k)
    # standardized residuals
    r_std <- .std_residuals(object)
    # check whether auxiliary matches with the regression estimate
    x <- .check_auxiliary(object, auxiliary, "mean", check.names, na.rm = na.rm)
    # linear predictor
    est <- sum(x * object$estimate)
    names(est) <- object$model$yname
    # GREG bias correction
    tmp <- .bias_correction(object, r_std, type, k)
    zi <- tmp$zi / sum(object$model$w)
    correction <- tmp$correction / sum(object$model$w)
    est <- est + correction
    # variance estimate
    design <- object$design
    v <- survey::svyrecvar(zi, design$cluster, design$strata, design$fpc,
        postStrata = design$postStrata)
    res <- structure(list(characteristic = "mean",
        estimator = list(string = .method_name(type, k),
        type = type, psi = 0, psi_fun = "Huber", k = k), estimate = est,
        robust = list(robweights = tmp$ui), residuals = r_std,
        model = list(object$model, coef = object$estimate,
        auxiliary = auxiliary), design = design, call = match.call(),
        variance = v, bias = correction), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    res
}
# robust GREG for total
svytotal_reg <- function(object, auxiliary, type, k = NULL, check.names = TRUE,
    na.rm = FALSE, keep_object = TRUE, verbose = TRUE)
{
    if (verbose)
        message("Note: This functions is experimental!")

    if (!inherits(object, "svyreg_rob"))
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))
    # check whether 'type' and 'k' are correctly specified
    type <- .check_k(type, k)
    # standardized residuals
    r_std <- .std_residuals(object)
    # check whether auxiliary matches with the regression estimate
    x <- .check_auxiliary(object, auxiliary, "total", check.names,
        na.rm = na.rm)
    # linear predictor
    est <- sum(x * object$estimate)
    names(est) <- object$model$yname
    # GREG bias correction
    tmp <- .bias_correction(object, r_std, type, k)
    zi <- tmp$zi
    est <- est + tmp$correction
    # variance estimate
    design <- object$design
    v <- survey::svyrecvar(zi, design$cluster, design$strata, design$fpc,
        postStrata = design$postStrata)
    res <- structure(list(characteristic = "total",
        estimator = list(string = .method_name(type, k),
        type = type, psi = 0, psi_fun = "Huber", k = k), estimate = est,
        robust = list(robweights = tmp$ui), residuals = r_std,
        model = list(object$model, coef = object$estimate,
        auxiliary = auxiliary), design = design, call = match.call(),
        variance = v, bias = tmp$correction), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    res
}
# check auxiliary totals and means
.check_auxiliary <- function(object, data, est = "mean", check.names = TRUE,
    na.rm = FALSE)
{
    names_beta <- names(object$estimate)
    names_data <- names(data)
    if (is.vector(data))
        if (NCOL(data) < NROW(data))
	        data <- t(data)

    # vector of population x-means or -totals
    if (NROW(data) == 1) {
        # drop intercept (if there is one)
        if (attr(object$terms, "intercept")) {
            names_data <- names_data[-1]
            names_beta <- names_beta[-1]
        }
        if (length(data) != object$model$p)
            stop("Length of auxiliary data does not match\n", call. = FALSE)

        if (check.names && !is.null(names_data)) {
            if (!all(names_data == names_beta))
                stop("Variable names do not match (check.names = TRUE)",
                    call. = FALSE)
        }
        x <- data
    # compute mean/total based on design matrix
    } else {
        mf <- stats::model.frame(object$call$formula, data,
            na.action = stats::na.pass)
        xmat <- stats::model.matrix(stats::terms(mf), mf)
        # check for missing values
        is_not_NA <- stats::complete.cases(xmat)
        if (!all(is_not_NA) && na.rm) {
            xmat <- xmat[is_not_NA, ]
        }
        x <- switch(est,
            "mean" = colMeans(xmat),
            "total" = colSums(xmat))
    }
    unname(x)
}
# Check whether 'type' and 'k' are correctly specified
.check_k <- function(type, k)
{
    type <- match.arg(type, c("projective", "ADU", "robust", "lee", "BR",
        "duchesne"))
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
        "robust" = stopifnot(is.numeric(k), length(k) == 1, k > 0),
        "lee" = stopifnot(is.numeric(k), length(k) == 1, k >= 0, k <= 1),
        "BR" = stopifnot(is.numeric(k), length(k) == 1, k > 0),
        "duchesne" = stopifnot(is.numeric(k),
            "k must be a 2-vector of positive real values" = length(k) == 2,
            all(k > 0)))
    type
}
# Standardized residuals
.std_residuals <- function(object)
{
    # scale
    r_std <- object$residuals / object$scale
    # account for heteroscedasticity
    if (!is.null(object$model$var))
        r_std <- r_std / sqrt(object$model$var)
    # Schweppe type GM-estimator
    if (object$estimator$type == 2)
        r_std <- r_std / object$model$xwgt
    r_std
}
# Name of the method
.method_name <- function(type, k)
{
    switch(type,
        "projective" = "projective",
        "ADU" = "ADU",
        "robust" = paste0("robust, Huber psi, k = ", k),
        "lee" = paste0("Lee, k = ", k),
        "BR" = paste0("BR, Huber psi, k = ", k),
        "duchesne" = paste0("Duchesne, k = (", k[1], ", ", k[2], ")"))
}
# GREG bias correction
.bias_correction <- function(object, r_std, type, k)
{
    w <- object$model$w
    n <- object$model$n
    # robustness weights
    ui <- switch(type,
        "projective" = rep(1, n),
        "ADU" = rep(1, n),
        "robust" = .psi_wgt_function(r_std, k, "Huber"),
        "lee" = rep(k, n),
        "BR" = .psi_wgt_beaumont_rivest(r_std, w, k),
        "duchesne" = .psi_duchesne(r_std, k[1], k[2]) / r_std)
    # multiplicative factor for Mallows type GM-estimator
    if (object$estimator$type == 1)
        ui <- ui * object$model$xwgt
    # linearization (influence function)
    zi <- w * ui * object$residuals
    # bias-correction term of GREG
    correction <- ifelse(type == "projective", 0, sum(zi))
    list(zi = zi, correction = correction, ui = ui)
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
