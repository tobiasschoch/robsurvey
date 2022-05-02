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
# robust GREG for total
#FIXME: Mallows + Schweppe est.
svytotal_reg <- function(object, auxiliary, type, k = NULL, check.names = TRUE,
    na.rm = FALSE, keep_object = TRUE)
{
    type <- match.arg(type, c("projective", "ADU", "robust", "lee", "BR",
        "duchesne"))
    if (is.null(k) && type %in% c("robust", "lee", "BR", "duchesne"))
        stop("Robustness tuning must not be NULL\n")
    if (class(object) != "svyreg_rob")
        stop(paste0("The function cannot be used for an object of class '",
            class(object), "'\n"))
    if (type == "lee" && (k < 0 || k > 1))
        stop("For type = 'lee', argument 'k' must satisfy 0 <= k <= 1\n")
    if (type %in% c("robust", "BR", "duchesne") && k <= 0)
        stop("Argument 'k' must be positive\n")
    if (type == "duchesne" && length(k) != 2)
        stop("For type = 'duchesene', argument 'k' must be a vector: c(a, b)\n")
    # standardized residuals
    r <- object$residuals / object$scale
    # account for heteroscedasticity
    if (!is.null(object$model$var))
        r <- r / sqrt(object$model$var)
    # Schweppe type GM-estimator
    if (object$estimator$type == 2)
        r <- r / object$model$xwgt
    # check whether auxiliary matches with the regression estimate
    x <- robsurvey:::.checkauxiliary(object, auxiliary, "total", check.names,
            na.rm = na.rm)
    # linear predictor
    est <- sum(x * object$estimate)
    names(est) <- as.character(attributes(object$terms)$variables[[2]])
    # bias correction factor
    w <- object$model$w
    switch(type,
        "projective" = {
            bi <- rep(0, object$model$n)
            string <- "projective"
        },
        "ADU" = {
            bi <- w
            string <- "ADU"
        },
        "robust" = {
            bi <- w * .psi_wgt_function(r, k, "Huber")
            string <- paste0("robust, Huber psi, k = ", k)
        },
        "lee" = {
            bi <- k * w
            string <- paste0("Lee, k = ", k)
        },
        "BR" = {
            bi <- w * .psi_wgt_beaumont_rivest(r, w, k)
            string <- paste0("BR, Huber psi, k = ", k)

        },
        "duchesne" = {
            bi <- .psi_duchesne(r, k[1], k[2]) / r
            string <- paste0("Duchesne, k = (", k[1], ", ", k[2], ")")
        })
    string <- paste0("Robust GREG estimator (", string, ")")
    # multiplicative factor for Mallows type GM-estimator
    if (object$estimator$type == 1)
        bi <- bi * object$model$xwgt
    # bias-correction term
    zi <- bi * object$residuals
    est <- est + sum(zi)
    # variance estimate
    design <- object$design
    v <- survey::svyrecvar(zi, design$cluster, design$strata, design$fpc,
        postStrata = design$postStrata)
    res <- structure(list(characteristic = "total",
        estimator = list(string = string,
        type = type, psi = 0, psi_fun = "Huber", k = k), estimate = est,
        robust = list(robweights = bi / w), residuals = r,
        model = list(object$model, coef = object$estimate,
        auxiliary = auxiliary), design = design, call = match.call(),
        variance = v), class = "svystat_rob")
    if (keep_object)
        res$object <- object
    res
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
