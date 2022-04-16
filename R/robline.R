# quantiles as implemented in line() but with weights x and w sorted
# according to x
.line_quantile <- function(x, w, prob)
{
    n <- length(x)
    cumw <- n * cumsum(w) / sum(w)
    low <- min(which(cumw >= (n - 1) * prob))
    high <- max(which(cumw <= (n - 1) * prob))
    return((x[low + 1] + x[high + 1]) /2)
}
# check data for missing values
.check_line_data <- function(x, y, w, na.rm = FALSE)
{
    ok <- stats::complete.cases(x, y, w)
    if (sum(ok) < length(x)) {
        if (na.rm) {
	        x <- x[ok]
            y <- y[ok]
            w <- w[ok]
        } else {
	        stop("There are missing values in 'x', 'y' or 'w'.\n",
                call. = FALSE)
        }
    }
    list(x = x, y = y, w = w)
}
# Weighted line estimator
weighted_line <- function(x, y = NULL, w, na.rm = FALSE, iter = 1)
{
    if(missing(w))
        stop("Argument 'w' (weights) is missing, with no default.\n")

    if (inherits(x, "formula")) {
        mf <- stats::model.frame(x)
        y <- mf[, 1]
        x <- mf[, -1]
    }
    if (NCOL(x) > 1)
        stop("Argument 'x' contains more than 1 explanatory variables.\n",
            call. = FALSE)

    # check for missing values
    tmp <- .check_line_data(x, y, w, na.rm)
    x <- tmp$x; y <- tmp$y; w <- tmp$w; n <- length(x)

    # Sort data according to x
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]

    # standardise weights to n
    w <- n * w[ord] / sum(w)

    # Groups
    lim1 <- .line_quantile(x, w, 1 / 3)
    lim2 <- .line_quantile(x, w, 2 / 3)
    groups <- (x <= lim1) * 1 + (x > lim1 & x < lim2) * 2 + (x >= lim2) * 3

    # Medians for x
    wmedx <- c(weighted_median(x[groups == 1], w[groups == 1], na.rm),
	    weighted_median(x[groups == 3], w[groups == 3], na.rm))

    # polishing (affects only the slope)
    slope <- 0; r <- y; j <- 0
    while (j <= iter - 1) {
        wmedr <- c(weighted_median(r[groups == 1], w[groups == 1], na.rm),
	        weighted_median(r[groups == 3], w[groups == 3], na.rm))
        slope <- slope + (wmedr[2] - wmedr[1]) / (wmedx[2] - wmedx[1])
        r <- y - slope * x
        j <- j + 1
    }

    # intercept and predicted values
    intercept <- weighted_median(r, w, na.rm)
    yhat <- intercept + slope * x

    structure(list(call = match.call(), coefficients = c(intercept, slope),
        residuals = y - yhat, fitted.values = yhat), class = "tukeyline")
}
# Weighted median line estimator
weighted_median_line <- function(x, y = NULL, w, type = "slopes",
    na.rm = FALSE)
{
    if (missing(w))
        stop("Argument 'w' (weights) is missing, with no default.\n",
            call. = FALSE)

    if (inherits(x, "formula")) {
        dat <- grDevices::xy.coords(x)
        x <- dat$x
        y <- dat$y
    }
    if (NCOL(x) > 1)
        stop("Argument 'x' contains more than 1 explanatory variables.\n",
            call. = FALSE)

    # check for missing values
    tmp <- .check_line_data(x, y, w, na.rm)
    x <- tmp$x; y <- tmp$y; w <- tmp$w

    # univariate medians
    wmedx <- weighted_median(x, w, na.rm)
    wmedy <- weighted_median(y, w, na.rm)

    # type (match.arg is used for partial matching)
    slope <- switch(match.arg(type, c("slopes", "products")),
        "slopes" = {
            tmp <- (y - wmedy) / (x - wmedx)
            at <- is.finite(tmp)
            weighted_median(tmp[at], w[at], na.rm = TRUE)
        },
        "products" = weighted_median((y - wmedy) * (x - wmedx), w,
            na.rm = TRUE) / weighted_median((x - wmedx)^2, w, na.rm = TRUE),
        stop("Argument 'type' must be either 'slopes' or 'products'\n")
   )

    # residuals
    r <- y - slope * x

    # intercept
    intercept <- weighted_median(r, w, na.rm)
    yhat <- intercept + slope * x

    structure(list(call = sys.call(), coefficients = c(intercept, slope),
        residuals = y - yhat, fitted.values = yhat), class = "medline")
}
# Weighted median ratio estimator
weighted_median_ratio <- function(x, y = NULL, w, na.rm = FALSE)
{
    if (missing(w))
        stop("Argument 'w' (weights) is missing, with no default.\n",
            call. = FALSE)

    if (inherits(x, "formula")) {
        dat <- grDevices::xy.coords(x)
        y <- dat$y
        x <- dat$x
    }
    if (NCOL(x) > 1)
        stop("Argument 'x' contains more than 1 explanatory variables.\n",
            call. = FALSE)

    # check for missing values
    tmp <- .check_line_data(x, y, w, na.rm)
    x <- tmp$x; y <- tmp$y; w <- tmp$w

    # ratio
    ratio <- weighted_median(y / x, w, na.rm)

    # fitted.varlues and residuals
    yhat <- ratio * x
    r <- y - yhat

    structure(list(call = match.call(), coefficients = ratio,
        residuals = y - yhat, fitted.values = yhat), class = "medline")
}
# some S3 functions for objects of class 'medline'
print.medline <- function(x, ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    cat("\nCoefficients:\n")
    print(x$coefficients)
}
coef.medline <- function(object, ...)
{
    object$coefficients
}
residuals.medline <- function(object, ...)
{
    object$residuals
}
fitted.medline <- function(object, ...)
{
    object$fitted.values
}
