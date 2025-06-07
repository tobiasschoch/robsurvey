# Within tolerance limits; observations that fall outside are declared outliers
within_tolerance <- function(x, w, method = c("quartile", "modified",
    "boxplot"), constants, lambda = 0.05, info = FALSE, boxplot_coef = 1.5)
{
    stopifnot(is.numeric(x), is.numeric(w), is.numeric(lambda))
    if (is.matrix(x) | is.data.frame(x))
        stop("Argument x must be a numeric vector\n", call. = FALSE)

    method <- match.arg(method, c("quartile", "modified", "boxplot"))
    if (method == "boxplot")
        constants <- c(1, 1)       # dummy argument; not used
    else if (missing(constants))
        stop("Argument constants must be defined\n", call. = FALSE)

    if (length(constants) != 2)
        stop("Cutoff must be a numeric vector of length 2\n",
            call. = FALSE)
    if (any(constants < 0))
        stop("Cutoff values must be nonnegative\n", call. = FALSE)

    if (lambda < 0 | lambda > 1)
        stop("Lambda must be in the interval [0, 1]\n", call. = FALSE)

    # trivial case
    if (length(x) == 1) {
        if (info)
            cat("[NA]\n")
        return(TRUE)
    }

    # regular cases
    q <- unname(weighted_quantile(x, w, probs = c(0.25, 0.5, 0.75), na.rm =
        TRUE))

    if (method == "quartile") {
        lb <- q[2] - constants[1] * (q[2] - q[1])
        ub <- q[2] + constants[2] * (q[3] - q[2])
    } else if (method == "modified"){
        ll <- lambda * abs(q[2])
        lb <- q[2] - constants[1] * max((q[2] - q[1]), ll)
        ub <- q[2] + constants[2] * max((q[3] - q[2]), ll)

    } else if (method == "boxplot") {
        iqr <- q[3] - q[1]
        lb <- q[1] - boxplot_coef * iqr
        ub <- q[3] + boxplot_coef * iqr
    }

    # is potential outlier
    out <- x < lb | x > ub
    out[is.na(x)] <- FALSE

    if (info) {
        min_x <- min(x, na.rm = TRUE)
        if (lb < min_x)
            lb <- "min."
        max_x <- max(x, na.rm = TRUE)
        if (ub > max_x)
            ub <- "max."
        cat(paste0("[", lb, ", ", ub, "]\n"))
    }

    !out
}
