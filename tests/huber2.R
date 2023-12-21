library("robsurvey", quietly = TRUE)

# test on floating point equality
all_equal <- function(target, current, label,
    tolerance = sqrt(.Machine$double.eps), scale = NULL,
    check.attributes = FALSE)
{
    if (missing(label))
        stop("Argument 'label' is missing\n")
    res <- all.equal(target, current, tolerance, scale,
        check.attributes = check.attributes)
    if (is.character(res))
        cat(paste0(label, ": ", res, "\n"))
}

# test against MASS::hubers
if (requireNamespace("MASS", quietly = TRUE)) {
    library(MASS)

    # convergence tolerance
    TOLERANCE <- 1e-8

    # make a copy of the function MASS::hubers
    hubers_mod <- hubers
    # replace the mad by the (scaled) IQR as initial scale estimator
    body(hubers_mod)[[7]][[3]][[2]] <- substitute(s0 <- IQR(y, type = 2) *
        0.7413)

    # payroll data
    attach(workplace)
    all_equal(huber2(payroll, rep(1, length(payroll)), tol = TOLERANCE),
        hubers_mod(payroll, tol = 1e-8)$mu,
        label = "huber2: payroll data: test against MASS::hubers failed")

    # x11 data
    x11 <- c(1:10, 1e99); w11 <- rep(1, 11)
    all_equal(huber2(x11, w11, tol = TOLERANCE),
        hubers_mod(x11, tol = 1e-8)$mu,
        label = "huber2: x11 data: test against MASS::hubers failed")

    # x11 data (scaled)
    x11 <- c(1:10, 1e99) / 1e20; w11 <- rep(1, 11)
    all_equal(huber2(x11, w11, tol = TOLERANCE),
        hubers_mod(x11, tol = 1e-8)$mu,
        label = "huber2: x11 data (scale): test against MASS::hubers failed")
}
