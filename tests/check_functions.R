# functions that are required in the testing suite

# utility function
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

# check functions for coef() and SE()
check <- function(ref, est, what, characteristic)
{
    all_equal(as.numeric(coef(ref)), as.numeric(coef(est)),
        paste0(what, ": ", characteristic, ": coef"))
    all_equal(as.numeric(SE(ref)), as.numeric(SE(est)),
        paste0(what,": ", characteristic, ": SE"))
}
