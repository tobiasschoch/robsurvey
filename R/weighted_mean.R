# weighted mean (Hajek estimator)
weighted_mean <- function(x, w, na.rm = FALSE)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    else
        return(sum(dat$x * dat$w) / sum(dat$w))
}
# weighted total
weighted_total <- function(x, w, na.rm = FALSE)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    else
        return(sum(dat$w * dat$x))
}
