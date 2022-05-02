# weighted interquartile range
weighted_IQR <- function(x, w, na.rm = FALSE, constant = 0.7413)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    qs <- weighted_quantile(dat$x, dat$w, probs = c(0.25, 0.75), na.rm)
    unname((qs[2] - qs[1]) * constant)
}
