# weighted median absolute deviations about the weighted median
weighted_mad <- function(x, w, na.rm = FALSE, constant = 1.482602)
{
    dat <- .check_data_weights(x, w, na.rm)
    if (is.null(dat))
        return(NA)
    med <- weighted_quantile(dat$x, dat$w, probs = 0.5, na.rm)
    mad <- weighted_quantile(abs(dat$x - med), dat$w, probs = 0.5, na.rm)
    unname(mad * constant)
}
