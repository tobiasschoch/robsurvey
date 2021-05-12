# weighted median
weighted_median <- function(x, w, na.rm = FALSE)
{
    weighted_quantile(x, w, 0.5, na.rm)
}
