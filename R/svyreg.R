# Regression estimator of the mean (depends on pkg survey)
svyreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    res <- svyreg_huber(formula, design, var, k = Inf, na.rm)
    res$estimator$string <- "Survey regression estimator"
    res$estimator$k <- Inf
    res$call <- match.call()
    class(res) <- "svyreg_rob"
    res
}
