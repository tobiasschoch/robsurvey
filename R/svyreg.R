svyreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
    res <- svyreg_huber(formula, design, var, k = Inf, na.rm)
    res$estimator <- "Survey regression estimator"
    res$call <- match.call()
    res$robust <- NULL
    res$optim <- NULL
    class(res) <- "svyreg_rob"
    res
}
