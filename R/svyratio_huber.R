# robust ratio estimator
svyratio_huber <- function(numerator, denominator, design, k,
                           var = denominator, na.rm = FALSE, asym = FALSE,
                           verbose = TRUE, ...)
{
    stopifnot(inherits(numerator, "formula"), inherits(denominator, "formula"))
    res <- svyreg_huberM(as.formula(paste0(numerator[[2]], "~ -1 +",
        denominator[[2]])), design, k, var, na.rm, asym, verbose, ...)
    names(res$estimate) <- paste0(numerator[[2]], "/", denominator[[2]])
    res$call <- match.call()
    res$characteristic <- "ratio"
    res$estimator$string <- paste0("Survey ratio M-estimator (Huber psi, k = ",
        k,")")
    class(res) <- c(class(res), "ratio")
    res
}
svyratio_tukey <- function(numerator, denominator, design, k,
    var = denominator, na.rm = FALSE, verbose = TRUE, ...)
{
    stopifnot(inherits(numerator, "formula"), inherits(denominator, "formula"))
    res <- svyreg_tukeyM(as.formula(paste0(numerator[[2]], "~ -1 +",
        denominator[[2]])), design, k, var, na.rm, verbose, ...)
    names(res$estimate) <- paste0(numerator[[2]], "/", denominator[[2]])
    res$call <- match.call()
    res$characteristic <- "ratio"
    res$estimator$string <- paste0("Survey ratio M-estimator (Tukey psi, k = ",
        k,")")
    class(res) <- c(class(res), "ratio")
    res
}
