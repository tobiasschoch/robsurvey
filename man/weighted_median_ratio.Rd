\name{weighted_median_ratio}
\alias{weighted_median_ratio}
\title{Weighted Robust Ratio Estimator Based on Median}
\usage{
weighted_median_ratio(x, y = NULL, w, na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} explanatory variable.}
    \item{y}{\code{[numeric vector]} response variable (default: \code{NULL}).}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
    A weighted median of the ratios y/x determines the slope of a regression
    through the origin.
}
\value{
    A vector with two components: intercept and slope
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link[stats]{line}}, \code{\link{weighted_line}} and
\code{\link{weighted_median_line}}
}
\examples{
x <- c(1,2,4,5)
y <- c(1,0,5,2)
m <- weighted_median_ratio(y ~ x, w = rep(1, length(y)))
m
coef(m)
fitted(m)
residuals(m)
}
