\name{weighted_line}
\alias{weighted_line}
\alias{print.medline}
\alias{coef.medline}
\alias{residuals.medline}
\alias{fitted.medline}
\title{Weighted Robust Line Fitting}
\usage{
weighted_line(x, y = NULL, w, na.rm = FALSE, iter = 1)

\method{print}{medline}(x, ...)
\method{coef}{medline}(object, ...)
\method{residuals}{medline}(object, ...)
\method{fitted}{medline}(object, ...)
}
\arguments{
    \item{x}{\code{[numeric vector]} explanatory variable.}
    \item{y}{\code{[numeric vector]} response variable (default: \code{NULL}).}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{iter}{\code{[integer]} number of iterations to use (default:
        \code{1}).}
    \item{object}{object of class \code{medline}.}
    \item{\dots}{additional arguments passed to the method.}
}
\description{
    \code{weighted_line} fits a robust line and allows weights.
}
\details{
    \code{weighted_line} uses different quantiles for splitting the sample
    than \code{stats::line()}.
}
\value{
    intercept and slope of the fitted line
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link[stats]{line}}
}
\examples{
head(cars)

# compute weighted line
weighted_line(cars$speed, cars$dist, w = rep(1, length(cars$speed)))
m <- weighted_line(cars$speed, cars$dist, w = rep(1:10, each = 5))
m
coef(m)
residuals(m)
fitted(m)
}
