\name{weighted_line}
\alias{weighted_line}
\title{Weighted Robust Line Fitting}
\usage{
weighted_line(x, y = NULL, w, na.rm = FALSE, iter = 1)
}
\arguments{
\item{x}{\code{[numeric vector]} explanatory variable.}
\item{y}{\code{[numeric vector]} response variable (default: \code{NULL}).}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
\item{iter}{\code{[integer]} number of iterations to use (default: \code{1}).}
}
\description{
\code{weighted_line} fits a robust line and allows weights.
}
\details{
\code{weighted_line} uses different quantiles for splitting the sample than
\code{stats::line()}.
}
\value{
intercept and slope of the fitted line
}
\examples{
data(cars)

# compute weighted line
weighted_line(cars$speed, cars$dist, w = rep(1, length(cars$speed)))
weighted_line(cars$speed, cars$dist, w = rep(1:10, each = 5))
}
\seealso{
\code{\link[stats]{line}}
}
