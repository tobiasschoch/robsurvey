\name{weighted_median}
\alias{weighted_median}
\title{Weighted Median}
\usage{
weighted_median(x, w, na.rm = FALSE)
}
\arguments{
\item{x}{\code{[numeric vector]} observations.}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
Weighted estimate of the population median
}
\description{
\code{weighted_median} computes the weighted sample quantile
}
\details{
Weighted sample median; see \code{\link{weighted_quantile}} for more
information.
}
\examples{
data(workplace)

weighted_median(workplace$employment, workplace$weight)
}
\seealso{
\code{\link{weighted_quantile}}
}
