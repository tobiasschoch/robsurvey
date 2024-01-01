\name{weighted_median}
\alias{weighted_median}
\title{Weighted Median}
\usage{
weighted_median(x, w, na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
    Weighted population median.
}
\details{
    Weighted sample median; see \code{\link{weighted_quantile}} for more
    information.
}
\value{
    Weighted estimate of the population median
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{weighted_quantile}}
}
\examples{
head(workplace)

weighted_median(workplace$employment, workplace$weight)
}
