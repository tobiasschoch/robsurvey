\name{weighted_IQR}
\alias{weighted_IQR}
\title{Weighted Interquartile Range (IQR)}
\usage{
weighted_IQR(x, w, na.rm = FALSE, constant = 0.7413)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values
        should be removed before the computation proceeds
        (default: \code{FALSE}).}
    \item{constant}{\code{[double]} constant scaling factor to make the
        weighted IQR a consistent estimator of the scale
        (default: \code{0.7413}).}
}
\description{
    Weighted (normalized) interquartile range
}
\details{
    By default, the weighted IQR is normalized to be an unbiased estimate of
    scale at the Gaussian core model. If normalization is not wanted, put
    \code{constant = 1}.
}
\value{
    Weighted IQR
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\examples{
head(workplace)

# normalized weighted IQR (default)
weighted_IQR(workplace$employment, workplace$weight)

# weighted IQR (without normalization)
weighted_IQR(workplace$employment, workplace$weight, constant = 1)
}
