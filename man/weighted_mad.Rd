\name{weighted_mad}
\alias{weighted_mad}
\title{Weighted Median Absolute Deviation from the Median (MAD)}
\usage{
weighted_mad(x, w, na.rm = FALSE, constant = 1.482602)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values
        should be removed before the computation proceeds
        (default: \code{FALSE}).}
    \item{constant}{\code{[double]} constant scaling factor to make the
        MAD a consistent estimator of the scale (default: \code{1.4826}).}
    }
\description{
    Weighted median of the absolute deviations from the weighted median
}
\details{
    The weighted MAD is computed as the (normalized) weighted median of the
    absolute deviation from the weighted median; see
    \code{\link{weighted_median}}.  The weighted MAD is normalized to be
    an unbiased estimate of scale at the Gaussian core model. If
    normalization is not wanted, put \code{constant = 1}.
}
\value{
    Weighted median absolute deviation from the (weighted) median
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\examples{
head(workplace)

# normalized weighted MAD (default)
weighted_mad(workplace$employment, workplace$weight)

# weighted MAD (without normalization)
weighted_mad(workplace$employment, workplace$weight, constant = 1)
}
