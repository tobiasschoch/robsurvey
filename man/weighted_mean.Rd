\name{weighted_mean}
\alias{weighted_mean}
\alias{weighted_total}
\title{Weighted Total and Mean (Horvitz-Thompson and Hajek Estimators)}
\usage{
weighted_mean(x, w, na.rm = FALSE)
weighted_total(x, w, na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values
        should be removed before the computation proceeds
        (default: \code{FALSE}).}
}
\description{
    Weighted total and mean (Horvitz-Thompson and Hajek estimators)
}
\details{
    \code{weighted_total} and \code{weighted_mean} compute, respectively,
    the Horvitz-Thompson estimator of the population total and the Hajek
    estimator of the population mean.
}
\value{
    Estimated population mean or total
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\examples{
head(workplace)

# Horvitz-Thompson estimator of the total
weighted_total(workplace$employment, workplace$weight)

# Hajek estimator of the mean
weighted_mean(workplace$employment, workplace$weight)
}
