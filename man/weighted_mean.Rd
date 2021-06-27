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
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
Estimated population mean or total
}
\description{
Weighted total and mean (Horvitz-Thompson and Hajek estimators)
}
\details{
\code{weighted_total} and \code{weighted_mean} compute, respectively, the
Horvitz-Thompson estimator of the total and the Hajek estimator of the mean.
}
\examples{
data(workplace)

# Horvitz-Thompson estimator of the total
weighted_total(workplace$employment, workplace$weight)

# Hajek estimator of the mean
weighted_mean(workplace$employment, workplace$weight)
}
