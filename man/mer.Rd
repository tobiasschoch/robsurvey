\name{mer}
\alias{mer}
\title{Minimum Estimated Risk (MER) M-Estimator}
\usage{
mer(object, verbose = TRUE, k = c(1, 50))
}
\arguments{
\item{object}{an M-estimate of the total or mean (object of class
    \code{svystat_rob}).}
\item{verbose}{\code{[logical]} indicating whether additional information
    is printed to the console (default: \code{TRUE}).}
\item{k}{\code{[numeric vector]} of size 2: the interval of robustness
    tuning constant within which the MER is supposed to lie.}
}
\value{
The tuning constant that minimizes the estimated mean square error of the
estimator
}
\description{
\code{mer} is an adaptive M-estimator of the weighted mean or the total. It
is defined as the estimator that minimizes the estimates mean square error
of the estimator under consideration.
}
\details{
MER-estimators are available for the methods \code{\link{svymean_huber}},
\code{\link{svytotal_huber}}, \code{\link{svymean_tukey}} and
\code{\link{svytotal_tukey}}
}
\examples{
library(survey)
data(losdata)
dn <- svydesign(ids = ~1, fpc = ~fpc, weights = ~weight, data = losdata)

# M-estimator of the total with tuning constant k = 8
m <- svymean_huber(~los, dn, type = "rwm", k = 8)

# mer-estimator
mer(m)
}
\references{
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators,
\emph{Surv. Methodol.} 21, pp. 79-87.
}

