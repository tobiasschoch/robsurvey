\name{mer}
\alias{mer}
\title{Minimum Estimated Risk (MER) M-Estimator}
\usage{
mer(object, verbose = TRUE, max_k = 1000, optim_args = list())
}
\arguments{
\item{object}{an M-estimate of the total or mean (object of class
    \code{svystat_rob}).}
\item{verbose}{\code{[logical]} indicating whether additional information
    is printed to the console (default: \code{TRUE}).}
\item{max_k}{\code{[numeric vector]} defines the right boundary of the
    search interval (default: \code{max_k = 1000})}
\item{optim_args}{\code{[list]}: arguments passed on to \code{\link{optim}}.}
}
\value{
The tuning constant that minimizes the estimated mean square error of the
estimator
}
\description{
\code{mer} is an adaptive M-estimator of the weighted mean or total. It
is defined as the estimator that minimizes the estimated mean square error
of the estimator under consideration.
}
\details{
MER-estimators are available for the methods \code{\link{svymean_huber}},
\code{\link{svytotal_huber}}, \code{\link{svymean_tukey}} and
\code{\link{svytotal_tukey}}.
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
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators.
\emph{Survey Methodology} 21, pp. 79-87.
}
