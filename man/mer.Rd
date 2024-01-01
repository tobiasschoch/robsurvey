\name{mer}
\alias{mer}
\title{Minimum Estimated Risk (MER) M-Estimator}
\usage{
mer(object, verbose = TRUE, max_k = 10, init = 1, method = "Brent",
    optim_args = list())
}
\arguments{
    \item{object}{an object of class \code{svystat_rob}.}
    \item{verbose}{\code{[logical]} indicating whether additional
        information is printed to the console (default: \code{TRUE}).}
    \item{init}{\code{[numeric]} determines the left boundary value of the
        search interval and the initial value of the search; we must have
        \code{init < max_k}.}
    \item{method}{\code{[character]} the method of \code{\link{optim}} to
        be used.}
    \item{max_k}{\code{[numeric vector]} defines the right boundary value of
        the search interval (default: \code{max_k = 1000})}
    \item{optim_args}{\code{[list]}: arguments passed on to
        \code{\link{optim}}.}
}
\description{
    \code{mer} is an adaptive M-estimator of the weighted mean or total. It
    is defined as the estimator that minimizes the estimated mean square error,
    \code{\link{mse}}, of the estimator under consideration.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    MER-estimators are available for the methods \code{\link{svymean_huber}},
    \code{\link{svytotal_huber}}, \code{\link{svymean_tukey}} and
    \code{\link{svytotal_tukey}}.
}
\value{
Object of class \code{\link{svystat_rob}}
}
\references{
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators.
    \emph{Survey Methodology} \bold{21}, 79--87.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\examples{
head(losdata)

library(survey)
# Survey design for simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~1, fpc = ~fpc, weights = ~weight, data = losdata,
                  calibrate.formula = ~1)
    } else {
        # legacy mode
        svydesign(ids = ~1, fpc = ~fpc, weights = ~weight, data = losdata)
    }

# M-estimator of the total with tuning constant k = 8
m <- svymean_huber(~los, dn, type = "rhj", k = 8)

# MER estimator
mer(m)
}
