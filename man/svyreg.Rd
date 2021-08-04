\name{svyreg}
\alias{svyreg}
\title{Survey Regression Estimator}
\usage{
svyreg(formula, design, var = NULL, na.rm = FALSE)
}
\arguments{
\item{formula}{a \code{[formula]} object (i.e., symbolic description of the
    model)}
\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}
\item{var}{\code{[numeric vector]} heteroscedastic variance (default:
    \code{NULL}, i.e., homoscedastic variance).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
Object of class \code{\link[=svyreg_rob]{svyreg_rob}}
}
\description{
Weighted survey regression estimator
}
\details{
\code{\link[=svyreg_rob]{residuals}}
\code{\link[=svyreg_rob]{coef}}
\code{\link[=svyreg_rob]{summary}}
\code{\link[=svyreg_rob]{vcov}}
\code{\link[=svyreg_rob]{fitted}}
\code{\link[=svyreg_rob]{plot}}
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)
# Compute the regression estimate
m <- svyreg(payroll ~ employment, dn)
# Regression inference
summary(m)
}
\seealso{
\code{\link{svyreg_huber}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukey}}, and \code{\link{svyreg_tukeyGM}}
}
