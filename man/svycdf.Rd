\name{svycdf_reg}
\alias{svycdf_reg}
\title{Regression Estimator of the Empirical Cumulative Distribution Function}
\usage{
svycdf_reg(object, at, non_sampled, type = "cd", ...)
}
\arguments{
    \item{object}{an object of class \code{[svyreg_rob]}, e.g., result of
        the (weighted) least squares estimator \code{svyreg} or
        the Huber regression \eqn{M}-estimator \code{svyreg_huberM}.}
    \item{at}{\code{[numeric]} vector of values at which the empirical
        cumulative distribution function (ECDF) should be evaluated.}
    \item{non_sampled}{\code{data.frame} with data of the \bold{non-sampled}
        part of the population.}
    \item{type}{\code{[character]} type of predictor; see Details Section.}
    \item{\dots}{additional arguments (currently not used).}
}
\description{
Compute an empirical cumulative distribution function (ECDF) using a
model-based or model-assisted estimation approach.
}
\details{
    First, we must fit the regression model using one of the methods:
    \code{\link{svyreg}}, \code{\link{svyreg_huberM}},
    \code{\link{svyreg_huberGM}}, \code{\link{svyreg_tukeyM}} or
    \code{\link{svyreg_tukeyGM}}. Second, we can compute an estimate of the
    ECDF using one of the following estimator:
    \itemize{
        \item Chambers and Dunstan (1986)
    }
}
\value{
    An estimate of the ECDF.
}
\references{
Chambers, R. L. and Dunstan, R. (1986). Estimating Distribution Functions from
    Survey Data. \emph{Biometrika}, \bold{73}, 597--604.
    \doi{10.1093/biomet/73.3.597}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svyreg}}, \code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}} and \code{\link{svyreg_tukeyGM}} for robust
regression \eqn{M}- and \eqn{GM}-estimators
}
\examples{
head(workplace)

library("survey")

# We define a design that mimics a simple random sample without replacement
workplace$N <- 86514
dn <- svydesign(~1, fpc = ~N, data = workplace)

# Fitted regression model with heteroscedastic variance structure
m <- svyreg(payroll ~ employment, var = ~employment, dn)

# Estimate of the ECDF at 1e6 with the Chambers-Dunstan estimator
svycdf_reg(m, 1e6, workplace)
}
