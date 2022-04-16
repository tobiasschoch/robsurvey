\name{svymean_reg}
\alias{svymean_reg}
\alias{svytotal_reg}
\title{Generalized Regression Estimator (GREG) of the Mean and Total}
\usage{
svymean_reg(object, auxiliary, k = Inf, psi = c("Huber", "Huberasym", "Tukey"),
    check.names = TRUE, na.rm = FALSE)
}
\arguments{
\item{object}{of class \code{[svyreg_rob]}, e.g., result of the Huber
    regression \eqn{M}-estimator \code{svyreg_huberM}.}
\item{auxiliary}{\code{[numeric]} mean or total of the auxiliary variables.}
\item{k}{\code{[numeric]} robustness tuning constant of the
    \code{psi}-function used in the bias-correction term in the GREG
    (default: \code{Inf}, which implies a non-robust GREG).}
\item{psi}{\code{[character]} \eqn{\psi}{\psi}-function.}
\item{check.names}{\code{[logical]} if \code{TRUE}, the names of
    \code{auxiliary} is checked against the independent variables of the
    fitted model \code{object} (default: \code{TRUE}.}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
(Robust) Generalized Regression Estimator (GREG) of the Mean and Total
}
\details{
A
}
\value{
Object of class \code{\link{svystat_rob}}
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Robust regression M-estimator with Huber psi-function
reg <- svyreg_huberM(payroll ~ employment, dn, k = 3)

# Robust GREG estimator of the mean; the population means of the auxiliary
# variables are from a register
m <- svymean_reg(reg, auxiliary = c(1, 11.03))
m

# Summarize
summary(m)

# Extract estimate
coef(m)

# Extract estimate of scale
scale(m)

# Extract estimated standard error
SE(m)
}
\seealso{
\code{\link{svyreg}}, \code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}}, and \code{\link{svyreg_tukeyGM}}
}
