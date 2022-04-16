\name{svyreg}
\alias{svyreg}
\title{Survey Regression Estimator}
\usage{
svyreg(formula, design, var = NULL, na.rm = FALSE)
}
\arguments{
\item{formula}{a \code{[formula]} object (i.e., symbolic description of the
    model)}
\item{design}{an object of class \code{survey.design}; see
    \code{\link[survey]{svydesign}}.}
\item{var}{\code{[numeric vector]} heteroscedastic variance (default:
    \code{NULL}, i.e., homoscedastic variance).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
\code{svyreg} is used to fit survey weighted linear models.
}
\details{
\code{svyreg} computes the regression coefficients by weighted least squares.

Models for \code{svyreg_rob} are specified symbolically. A typical model has
the form \code{response ~ terms} where \code{response} is the (numeric)
response vector and \code{terms} is a series of terms which specifies a
linear predictor for \code{response}; see \code{\link{formula}} and
\code{\link{lm}}.

A formula has an implied intercept term. To remove this use either
\code{y ~ x - 1} or \code{y ~ 0 + x}; see \code{\link{formula}} for more
details of allowed formulae.
}
\value{
Object of class \code{\link{svyreg_rob}}.
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Compute the regression estimate (weighted least squares)
m <- svyreg(payroll ~ employment, dn)

# Regression inference
summary(m)

# Extract the coefficients
coef(m)

# Extract variance/ covariance matrix
vcov(m)
}
\seealso{
\code{\link[=svyreg_rob]{summary}} for summaries.

The generic functions \code{\link[=svyreg_rob]{coef}},
\code{\link[=svyreg_rob]{residuals}}, \code{\link[=svyreg_rob]{fitted}}, and
\code{\link[=svyreg_rob]{vcov}}.

\code{\link[=svyreg_rob]{plot}} for regression diagnostic plot methods.

Robust estimating methods \code{\link{svyreg_huberM}},
\code{\link{svyreg_huberGM}}, \code{\link{svyreg_tukeyM}}, and
\code{\link{svyreg_tukeyGM}}.
}
