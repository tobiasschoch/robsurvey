\name{svyratio_huber}
\alias{svyratio_huber}
\alias{svyratio_tukey}
\title{Robust Survey Ratio M-Estimator}
\usage{
svyratio_huber(numerator, denominator, design, k, var = denominator,
               na.rm = FALSE, asym = FALSE, verbose = TRUE, ...)
svyratio_tukey(numerator, denominator, design, k, var = denominator,
               na.rm = FALSE, verbose = TRUE, ...)
}
\arguments{
    \item{numerator}{a one-sided \code{[formula]} object (i.e., symbolic
        description, e.g., \code{~payroll}).}
    \item{denominator}{a one-sided \code{[formula]} object (i.e., symbolic
        description, e.g., \code{~employment}).}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{var}{a \code{[formula]} object that defines the heteroscedastic
        variance (default: \code{numerator}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{asym}{\code{[logical]} toggle for asymmetric Huber psi-function
        (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
        maxit number of iterations, etc.).}
}
\description{
    \code{svyratio_huber} and \code{svyratio_tukey} compute the robust
    \eqn{M}-estimator of the ratio of two variables with, respectively,
    Huber and Tukey biweight (bisquare) psi-function.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    The functions \code{svyratio_huber} and \code{svyratio_tukey} are
    implemented as wrapper functions of the regression estimators
    \code{\link{svyreg_huberM}} and \code{\link{svyreg_tukeyM}}. See
    the help files of these functions (e.g., on how additional
    parameters can be passed via \code{\dots} or on the usage of the
    \code{var} argument).
}
\value{
    Object of class \code{svyreg.rob} and \code{ratio}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link[=svyreg_rob]{summary}}, \code{\link[=svyreg_rob]{coef}},
\code{\link[=svyreg_rob]{residuals}}, \code{\link[=svyreg_rob]{fitted}},
\code{\link[=svyreg_rob]{SE}} and \code{\link[=svyreg_rob]{vcov}}

\code{\link[=svyreg_rob]{plot}} for regression diagnostic plot methods

\code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}} and \code{\link{svyreg_tukeyGM}} for robust
regression estimators
}
\examples{
head(workplace)

library(survey)
# Survey design for stratified simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace, calibrate.formula = ~-1 + strat)
    } else {
        # legacy mode
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace)
    }

# Compute regression M-estimate with Huber psi-function
m <- svyratio_huber(~payroll, ~employment, dn, k = 8)

# Regression inference
summary(m)

# Extract the coefficients
coef(m)

# Extract estimated standard error
SE(m)

# Extract variance/ covariance matrix
vcov(m)

# Diagnostic plots (e.g., standardized residuals against fitted values)
plot(m, which = 1L)

# Plot of the robustness weights of the M-estimate against its residuals
plot(residuals(m), robweights(m))
}
