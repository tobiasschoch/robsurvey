\name{svyreg_tukey}
\alias{svyreg_tukey}
\alias{svyreg_tukeyM}
\alias{svyreg_tukeyGM}
\title{Tukey Biweight Robust Survey Regression M- and GM-Estimator}
\usage{
svyreg_tukeyM(formula, design, k, var = NULL, na.rm = FALSE, verbose = TRUE,
    ...)
svyreg_tukeyGM(formula, design, k, type = c("Mallows", "Schweppe"),
    xwgt, var = NULL, na.rm = FALSE, verbose = TRUE, ...)

# deprecated
svyreg_tukey(formula, design, k, var = NULL, na.rm = FALSE, verbose = TRUE,
    ...)
}
\arguments{
\item{formula}{a \code{[formula]} object (i.e., symbolic description of the
    model)}
\item{design}{an object of class \code{survey.design}; see
    \code{\link[survey]{svydesign}}.}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
\item{var}{\code{[numeric vector]} heteroscedastic variance (default:
    \code{NULL}, i.e., homoscedastic variance).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
\item{verbose}{\code{[logical]} indicating whether additional information
		is printed to the console (default: \code{TRUE}).}
\item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
    maxit number of iterations, etc.).}
\item{type}{\code{[character]} \code{"Mallows"} or \code{"Schweppe"}.}
\item{xwgt}{\code{[numerical vector]} of weights in the design space
    (default: \code{NULL}); \code{xwgt} is only relevant for
    \code{type = "Mallows"} or \code{type = "Schweppe"}.}
}
\description{
\code{svyreg_tukeyM} and \code{svyreg_tukeyGM} compute, respectively,
a survey weighted M- and GM-estimator of regression using the biweight
Tukey psi-function. Use \code{svyreg_tukeyM} not \code{svyreg_tukey}
because the latter is deprecated but is kept for computability reasons.
}
\details{
\code{svyreg_tukeyM} and \code{svyreg_tukeyGM} compute, respectively, M- and
GM-estimates of regression by iteratively re-weighted least squares (IRWLS).
The estimate of regression scale is (by default) computed by the (normalized)
weighted median of absolute deviations from the weighted median (MAD; see
\code{\link{weighted_mad}}) for each IRWLS iteration. If the weighted MAD
is zero (or nearly so), the scale is computed by the (normalized) weighted
interquartile range (IQR).
    \describe{
        \item{M-estimator}{
            The regression M-estimator \code{svyreg_tukeyM} is robust against
            residual outliers (granted that the tuning constant \code{k} is
            chosen appropriately).
        }
        \item{GM-estimator}{
            \code{svyreg_huberGM} implements the Mallows and Schweppe
            regression GM-estimator (see argument \code{type}).
            The regression GM-estimator are robust against residual outliers
            \emph{and} outliers in the model's design space (leverage
            observations; see argument \code{xwgt}).
        }
        \item{Numerical optimization}{
            See \code{\link{svyreg_control}}.
        }
        \item{Models}{
            Models for \code{svyreg_rob} are specified symbolically. A typical
            model has the form \code{response ~ terms} where \code{response}
            is the (numeric) response vector and \code{terms} is a series of
            terms which specifies a linear predictor for \code{response}; see
            \code{\link{formula}} and \code{\link{lm}}.

            A formula has an implied intercept term. To remove this use
            either \code{y ~ x - 1} or \code{y ~ 0 + x}; see
            \code{\link{formula}} for more details of allowed formulae.
        }
    }
}
\section{Failure of convergence}{
By default, the method assumes a maximum number of \code{maxit = 100}
iterations and a numerical tolerance criterion to stop the iterations of
\code{tol = 1e-05}. You can run the code with specifications other than
the default values by specifying the arguments \code{maxit} and/or
\code{tol} in the function call; see also \code{\link{svyreg_control}}.
}
\value{
Object of class \code{svyreg.rob}
}
\seealso{
\code{\link[=svyreg_rob]{summary}} for summaries.

The generic functions \code{\link[=svyreg_rob]{coef}},
\code{\link[=svyreg_rob]{residuals}}, \code{\link[=svyreg_rob]{fitted}}, and
\code{\link[=svyreg_rob]{vcov}}.

\code{\link[=svyreg_rob]{plot}} for regression diagnostic plot methods.

Other robust estimating methods \code{\link{svyreg_huberM}} and
\code{\link{svyreg_huberGM}}.
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Compute regression M-estimate with Tukey bisquare psi-function
m <- svyreg_tukeyM(payroll ~ employment, dn, k = 8)

# Regression inference
summary(m)

# Extract the coefficients
coef(m)

# Extract variance/ covariance matrix
vcov(m)
}
