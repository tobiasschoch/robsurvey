\name{svyreg_huberM}
\alias{svyreg_huberM}
\alias{svyreg_huberGM}
\title{Huber Robust Survey Regression M- and GM-Estimator}
\usage{
svyreg_huberM(formula, design, k, var = NULL, na.rm = FALSE, asym = FALSE,
              verbose = TRUE, ...)
svyreg_huberGM(formula, design, k, type = c("Mallows", "Schweppe"),
               xwgt, var = NULL, na.rm = FALSE, asym = FALSE, verbose = TRUE,
               ...)
}
\arguments{
    \item{formula}{a \code{[formula]} object (i.e., symbolic description
        of the model)}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{var}{a one-sided \code{[formula]} object or variable name
        (\code{[character]}) that defines the heteroscedastic variance
        or \code{[NULL]} indicating homoscedastic variance (default:
        \code{NULL}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{asym}{\code{[logical]} toggle for asymmetric Huber psi-function
        (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{type}{\code{[character]} \code{"Mallows"} or \code{"Schweppe"}.}
    \item{xwgt}{\code{[numerical vector]} or \code{[NULL]} of weights in
        the design space (default: \code{NULL}); \code{xwgt} is only
        relevant for \code{type = "Mallows"} or \code{type = "Schweppe"}.}
    \item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
        maxit number of iterations, etc.).}
}
\description{
    \code{svyreg_huberM} and \code{svyreg_huberGM} compute, respectively,
    a survey weighted \emph{M}- and \emph{GM}-estimator of regression using
    the Huber psi-function.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    \code{svyreg_huberM} and \code{svyreg_huberGM} compute, respectively,
    \emph{M}- and \emph{GM}-estimates of regression by iteratively re-weighted
    least squares (IRWLS).  The estimate of regression scale is (by default)
    computed as the (normalized) weighted median of absolute deviations
    from the weighted median (MAD; see \code{\link{weighted_mad}}) for
    each IRWLS iteration. If the weighted MAD is zero (or nearly so),
    the scale is computed as the (normalized) weighted interquartile
    range (IQR).
    \describe{
        \item{M-estimator}{
            The regression M-estimator \code{svyreg_huberM} is robust against
            residual outliers (granted that the tuning constant \code{k} is
            chosen appropriately).
        }
        \item{GM-estimator}{
            Function \code{svyreg_huberGM} implements the Mallows and Schweppe
            regression GM-estimator (see argument \code{type}).
            The regression GM-estimators are robust against residual outliers
            \emph{and} outliers in the model's design space (leverage
            observations; see argument \code{xwgt}).
        }
        \item{Numerical optimization}{
            See \code{\link{svyreg_control}}.
        }
        \item{Models}{
            Models for \code{svyreg_rob} are specified symbolically. A typical
            model has the form \code{response ~ terms}, where \code{response}
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
    \code{tol = 1e-05}. If the algorithm fails to converge, you may
    consider changing the default values; see \code{\link{svyreg_control}}.
}
\value{
    Object of class \code{svyreg.rob}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link[=svyreg_rob]{summary}}, \code{\link[=svyreg_rob]{coef}},
\code{\link[=svyreg_rob]{residuals}}, \code{\link[=svyreg_rob]{fitted}},
\code{\link[=svyreg_rob]{SE}} and \code{\link[=svyreg_rob]{vcov}}

\code{\link[=svyreg_rob]{plot}} for regression diagnostic plot methods

Other robust estimating methods \code{\link{svyreg_tukeyM}} and
\code{\link{svyreg_tukeyGM}}
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
m <- svyreg_huberM(payroll ~ employment, dn, k = 8)

# Regression inference
summary(m)

# Extract the coefficients
coef(m)

# Extract variance/ covariance matrix
vcov(m)

# Diagnostic plots (e.g., standardized residuals against fitted values)
plot(m, which = 1L)

# Plot of the robustness weights of the M-estimate against its residuals
plot(residuals(m), robweights(m))
}
