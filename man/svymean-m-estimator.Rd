\name{svymean-m-estimator}
\alias{svymean_huber}
\alias{svytotal_huber}
\alias{svymean_tukey}
\alias{svytotal_tukey}
\title{Weighted Huber and Tukey Mean and Total (M-Estimator) -- Robust
    Horvitz-Thompson Estimator}
\usage{
svymean_huber(x, design, k, type = "rwm", asym = FALSE, na.rm = FALSE,
              verbose = TRUE, ...)
svytotal_huber(x, design, k, type = "rwm", asym = FALSE, na.rm = FALSE,
               verbose = TRUE, ...)
svymean_tukey(x, design, k, type = "rwm", na.rm = FALSE, verbose = TRUE, ...)
svytotal_tukey(x, design, k, type = "rwm", na.rm = FALSE, verbose = TRUE, ...)
}
\arguments{
    \item{x}{a one-sided \code{[formula]}, e.g., \code{~myVariable}.}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{type}{\code{[character]} type of method: \code{"rwm"} or
        \code{"rht"}.}
    \item{asym}{\code{[logical]} if \code{TRUE}, an asymmetric Huber
        psi-function is used (default: \code{FALSE}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
        maxit number of iterations, etc.; see \code{\link{svyreg_control}}).}
}
\description{
    Weighted Huber and Tukey \emph{M}-estimator of the population mean and
    total (robust Horvitz-Thompson estimator)
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    \describe{
        \item{Methods/ types}{\code{type = "rht"} or \code{type = "rwm"}; see
            \code{\link{weighted_mean_huber}} or
            \code{\link{weighted_mean_tukey}} for more details.}
        \item{Variance estimation.}{Taylor linearization (residual variance
            estimator).}
        \item{Utility functions}{\code{\link[=svystat_rob]{summary}},
            \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
            \code{\link[=svystat_rob]{vcov}},
            \code{\link[=svystat_rob]{residuals}},
            \code{\link[=svystat_rob]{fitted}},
            \code{\link[=svystat_rob]{robweights}}.}
        \item{Bare-bone functions}{See \code{\link{weighted_mean_huber}}
           \code{\link{weighted_mean_tukey}},
           \code{\link{weighted_total_huber}}, and
           \code{\link{weighted_total_tukey}}.}
    }
}
\section{Failure of convergence}{
    By default, the method assumes a maximum number of \code{maxit = 100}
    iterations and a numerical tolerance criterion to stop the iterations of
    \code{tol = 1e-05}. If the algorithm fails to converge, you may
    consider changing the default values; see \code{\link{svyreg_control}}.
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

# Robust Horvitz-Thompson M-estimator of the population total
svytotal_huber(~employment, dn, k = 9, type = "rht")

# Robust weighted M-estimator of the population mean
m <- svymean_huber(~employment, dn, k = 12, type = "rwm")

# Summary statistic
summary(m)

# Plot of the robustness weights of the M-estimate against its residuals
plot(residuals(m), robweights(m))

# Extract estimate
coef(m)

# Extract estimate of scale
scale(m)

# Extract estimated standard error
SE(m)
}
