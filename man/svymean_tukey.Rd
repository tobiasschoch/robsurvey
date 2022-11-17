\name{svymean_tukey}
\alias{svymean_tukey}
\alias{svytotal_tukey}
\title{Weighted Tukey Mean and Total -- Robust Horvitz-Thompson
Estimator}
\usage{
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
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method (e.g.,
        \code{maxit}: maxit number of iterations, etc.; see
        \code{\link{svyreg_control}}).}
}
\description{
    Weighted Tukey biweight (or bisquare) \emph{M}-estimator of the population
    mean and total (robust Horvitz-Thompson estimator)
}
\details{
    Package \pkg{survey} must be loaded in order to use the functions.

    \describe{
        \item{Methods/ types}{\code{type = "rht"} or \code{type = "rwm"}; see
            \code{\link{weighted_mean_tukey}} for more details.}
        \item{Variance estimation.}{Taylor linearization (residual variance
            estimator).}
        \item{Utility functions}{\code{\link[=svystat_rob]{summary}},
            \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
            \code{\link[=svystat_rob]{vcov}},
            \code{\link[=svystat_rob]{residuals}},
            \code{\link[=svystat_rob]{fitted}},
            \code{\link[=svystat_rob]{robweights}}.}
       \item{Bare-bone functions}{See \code{\link{weighted_mean_tukey}} and
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

\code{\link{svymean_huber}} and \code{\link{svytotal_huber}}
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Robust Horvitz-Thompson M-estimator of the population total
svytotal_tukey(~employment, dn, k = 9, type = "rht")

# Robust weighted M-estimator of the population mean
m <- svymean_tukey(~employment, dn, k = 12, type = "rwm")

# Summarize
summary(m)

# Extract estimate
coef(m)

# Extract estimate of scale
scale(m)

# Extract estimated standard error
SE(m)
}
