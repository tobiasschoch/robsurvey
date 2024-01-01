\name{svymean_trimmed}
\alias{svymean_trimmed}
\alias{svytotal_trimmed}
\title{Weighted Trimmed Mean and Total}
\usage{
svymean_trimmed(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
svytotal_trimmed(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
}
\arguments{
\item{x}{a one-sided \code{[formula]}, e.g., \code{~myVariable}.}
\item{design}{an object of class \code{survey.design}; see
    \code{\link[survey]{svydesign}}.}
\item{LB}{\code{[double]} lower bound of trimming such that
    \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
\item{UB}{\code{[double]} upper bound of trimming such that
    \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
    Weighted trimmed population mean and total.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    \describe{
    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote
        the estimated trimmed population mean; then, the estimated trimmed
        total is given by \eqn{\hat{N} \mu}{Nhat \mu} with
        \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where
        summation is over all observations in the sample.}
    \item{Trimming.}{The methods trims the \code{LB}\eqn{~\cdot 100\%}
        of the smallest observations and the (1 - \code{UB})\eqn{~\cdot 100\%}
        of the largest observations from the data.}
    \item{Variance estimation.}{Large-sample approximation based on the
        influence function; see Huber and Ronchetti (2009, Chap. 3.3) and
        Shao (1994).}
    \item{Utility functions.}{\code{\link[=svystat_rob]{summary}},
        \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
        \code{\link[=svystat_rob]{vcov}}, \code{\link[=svystat_rob]{residuals}},
        \code{\link[=svystat_rob]{fitted}},
        \code{\link[=svystat_rob]{robweights}}.}
    \item{Bare-bone functions.}{See \code{\link{weighted_mean_trimmed}} and
        \code{\link{weighted_total_trimmed}}.}
}
}
\value{
    Object of class \code{\link{svystat_rob}}
}
\references{
Huber, P. J. and Ronchetti, E. (2009). \emph{Robust Statistics},
    New York: John Wiley and Sons, 2nd edition.
    \doi{10.1002/9780470434697}

Shao, J. (1994). L-Statistics in Complex Survey Problems.
    \emph{The Annals of Statistics} \bold{22}, 976--967.
    \doi{10.1214/aos/1176325505}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{weighted_mean_trimmed}} and \code{\link{weighted_total_trimmed}}
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

# Estimated trimmed population total (5\% symmetric trimming)
svytotal_trimmed(~employment, dn, LB = 0.05, UB = 0.95)

# Estimated trimmed population mean (5\% trimming at the top of the distr.)
svymean_trimmed(~employment, dn, UB = 0.95)
}
